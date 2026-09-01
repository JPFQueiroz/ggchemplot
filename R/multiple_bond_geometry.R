# ================================================================
# MULTIPLE-BOND GEOMETRY ENGINE
# ================================================================

# Small numerical tolerance
.mb_eps <- 1e-8


# ------------------------------------------------
# Basic 2-D geometry
# ------------------------------------------------

.mb_cross <- function(ax, ay, bx, by) {
  ax * by - ay * bx
}


.mb_side_of_bond <- function(x1, y1, x2, y2, x, y) {
  .mb_cross(
    x2 - x1,
    y2 - y1,
    x - x1,
    y - y1
  )
}


.mb_distance <- function(x1, y1, x2, y2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2)
}


# ------------------------------------------------
# Direct heavy-atom neighbours
#
# Uses the existing bonds/atoms structure.
# No modification of either object.
# ------------------------------------------------

.mb_neighbors <- function(atom_id, bonds, atoms,
                          exclude = NULL,
                          heavy_only = TRUE) {

  nb <- c(
    bonds$to[bonds$from == atom_id],
    bonds$from[bonds$to == atom_id]
  )

  nb <- unique(nb)

  if (!is.null(exclude)) {
    nb <- nb[nb != exclude]
  }

  if (heavy_only) {
    nb <- nb[
      nb %in% atoms$atom_id[
        atoms$symbol != "H"
      ]
    ]
  }

  nb
}


# ------------------------------------------------
# Find one alternate path between two atoms while
# temporarily removing the bond between them.
#
# If a path exists, the bond is part of a ring.
# ------------------------------------------------

.mb_find_alternate_path <- function(start,
                                    goal,
                                    bonds) {

  # adjacency from the existing bond table
  adjacency <- split(
    c(bonds$to, bonds$from),
    c(bonds$from, bonds$to)
  )

  # Remove the target edge
  adjacency[[as.character(start)]] <-
    setdiff(
      adjacency[[as.character(start)]],
      goal
    )

  adjacency[[as.character(goal)]] <-
    setdiff(
      adjacency[[as.character(goal)]],
      start
    )

  queue <- list(
    list(
      node = start,
      path = start
    )
  )

  visited <- start

  while (length(queue) > 0) {

    current <- queue[[1]]
    queue <- queue[-1]

    node <- current$node
    path <- current$path

    neighbours <- adjacency[[as.character(node)]]

    if (is.null(neighbours)) {
      next
    }

    for (nb in neighbours) {

      if (nb %in% visited) {
        next
      }

      new_path <- c(path, nb)

      if (nb == goal) {
        return(new_path)
      }

      visited <- c(visited, nb)

      queue[[length(queue) + 1]] <- list(
        node = nb,
        path = new_path
      )
    }
  }

  NULL
}


# ------------------------------------------------
# Determine whether a bond is in a ring and, if so,
# return the alternate path.
# ------------------------------------------------

.mb_ring_info <- function(from, to, bonds) {

  path <- .mb_find_alternate_path(
    start = from,
    goal = to,
    bonds = bonds
  )

  if (is.null(path)) {
    return(list(
      in_ring = FALSE,
      path = NULL,
      ring_atoms = integer(0)
    ))
  }

  list(
    in_ring = TRUE,
    path = path,
    ring_atoms = unique(path)
  )
}


# ------------------------------------------------
# Ring centroid
#
# For an endocyclic multiple bond, this gives us
# the direction toward the ring interior.
# ------------------------------------------------

.mb_ring_centroid <- function(ring_atoms, atoms) {

  pts <- atoms[
    atoms$atom_id %in% ring_atoms,
    c("x", "y")
  ]

  if (nrow(pts) == 0) {
    return(c(NA_real_, NA_real_))
  }

  c(
    mean(pts$x),
    mean(pts$y)
  )
}


# ------------------------------------------------
# Count multiple bonds belonging to a ring path.
# ------------------------------------------------

.mb_ring_multiple_bonds <- function(ring_atoms, bonds) {

  if (length(ring_atoms) < 2) {
    return(0L)
  }

  ring_atoms <- unique(ring_atoms)

  sum(
    bonds$order > 1 &
      bonds$from %in% ring_atoms &
      bonds$to %in% ring_atoms
  )
}


# ------------------------------------------------
# Obtain local substituents at each end.
#
# A substituent is a heavy-atom neighbour other
# than the other atom of the multiple bond.
# ------------------------------------------------

.mb_substituents <- function(from,
                             to,
                             bonds,
                             atoms) {

  s_from <- .mb_neighbors(
    atom_id = from,
    bonds = bonds,
    atoms = atoms,
    exclude = to,
    heavy_only = TRUE
  )

  s_to <- .mb_neighbors(
    atom_id = to,
    bonds = bonds,
    atoms = atoms,
    exclude = from,
    heavy_only = TRUE
  )

  list(
    from = s_from,
    to = s_to
  )
}


# ------------------------------------------------
# Determine the geometric side of a substituent.
# ------------------------------------------------

.mb_substituent_sides <- function(from,
                                  to,
                                  substituents,
                                  atoms) {

  a <- atoms[atoms$atom_id == from, ]
  b <- atoms[atoms$atom_id == to, ]

  if (nrow(a) == 0 || nrow(b) == 0) {
    return(numeric(0))
  }

  s <- atoms[
    atoms$atom_id %in% substituents,
    ,
    drop = FALSE
  ]

  if (nrow(s) == 0) {
    return(numeric(0))
  }

  .mb_side_of_bond(
    a$x,
    a$y,
    b$x,
    b$y,
    s$x,
    s$y
  )
}


# ------------------------------------------------
# Choose a side from a point.
#
# Returns +1 or -1.
# ------------------------------------------------

.mb_side_toward_point <- function(x1, y1,
                                  x2, y2,
                                  px, py) {

  value <- .mb_side_of_bond(
    x1, y1,
    x2, y2,
    px, py
  )

  if (value >= 0) {
    1L
  } else {
    -1L
  }
}


# ------------------------------------------------
# Determine whether substituents on opposite ends
# are trans.
#
# For one substituent on each end:
#
#   same side  -> cis-like
#   opposite   -> trans-like
# ------------------------------------------------

.mb_is_trans <- function(from,
                         to,
                         s_from,
                         s_to,
                         atoms) {

  if (length(s_from) != 1 ||
      length(s_to) != 1) {
    return(FALSE)
  }

  sf <- atoms[
    atoms$atom_id == s_from,
    ,
    drop = FALSE
  ]

  st <- atoms[
    atoms$atom_id == s_to,
    ,
    drop = FALSE
  ]

  a <- atoms[
    atoms$atom_id == from,
    ,
    drop = FALSE
  ]

  b <- atoms[
    atoms$atom_id == to,
    ,
    drop = FALSE
  ]

  if (nrow(sf) == 0 ||
      nrow(st) == 0 ||
      nrow(a) == 0 ||
      nrow(b) == 0) {
    return(FALSE)
  }

  side_f <- .mb_side_of_bond(
    a$x, a$y,
    b$x, b$y,
    sf$x, sf$y
  )

  side_t <- .mb_side_of_bond(
    a$x, a$y,
    b$x, b$y,
    st$x, st$y
  )

  side_f * side_t < 0
}


# ------------------------------------------------
# Decide the side of the secondary segment.
#
# Implements the local geometry corresponding to
# GR-1.10.1 through GR-1.10.4.
# ------------------------------------------------

.mb_choose_double_bond_side <- function(from,
                                        to,
                                        bonds,
                                        atoms) {

  a <- atoms[
    atoms$atom_id == from,
    ,
    drop = FALSE
  ]

  b <- atoms[
    atoms$atom_id == to,
    ,
    drop = FALSE
  ]

  if (nrow(a) == 0 || nrow(b) == 0) {
    return(1L)
  }

  # ----------------------------------------------
  # Local substituent information
  # ----------------------------------------------

  subs <- .mb_substituents(
    from = from,
    to = to,
    bonds = bonds,
    atoms = atoms
  )

  s_from <- subs$from
  s_to   <- subs$to

  n_from <- length(s_from)
  n_to   <- length(s_to)

  # ----------------------------------------------
  # Ring information
  # ----------------------------------------------

  ring <- .mb_ring_info(
    from = from,
    to = to,
    bonds = bonds
  )

  # ------------------------------------------------
  # GR-1.10.2
  #
  # Two or more substituents on one end and none
  # on the other.
  #
  # IUPAC recommends a centered representation.
  # ggchemplot intentionally uses an offset style
  # throughout, so choose the offset direction
  # toward the substituent-bearing end.
  # ------------------------------------------------

  if (n_from >= 2 && n_to == 0) {

    sides <- .mb_substituent_sides(
      from = from,
      to = to,
      substituents = s_from,
      atoms = atoms
    )

    if (length(sides) > 0) {
      score <- sum(sign(sides))

      if (score > 0) return(1L)
      if (score < 0) return(-1L)
    }
  }

  if (n_to >= 2 && n_from == 0) {

    sides <- .mb_substituent_sides(
      from = from,
      to = to,
      substituents = s_to,
      atoms = atoms
    )

    if (length(sides) > 0) {
      score <- sum(sign(sides))

      if (score > 0) return(1L)
      if (score < 0) return(-1L)
    }
  }

  # ----------------------------------------------
  # GR-1.10.4 / ring case
  #
  # If the multiple bond is endocyclic, the
  # secondary line should normally point toward
  # the ring interior.
  # ----------------------------------------------

  if (ring$in_ring) {

    centroid <- .mb_ring_centroid(
      ring_atoms = ring$ring_atoms,
      atoms = atoms
    )

    if (all(is.finite(centroid))) {

      return(
        .mb_side_toward_point(
          a$x, a$y,
          b$x, b$y,
          centroid[1],
          centroid[2]
        )
      )
    }
  }

  # ----------------------------------------------
  # GR-1.10.1
  #
  # Asymmetric substitution:
  #
  # offset toward the side having more
  # substituents.
  # ----------------------------------------------

  if (n_from != n_to) {

    # More substituents at FROM
    if (n_from > n_to) {

      sides <- .mb_substituent_sides(
        from = from,
        to = to,
        substituents = s_from,
        atoms = atoms
      )

      if (length(sides) > 0) {

        # Sum signed local contributions.
        # The sign tells us which side the substituent
        # occupies.
        score <- sum(sign(sides))

        if (score > 0) return(1L)
        if (score < 0) return(-1L)
      }
    }

    # More substituents at TO
    if (n_to > n_from) {

      sides <- .mb_substituent_sides(
        from = from,
        to = to,
        substituents = s_to,
        atoms = atoms
      )

      if (length(sides) > 0) {

        score <- sum(sign(sides))

        if (score > 0) return(1L)
        if (score < 0) return(-1L)
      }
    }
  }

  # ----------------------------------------------
  # GR-1.10.3
  #
  # One substituent on each end.
  #
  # If trans, choose the side occupied by one
  # substituent. For an acyclic chain either
  # direction is acceptable; choosing the side
  # of the FROM substituent makes the result
  # deterministic.
  # ----------------------------------------------

  if (n_from == 1 && n_to == 1) {

    if (.mb_is_trans(
      from = from,
      to = to,
      s_from = s_from,
      s_to = s_to,
      atoms = atoms
    )) {

      sides <- .mb_substituent_sides(
        from = from,
        to = to,
        substituents = s_from,
        atoms = atoms
      )

      if (length(sides) == 1 &&
          abs(sides) > .mb_eps) {

        return(
          ifelse(sides > 0, 1L, -1L)
        )
      }
    }
  }

  # ----------------------------------------------
  # GR-1.10.4
  #
  # Four-substituent acyclic bond.
  #
  # IUPAC allows either direction for acyclic
  # systems. We nevertheless make the result
  # deterministic by examining all substituents.
  # ----------------------------------------------

  if (n_from >= 2 && n_to >= 2) {

    all_subs <- c(s_from, s_to)

    sides <- .mb_substituent_sides(
      from = from,
      to = to,
      substituents = all_subs,
      atoms = atoms
    )

    if (length(sides) > 0) {

      score <- sum(sign(sides))

      if (score > 0) return(1L)
      if (score < 0) return(-1L)
    }
  }

  # ----------------------------------------------
  # General local geometric fallback
  #
  # Prefer the side with the strongest local
  # substituent signal.
  # ----------------------------------------------

  all_subs <- c(s_from, s_to)

  if (length(all_subs) > 0) {

    sides <- .mb_substituent_sides(
      from = from,
      to = to,
      substituents = all_subs,
      atoms = atoms
    )

    if (length(sides) > 0) {

      score <- sum(sign(sides))

      if (score > 0) return(1L)
      if (score < 0) return(-1L)
    }
  }

  # ----------------------------------------------
  # Completely symmetric / ambiguous case.
  #
  # Either side is chemically equivalent.
  # Choose one deterministically.
  # ----------------------------------------------

  1L
}
