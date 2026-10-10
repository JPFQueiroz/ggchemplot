#' @title Sketch a structure by adding named groups
#'
#' @description
#' \code{sketch_molecule()} starts from a ring or a single atom.
#' \code{add_group()} attaches a named organic or biochemical group to an
#' atom id. Placement is outward from the existing neighbours unless
#' \code{orientation} is given in degrees (0 = right, 90 = up). Implicit
#' hydrogens are recalculated from valence.
#'
#' @param template Ring or atom template. See \code{\link{sketch_templates}}.
#' @param bond_length Bond length in data units. Default \code{1} for a new
#'   sketch, or the median heavy-atom bond when adding a group.
#' @param charge Integer. Charge placed on the first atom of the template.
#'   Default \code{0}.
#' @param result A sketch result from \code{\link{sketch_molecule}}.
#' @param atom_id Atom that receives the group.
#' @param group Group name. See \code{\link{sketch_groups}}.
#' @param orientation Direction in degrees, or \code{NULL} for outward.
#' @return A ggchemplot result list.
#' @export
sketch_molecule <- function(template = "benzene", bond_length = 1, charge = 0L) {
  template <- match.arg(template, sketch_templates())
  frag <- .sk_template(template, bond_length)
  if (charge != 0L && nrow(frag$atoms)) frag$atoms$charge[1] <- as.integer(charge)
  .sk_finish(frag$atoms, frag$bonds, bond_length)
}

#' @rdname sketch_molecule
#' @export
add_group <- function(result, atom_id, group, orientation = NULL, bond_length = NULL) {
  group <- match.arg(group, sketch_groups())
  len <- bond_length %||% .sk_length(result, if (group %in% c("hydrogen", "H")) "H" else "C")
  parent <- result$atoms[match(atom_id, result$atoms$atom_id), ]
  if (!nrow(parent)) stop("Atom ", atom_id, " not found.")
  frag <- .sk_group(group, len)
  opened <- integer()
  if (.sk_free(result, atom_id) < frag$order) {
    dbl <- result$bond_coords$order >= 2 &
      (result$bond_coords$from == atom_id | result$bond_coords$to == atom_id)
    if (!any(dbl)) stop("Atom ", atom_id, " has no free valence for ", group, ".")
    idx <- which(dbl)[1]
    opened <- ifelse(result$bond_coords$from[idx] == atom_id,
                     result$bond_coords$to[idx],
                     result$bond_coords$from[idx])
    result$bond_coords$order[idx] <- result$bond_coords$order[idx] - 1L
  }
  h_ids <- unique(c(
    result$bond_coords$to[result$bond_coords$from == atom_id],
    result$bond_coords$from[result$bond_coords$to == atom_id]
  ))
  h_ids <- h_ids[result$atoms$symbol[match(h_ids, result$atoms$atom_id)] == "H"]
  if (length(h_ids)) {
    drop_id <- h_ids[1]
    result$atoms <- result$atoms[result$atoms$atom_id != drop_id, ]
    result$bond_coords <- result$bond_coords[
      result$bond_coords$from != drop_id & result$bond_coords$to != drop_id, ]
  }
  dir <- .sk_direction(result, atom_id, orientation)
  xy <- cbind(
    dir[1] * frag$atoms$x - dir[2] * frag$atoms$y,
    dir[2] * frag$atoms$x + dir[1] * frag$atoms$y
  )
  frag$atoms$x <- parent$x + xy[, 1]
  frag$atoms$y <- parent$y + xy[, 2]
  if (nrow(frag$atoms) > 1L) {
    centre <- c(mean(result$atoms$x), mean(result$atoms$y))
    axis <- c(frag$atoms$x[1] - parent$x, frag$atoms$y[1] - parent$y)
    axis <- axis / sqrt(sum(axis^2))
    rel <- cbind(frag$atoms$x - parent$x, frag$atoms$y - parent$y)
    along <- as.numeric(rel %*% axis)
    perp <- rel[, 1] * axis[2] - rel[, 2] * axis[1]
    flipped <- cbind(
      parent$x + along * axis[1] + perp * axis[2],
      parent$y + along * axis[2] - perp * axis[1]
    )
    now <- sum((c(frag$atoms$x[2], frag$atoms$y[2]) - centre)^2)
    other <- sum((flipped[2, ] - centre)^2)
    if (other > now) {
      frag$atoms$x <- flipped[, 1]
      frag$atoms$y <- flipped[, 2]
    }
  }
  shift <- max(result$atoms$atom_id)
  frag$atoms$atom_id <- frag$atoms$atom_id + shift
  if (nrow(frag$bonds)) {
    frag$bonds$from <- frag$bonds$from + shift
    frag$bonds$to <- frag$bonds$to + shift
  }
  link <- .sk_bond(atom_id, frag$attach + shift, frag$order)
  atoms <- rbind(result$atoms, frag$atoms)
  bonds <- rbind(result$bond_coords[, .sk_bond_cols()], frag$bonds, link)
  if (length(opened) && !isTRUE(result$params$collapse_hydrogens)) {
    partner <- atoms[atoms$atom_id == opened, ]
    id <- max(atoms$atom_id) + 1L
    hdir <- .sk_h_direction(list(atoms = atoms, bond_coords = bonds), opened, 1L, 2L)
    atoms <- rbind(atoms, .sk_atom(id, partner$x + hdir[1] * len * 0.7,
                                   partner$y + hdir[2] * len * 0.7, "H"))
    bonds <- rbind(bonds, .sk_bond(opened, id, 1L))
  }
  .sk_finish(atoms, bonds, len, result$params)
}

#' @rdname sketch_molecule
#' @export
sketch_groups <- function() {
  c("methyl", "ethyl", "propyl", "isopropyl", "tertbutyl", "vinyl", "phenyl",
    "hydroxy", "methoxy", "amino", "methylamino", "thiol",
    "formyl", "acetyl", "carboxyl", "carboxylate", "amide",
    "nitro", "cyano", "fluoro", "chloro", "bromo", "iodo",
    "phosphate", "sulfate")
}

#' @rdname sketch_molecule
#' @export
sketch_templates <- function() {
  c("atom", "benzene", "pyridine", "pyrimidine", "pyrrole", "imidazole",
    "furan", "thiophene", "cyclohexane", "cyclopentane", "cyclobutane",
    "naphthalene", "indole", "purine")
}

.sk_template <- function(template, len) {
  ring <- function(sides, hetero = NULL, orders = NULL) {
    radius <- len / (2 * sin(pi / sides))
    ang <- pi / 2 + seq(0, sides - 1) * 2 * pi / sides
    sym <- rep("C", sides)
    if (!is.null(hetero)) sym[hetero$at] <- hetero$el
    atoms <- do.call(rbind, lapply(seq_len(sides), function(i) {
      .sk_atom(i, radius * cos(ang[i]), radius * sin(ang[i]), sym[i])
    }))
    bonds <- do.call(rbind, lapply(seq_len(sides), function(i) {
      j <- if (i == sides) 1L else i + 1L
      .sk_bond(i, j, if (is.null(orders)) 1L else orders[i])
    }))
    list(atoms = atoms, bonds = bonds)
  }
  if (template == "atom") return(list(atoms = .sk_atom(1L, 0, 0, "C"), bonds = .sk_bonds0()))
  if (template == "benzene") return(ring(6, orders = rep(c(2L, 1L), 3)))
  if (template == "pyridine") return(ring(6, list(at = 1, el = "N"), rep(c(2L, 1L), 3)))
  if (template == "pyrimidine") return(ring(6, list(at = c(1, 3), el = c("N", "N")), rep(c(2L, 1L), 3)))
  if (template == "pyrrole") return(ring(5, list(at = 1, el = "N"), c(1L, 2L, 1L, 2L, 1L)))
  if (template == "imidazole") return(ring(5, list(at = c(1, 3), el = c("N", "N")), c(1L, 2L, 1L, 2L, 1L)))
  if (template == "furan") return(ring(5, list(at = 1, el = "O"), c(1L, 2L, 1L, 2L, 1L)))
  if (template == "thiophene") return(ring(5, list(at = 1, el = "S"), c(1L, 2L, 1L, 2L, 1L)))
  if (template == "cyclohexane") return(ring(6))
  if (template == "cyclopentane") return(ring(5))
  if (template == "cyclobutane") return(ring(4))
  if (template == "naphthalene") return(.sk_naphthalene(len))
  if (template == "indole") return(.sk_indole(len))
  if (template == "purine") return(.sk_purine(len))
  .sk_fuse(
    ring(6, list(at = c(1, 3), el = c("N", "N")), rep(c(2L, 1L), 3)),
    ring(5, list(at = c(1, 3), el = c("N", "N")), c(1L, 2L, 1L, 2L, 1L)),
    len
  )
  a <- ring(6, list(at = c(1, 3), el = c("N", "N")), rep(c(2L, 1L), 3))
  b <- ring(5, list(at = c(1, 3), el = c("N", "N")), c(1L, 2L, 1L, 2L, 1L))
  b$atoms$x <- b$atoms$x + 2 * len
  b$atoms$atom_id <- b$atoms$atom_id + 6L
  b$bonds$from <- b$bonds$from + 6L
  b$bonds$to <- b$bonds$to + 6L
  list(atoms = rbind(a$atoms, b$atoms), bonds = rbind(a$bonds, b$bonds, .sk_bond(3L, 9L, 1L)))
}

.sk_naphthalene <- function(len) {
  a <- .sk_template("benzene", len)
  p6 <- c(a$atoms$x[6], a$atoms$y[6])
  p5 <- c(a$atoms$x[5], a$atoms$y[5])
  mid <- (p6 + p5) / 2
  out <- c(1, 0)
  centre <- mid + out * len * sqrt(3) / 2
  ang <- atan2(p6[2] - centre[2], p6[1] - centre[1])
  extra <- do.call(rbind, lapply(1:4, function(i) {
    th <- ang - i * pi / 3
    .sk_atom(6L + i, centre[1] + len * cos(th), centre[2] + len * sin(th), "C")
  }))
  bonds <- rbind(
    a$bonds,
    .sk_bond(6L, 7L, 1L), .sk_bond(7L, 8L, 2L), .sk_bond(8L, 9L, 1L),
    .sk_bond(9L, 10L, 2L), .sk_bond(10L, 5L, 1L)
  )
  list(atoms = rbind(a$atoms, extra), bonds = bonds)
}

.sk_indole <- function(len) {
  a <- .sk_template("benzene", len)
  p6 <- c(a$atoms$x[6], a$atoms$y[6])
  p5 <- c(a$atoms$x[5], a$atoms$y[5])
  mid <- (p6 + p5) / 2
  d <- len / (2 * tan(pi / 5))
  r <- len / (2 * sin(pi / 5))
  centre <- mid + c(d, 0)
  ang <- atan2(p6[2] - centre[2], p6[1] - centre[1])
  sym <- c("C", "N", "C")
  extra <- do.call(rbind, lapply(1:3, function(i) {
    th <- ang - i * 2 * pi / 5
    .sk_atom(6L + i, centre[1] + r * cos(th), centre[2] + r * sin(th), sym[i])
  }))
  bonds <- rbind(
    a$bonds,
    .sk_bond(6L, 7L, 1L), .sk_bond(7L, 8L, 2L),
    .sk_bond(8L, 9L, 1L), .sk_bond(9L, 5L, 1L)
  )
  list(atoms = rbind(a$atoms, extra), bonds = bonds)
}

.sk_purine <- function(len) {
  a <- .sk_template("pyrimidine", len)
  p6 <- c(a$atoms$x[6], a$atoms$y[6])
  p5 <- c(a$atoms$x[5], a$atoms$y[5])
  mid <- (p6 + p5) / 2
  d <- len / (2 * tan(pi / 5))
  r <- len / (2 * sin(pi / 5))
  centre <- mid + c(d, 0)
  ang <- atan2(p6[2] - centre[2], p6[1] - centre[1])
  sym <- c("N", "C", "N")
  extra <- do.call(rbind, lapply(1:3, function(i) {
    th <- ang - i * 2 * pi / 5
    .sk_atom(6L + i, centre[1] + r * cos(th), centre[2] + r * sin(th), sym[i])
  }))
  bonds <- rbind(
    a$bonds,
    .sk_bond(6L, 7L, 1L), .sk_bond(7L, 8L, 2L),
    .sk_bond(8L, 9L, 1L), .sk_bond(9L, 5L, 1L)
  )
  list(atoms = rbind(a$atoms, extra), bonds = bonds)
}

.sk_fuse <- function(a, b, len) {
  right <- order(a$atoms$x, decreasing = TRUE)[1:2]
  shared <- a$atoms[right, ]
  mid <- c(mean(shared$x), mean(shared$y))
  edge <- c(shared$x[1] - shared$x[2], shared$y[1] - shared$y[2])
  normal <- c(-edge[2], edge[1])
  normal <- normal / sqrt(sum(normal^2))
  if (sum(normal * mid) < 0) normal <- -normal
  gap <- len / (2 * tan(pi / nrow(b$atoms)))
  centre <- mid + normal * gap
  ang0 <- atan2(shared$y[1] - centre[2], shared$x[1] - centre[1])
  step <- 2 * pi / nrow(b$atoms)
  radius <- len / (2 * sin(pi / nrow(b$atoms)))
  keep <- 2:(nrow(b$atoms) - 1L)
  new_atoms <- do.call(rbind, lapply(keep, function(i) {
    th <- ang0 + (i - 1L) * step
    .sk_atom(i, centre[1] + radius * cos(th), centre[2] + radius * sin(th), b$atoms$symbol[i])
  }))
  shift <- max(a$atoms$atom_id)
  new_atoms$atom_id <- seq_len(nrow(new_atoms)) + shift
  ids <- c(shared$atom_id[1], new_atoms$atom_id, shared$atom_id[2])
  bonds <- do.call(rbind, lapply(seq_len(length(ids) - 1L), function(i) {
    .sk_bond(ids[i], ids[i + 1L], 1L)
  }))
  list(atoms = rbind(a$atoms, new_atoms), bonds = rbind(a$bonds, bonds))
}

.sk_group <- function(group, len) {
  carbon <- function(n, bend = 60) .sk_chain("C", n, len, bend)
  one <- function(el, order = 1L) {
    list(atoms = .sk_atom(1L, len, 0, el), bonds = .sk_bonds0(), attach = 1L, order = order)
  }
  if (group == "methyl") return(one("C"))
  if (group == "ethyl") return(c(carbon(2), list(attach = 1L, order = 1L)))
  if (group == "propyl") return(c(carbon(3), list(attach = 1L, order = 1L)))
  if (group == "isopropyl") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "C"),
      .sk_atom(3L, len + len * cos(-pi / 3), len * sin(-pi / 3), "C")
    )
    bonds <- rbind(.sk_bond(1L, 2L, 1L), .sk_bond(1L, 3L, 1L))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "tertbutyl") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, 2 * len, 0, "C"),
      .sk_atom(3L, len, len, "C"),
      .sk_atom(4L, len, -len, "C")
    )
    bonds <- rbind(
      .sk_bond(1L, 2L, 1L),
      .sk_bond(1L, 3L, 1L),
      .sk_bond(1L, 4L, 1L)
    )
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "vinyl") {
    g <- carbon(2)
    g$bonds$order[1] <- 2L
    return(c(g, list(attach = 1L, order = 1L)))
  }
  if (group == "phenyl") {
    radius <- len
    centre <- c(len + radius, 0)
    ang <- pi + seq(0, 5) * 2 * pi / 6
    atoms <- do.call(rbind, lapply(seq_len(6), function(i) {
      .sk_atom(i, centre[1] + radius * cos(ang[i]), centre[2] + radius * sin(ang[i]), "C")
    }))
    bonds <- do.call(rbind, lapply(seq_len(6), function(i) {
      j <- if (i == 6L) 1L else i + 1L
      .sk_bond(i, j, if (i %% 2L == 1L) 2L else 1L)
    }))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "hydroxy") return(one("O"))
  if (group == "methoxy") return(.sk_pair("O", "C", len))
  if (group == "amino") return(one("N"))
  if (group == "methylamino") return(.sk_pair("N", "C", len))
  if (group == "thiol") return(one("S"))
  if (group == "fluoro") return(one("F"))
  if (group == "chloro") return(one("Cl"))
  if (group == "bromo") return(one("Br"))
  if (group == "iodo") return(one("I"))
  if (group == "cyano") {
    g <- .sk_pair("C", "N", len)
    g$bonds$order[1] <- 3L
    g$order <- 1L
    return(g)
  }
  if (group == "formyl") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "O")
    )
    bonds <- .sk_bond(1L, 2L, 2L)
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "acetyl") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "O"),
      .sk_atom(3L, len + len * cos(-pi / 3), len * sin(-pi / 3), "C")
    )
    bonds <- rbind(.sk_bond(1L, 2L, 2L), .sk_bond(1L, 3L, 1L))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "carboxyl") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "O"),
      .sk_atom(3L, len + len * cos(-pi / 3), len * sin(-pi / 3), "O")
    )
    bonds <- rbind(.sk_bond(1L, 2L, 2L), .sk_bond(1L, 3L, 1L))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "carboxylate") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "O"),
      .sk_atom(3L, len + len * cos(-pi / 3), len * sin(-pi / 3), "O", -1L)
    )
    bonds <- rbind(.sk_bond(1L, 2L, 2L), .sk_bond(1L, 3L, 1L))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "amide") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "C"),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "O"),
      .sk_atom(3L, len + len * cos(-pi / 3), len * sin(-pi / 3), "N")
    )
    bonds <- rbind(.sk_bond(1L, 2L, 2L), .sk_bond(1L, 3L, 1L))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "nitro") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "N", 1L),
      .sk_atom(2L, len + len * cos(pi / 3), len * sin(pi / 3), "O"),
      .sk_atom(3L, len + len * cos(-pi / 3), len * sin(-pi / 3), "O", -1L)
    )
    bonds <- rbind(.sk_bond(1L, 2L, 2L), .sk_bond(1L, 3L, 1L))
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "phosphate") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "P"),
      .sk_atom(2L, len + len * cos(pi / 2), len * sin(pi / 2), "O"),
      .sk_atom(3L, len + len * cos(-pi / 2), len * sin(-pi / 2), "O", -1L),
      .sk_atom(4L, len + len, 0, "O", -1L)
    )
    bonds <- rbind(
      .sk_bond(1L, 2L, 2L),
      .sk_bond(1L, 3L, 1L),
      .sk_bond(1L, 4L, 1L)
    )
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
  if (group == "sulfate") {
    atoms <- rbind(
      .sk_atom(1L, len, 0, "S"),
      .sk_atom(2L, len + len * cos(pi / 2), len * sin(pi / 2), "O"),
      .sk_atom(3L, len + len * cos(-pi / 2), len * sin(-pi / 2), "O", -1L),
      .sk_atom(4L, len + len, 0, "O", -1L)
    )
    bonds <- rbind(
      .sk_bond(1L, 2L, 2L),
      .sk_bond(1L, 3L, 1L),
      .sk_bond(1L, 4L, 1L)
    )
    return(list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L))
  }
}

.sk_chain <- function(el, n, len, bend) {
  atoms <- vector("list", n)
  x <- len
  y <- 0
  heading <- 0
  for (i in seq_len(n)) {
    atoms[[i]] <- .sk_atom(i, x, y, el)
    heading <- heading + bend * pi / 180 * if (i %% 2L == 1L) 1 else -1
    x <- x + len * cos(heading)
    y <- y + len * sin(heading)
  }
  bonds <- if (n == 1L) .sk_bonds0() else {
    do.call(rbind, lapply(seq_len(n - 1L), function(i) .sk_bond(i, i + 1L, 1L)))
  }
  list(atoms = do.call(rbind, atoms), bonds = bonds)
}

.sk_branched <- function(len, tips) {
  atoms <- .sk_atom(1L, len, 0, "C")
  bonds <- .sk_bonds0()
  for (i in seq_along(tips)) {
    ang <- (i - 1) * 2 * pi / length(tips)
    atoms <- rbind(atoms, .sk_atom(i + 1L, len + len * cos(ang), len * sin(ang), tips[i]))
    bonds <- rbind(bonds, .sk_bond(1L, i + 1L, 1L))
  }
  list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L)
}

.sk_pair <- function(a, b, len) {
  list(atoms = rbind(.sk_atom(1L, len, 0, a), .sk_atom(2L, 2 * len, 0.4 * len, b)),
       bonds = .sk_bond(1L, 2L, 1L), attach = 1L, order = 1L)
}

.sk_carbonyl <- function(len, extra) {
  atoms <- rbind(.sk_atom(1L, len, 0, "C"), .sk_atom(2L, len, len, "O"))
  bonds <- .sk_bond(1L, 2L, 2L)
  if (!is.null(extra)) {
    atoms <- rbind(atoms, .sk_atom(3L, 2 * len, -0.4 * len, extra))
    bonds <- rbind(bonds, .sk_bond(1L, 3L, 1L))
  }
  list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L)
}

.sk_oxo_acid <- function(el, len, n_o) {
  atoms <- .sk_atom(1L, len, 0, el)
  bonds <- .sk_bonds0()
  for (i in seq_len(n_o)) {
    ang <- -pi / 2 + (i - 1) * pi / max(1, n_o - 1)
    atoms <- rbind(atoms, .sk_atom(i + 1L, len + len * cos(ang), len * sin(ang), "O"))
    bonds <- rbind(bonds, .sk_bond(1L, i + 1L, if (i == 1L) 2L else 1L))
  }
  list(atoms = atoms, bonds = bonds, attach = 1L, order = 1L)
}

.sk_free <- function(result, atom_id) {
  row <- result$atoms[match(atom_id, result$atoms$atom_id), ]
  bonds <- result$bond_coords
  if (!nrow(bonds)) return(.sk_valence(row$symbol))
  other <- ifelse(bonds$from == atom_id, bonds$to, ifelse(bonds$to == atom_id, bonds$from, NA))
  hit <- !is.na(other)
  other <- other[hit]
  orders <- bonds$order[hit]
  sym <- result$atoms$symbol[match(other, result$atoms$atom_id)]
  used <- sum(orders[sym != "H"], na.rm = TRUE)
  .sk_valence(row$symbol) - used
}

.sk_valence <- function(symbol) {
  switch(symbol, C = 4L, N = 3L, O = 2L, S = 2L, P = 5L, H = 1L, 1L)
}

.sk_used <- function(bonds, atom_id) {
  if (!nrow(bonds)) return(0L)
  hit <- bonds$from == atom_id | bonds$to == atom_id
  if (!any(hit)) 0L else sum(bonds$order[hit])
}

.sk_direction <- function(result, atom_id, orientation) {
  if (!is.null(orientation)) {
    th <- orientation * pi / 180
    return(c(cos(th), sin(th)))
  }
  parent <- result$atoms[match(atom_id, result$atoms$atom_id), ]
  heavy <- result$atoms[result$atoms$symbol != "H", ]
  bonds <- result$bond_coords
  nbs <- integer()
  if (nrow(bonds)) nbs <- c(bonds$to[bonds$from == atom_id], bonds$from[bonds$to == atom_id])
  nbs <- nbs[result$atoms$symbol[match(nbs, result$atoms$atom_id)] != "H"]
  nbs <- nbs[!is.na(nbs)]
  if (!length(nbs)) return(c(1, 0))
  angs <- vapply(nbs, function(nb) {
    other <- result$atoms[match(nb, result$atoms$atom_id), ]
    atan2(other$y - parent$y, other$x - parent$x)
  }, numeric(1))
  if (length(angs) == 1L) {
    turn <- if (atom_id %% 2L == 1L) pi / 3 else -pi / 3
    th <- angs[1] + pi + turn
    return(c(cos(th), sin(th)))
  }
  angs <- sort((angs + pi) %% (2 * pi) - pi)
  gaps <- c(diff(angs), angs[1] + 2 * pi - angs[length(angs)])
  k <- which.max(gaps)
  th <- if (k == length(angs)) angs[length(angs)] + gaps[k] / 2 else angs[k] + gaps[k] / 2
  c(cos(th), sin(th))
}

.sk_length <- function(result, symbol = "C") {
  atoms <- result$atoms
  bonds <- result$bond_coords
  if (!nrow(bonds)) return(result$params$target_bond_length %||% 1)
  len <- sqrt((bonds$x2 - bonds$x1)^2 + (bonds$y2 - bonds$y1)^2)
  s1 <- atoms$symbol[match(bonds$from, atoms$atom_id)]
  s2 <- atoms$symbol[match(bonds$to, atoms$atom_id)]
  hetero <- c("O", "N", "S", "P")
  h_het <- (s1 == "H" & s2 %in% hetero) | (s2 == "H" & s1 %in% hetero)
  use <- if (identical(symbol, "H")) h_het else !h_het & s1 != "H" & s2 != "H"
  vals <- len[use & is.finite(len) & len > 1e-6]
  if (!length(vals)) vals <- len[is.finite(len) & len > 1e-6]
  if (!length(vals)) return(1)
  med <- median(vals)
  vals <- vals[vals > 0.5 * med & vals < 1.5 * med]
  if (length(vals)) median(vals) else med
}

.sk_align_bonds <- function(new_bonds, template, atoms) {
  if (!nrow(new_bonds)) return(new_bonds)
  missing <- setdiff(names(template), names(new_bonds))
  for (col in missing) new_bonds[[col]] <- NA
  new_bonds <- new_bonds[, names(template), drop = FALSE]
  new_bonds$x1 <- atoms$x[match(new_bonds$from, atoms$atom_id)]
  new_bonds$y1 <- atoms$y[match(new_bonds$from, atoms$atom_id)]
  new_bonds$x2 <- atoms$x[match(new_bonds$to, atoms$atom_id)]
  new_bonds$y2 <- atoms$y[match(new_bonds$to, atoms$atom_id)]
  new_bonds$sym1 <- atoms$symbol[match(new_bonds$from, atoms$atom_id)]
  new_bonds$sym2 <- atoms$symbol[match(new_bonds$to, atoms$atom_id)]
  new_bonds
}

.sk_atom <- function(id, x, y, symbol, charge = 0L) {
  if (!length(symbol)) return(.sk_empty_atoms())
  symbol <- as.character(symbol)
  col <- vapply(symbol, function(s) {
    switch(s, C = "black", O = "#e31a1c", N = "#1f78b4",
           P = "violet", S = "#e6ab02", H = "gray40", "black")
  }, character(1), USE.NAMES = FALSE)
  data.frame(atom_id = as.integer(id), x = x, y = y, symbol = symbol,
             default_color = col, color = col, charge = as.integer(charge),
             radical = 0L, show_label = symbol != "C", stringsAsFactors = FALSE)
}

.sk_empty_atoms <- function() {
  data.frame(atom_id = integer(), x = numeric(), y = numeric(), symbol = character(),
             default_color = character(), color = character(), charge = integer(),
             radical = integer(), show_label = logical(), stringsAsFactors = FALSE)
}

.sk_bond <- function(from, to, order) {
  data.frame(from = as.integer(from), to = as.integer(to), order = as.integer(order),
             stereo = 0L, bond_type = NA_character_, shorten_start = NA_real_,
             shorten_end = NA_real_, width = NA_real_, wedge_thickness = NA_real_,
             n_hashes = NA_integer_, colour = NA_character_, stringsAsFactors = FALSE)
}

.sk_bonds0 <- function() {
  data.frame(
    from = integer(), to = integer(), order = integer(), stereo = integer(),
    bond_type = character(), shorten_start = numeric(), shorten_end = numeric(),
    width = numeric(), wedge_thickness = numeric(), n_hashes = integer(),
    colour = character(), stringsAsFactors = FALSE
  )
}
.sk_bond_cols <- function() names(.sk_bond(1L, 2L, 1L))

#' @title Collapse or restore sketch hydrogens
#'
#' @description
#' \code{collapse_sketch()} removes explicit hydrogens from the working
#' tables and stores them as collapsed labels. \code{uncollapse_sketch()}
#' restores them as atoms and bonds.
#'
#' @param result A sketch result from \code{\link{sketch_molecule}}.
#' @param atom_id Optional atom ids. If given, only hydrogens on those
#'   atoms are restored. Default \code{NULL} restores all.
#' @return The updated result.
#' @export
collapse_sketch <- function(result) {
  h_ids <- result$atoms$atom_id[result$atoms$symbol == "H"]
  result$atoms <- result$atoms[!result$atoms$atom_id %in% h_ids, ]
  result$bond_coords <- result$bond_coords[
    !result$bond_coords$from %in% h_ids & !result$bond_coords$to %in% h_ids, ]
  explicit <- .sk_explicit_h(result)
  h_atoms <- explicit$atoms
  h_atoms$parent_id <- NULL
  result$original_atoms <- rbind(result$atoms, h_atoms)
  result$original_bond_coords <- rbind(
    result$bond_coords,
    .sk_align_bonds(explicit$bonds, result$bond_coords, result$atoms)
  )
  result$h_labels <- explicit$labels
  result$params$collapse_hydrogens <- TRUE
  result$plot <- ggchemplot2(result)
  result
}


#' @rdname collapse_sketch
#' @export
uncollapse_sketch <- function(result, atom_id = NULL) {
  explicit <- .sk_explicit_h(result)
  if (!is.null(atom_id)) {
    keep <- explicit$atoms$parent_id %in% atom_id
    explicit$atoms <- explicit$atoms[keep, ]
    explicit$bonds <- explicit$bonds[explicit$from %in% atom_id | explicit$to %in% atom_id, ]
    if (!is.null(result$h_labels)) {
      result$h_labels <- result$h_labels[!result$h_labels$parent_id %in% atom_id, ]
    }
  } else {
    result$h_labels <- NULL
  }
  if (nrow(explicit$atoms)) {
    explicit$atoms$parent_id <- NULL
    result$atoms <- rbind(result$atoms, explicit$atoms)
  }
  if (nrow(explicit$bonds)) {
    result$bond_coords <- rbind(
      result$bond_coords,
      .sk_align_bonds(explicit$bonds, result$bond_coords, result$atoms)
    )
  }
  result$params$collapse_hydrogens <- FALSE
  result$plot <- ggchemplot2(result)
  result
}

.sk_align_bonds <- function(new_bonds, template, atoms) {
  if (!nrow(new_bonds)) return(new_bonds)
  missing <- setdiff(names(template), names(new_bonds))
  for (col in missing) new_bonds[[col]] <- NA
  new_bonds <- new_bonds[, names(template), drop = FALSE]
  new_bonds$x1 <- atoms$x[match(new_bonds$from, atoms$atom_id)]
  new_bonds$y1 <- atoms$y[match(new_bonds$from, atoms$atom_id)]
  new_bonds$x2 <- atoms$x[match(new_bonds$to, atoms$atom_id)]
  new_bonds$y2 <- atoms$y[match(new_bonds$to, atoms$atom_id)]
  new_bonds$sym1 <- atoms$symbol[match(new_bonds$from, atoms$atom_id)]
  new_bonds$sym2 <- atoms$symbol[match(new_bonds$to, atoms$atom_id)]
  new_bonds
}

.sk_explicit_h <- function(result) {
  atoms <- result$atoms
  bonds <- result$bond_coords
  len <- result$params$target_bond_length %||% 1
  rows <- list()
  links <- list()
  labels <- list()
  next_id <- if (nrow(atoms)) max(atoms$atom_id) else 0L
  for (i in seq_len(nrow(atoms))) {
    parent <- atoms[i, ]
    nH <- max(0L, .sk_valence(parent$symbol) - .sk_used(bonds, parent$atom_id) + min(parent$charge, 0L))
    if (!nH) next
    labels[[length(labels) + 1L]] <- data.frame(
      parent_id = parent$atom_id, nH = nH,
      h_text = if (nH == 1L) "H" else paste0("H[", nH, "]"),
      x = parent$x, y = parent$y, color = "gray40",
      side_force = NA_character_, stringsAsFactors = FALSE
    )
    for (k in seq_len(nH)) {
      next_id <- next_id + 1L
      dir <- .sk_h_direction(result, parent$atom_id, k, nH)
      rows[[length(rows) + 1L]] <- .sk_atom(
        next_id, parent$x + dir[1] * len * 0.7, parent$y + dir[2] * len * 0.7, "H"
      )
      rows[[length(rows)]]$parent_id <- parent$atom_id
      links[[length(links) + 1L]] <- .sk_bond(parent$atom_id, next_id, 1L)
    }
  }
  list(
    atoms = if (length(rows)) do.call(rbind, rows) else .sk_atom(integer(), numeric(), numeric(), character())[0, ],
    bonds = if (length(links)) do.call(rbind, links) else .sk_bonds0(),
    labels = if (length(labels)) do.call(rbind, labels) else NULL
  )
}

.sk_h_direction <- function(result, atom_id, k, nH) {
  bonds <- result$bond_coords
  nbs <- integer()
  if (nrow(bonds)) nbs <- c(bonds$to[bonds$from == atom_id], bonds$from[bonds$to == atom_id])
  nbs <- nbs[result$atoms$symbol[match(nbs, result$atoms$atom_id)] != "H"]
  nbs <- nbs[!is.na(nbs)]
  if (!length(nbs)) {
    th <- pi / 2 - (k - 1) * 2 * pi / nH
    return(c(cos(th), sin(th)))
  }
  if (length(nbs) == 1L) {
    other <- result$atoms[match(nbs, result$atoms$atom_id), ]
    parent <- result$atoms[match(atom_id, result$atoms$atom_id), ]
    incoming <- atan2(other$y - parent$y, other$x - parent$x)
    offsets <- if (nH == 3L) {
      c(pi, pi - pi / 3, pi + pi / 3)
    } else if (nH == 2L) {
      c(pi - pi / 3, pi + pi / 3)
    } else {
      pi
    }
    th <- incoming + offsets[k]
    return(c(cos(th), sin(th)))
  }
  base <- .sk_direction(result, atom_id, NULL)
  th <- atan2(base[2], base[1]) + (k - (nH + 1) / 2) * pi / 3
  c(cos(th), sin(th))
}

.sk_finish <- function(atoms, bonds, bond_length, params = NULL) {
  if (!nrow(bonds)) bonds <- .sk_bonds0()
  bonds <- bonds[, .sk_bond_cols()]
  bonds$x1 <- atoms$x[match(bonds$from, atoms$atom_id)]
  bonds$y1 <- atoms$y[match(bonds$from, atoms$atom_id)]
  bonds$sym1 <- atoms$symbol[match(bonds$from, atoms$atom_id)]
  bonds$col1 <- atoms$color[match(bonds$from, atoms$atom_id)]
  bonds$x2 <- atoms$x[match(bonds$to, atoms$atom_id)]
  bonds$y2 <- atoms$y[match(bonds$to, atoms$atom_id)]
  bonds$sym2 <- atoms$symbol[match(bonds$to, atoms$atom_id)]
  nH <- vapply(atoms$atom_id, function(id) {
  row <- atoms[atoms$atom_id == id, ]
  max(0L, .sk_valence(row$symbol) - .sk_used(bonds, id) + min(row$charge, 0L))
}, integer(1))
  h <- atoms[nH > 0, ]
  h_labels <- if (!nrow(h)) NULL else {
    data.frame(parent_id = h$atom_id, nH = nH[nH > 0],
               h_text = ifelse(nH[nH > 0] == 1L, "H", paste0("H", nH[nH > 0])),
               x = h$x, y = h$y, color = "gray40", side_force = NA_character_,
               stringsAsFactors = FALSE)
  }
  if (is.null(params)) {
    params <- list(title = NULL, collapse_hydrogens = FALSE, rotation = 0,
                   flip_horizontal = FALSE, flip_vertical = FALSE, label_padding = 0.6,
                   show_atom_circles = TRUE, hide_carbon_circles = TRUE, circle_stroke = 0,
                   show_atom_labels = TRUE, hide_carbon_labels = FALSE, bond_width = NULL,
                   atom_size = NULL, label_size = 12, double_bond_offset = 0.15,
                   custom_atom_colors = NULL, atom_circle_color = "transparent",
                   paint_it_black = FALSE, H_offset = NULL, label_fontface = "plain",
                   target_bond_length = bond_length, normalize = FALSE, pubchem_cid = NULL)
  }
  result <- list(atoms = atoms, bond_coords = bonds, h_labels = h_labels,
                 original_atoms = atoms, original_bond_coords = bonds, params = params)
  if (isTRUE(params$collapse_hydrogens)) {
    result <- collapse_sketch(result)
  } else {
    result$h_labels <- NULL
    result$params$collapse_hydrogens <- FALSE
    explicit <- .sk_explicit_h(result)
    if (nrow(explicit$atoms)) {
      explicit$atoms$parent_id <- NULL
      result$atoms <- rbind(result$atoms, explicit$atoms)
      result$bond_coords <- rbind(
        result$bond_coords,
        .sk_align_bonds(explicit$bonds, result$bond_coords, result$atoms)
      )
    }
    result$plot <- ggchemplot2(result)
  }
  result
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
