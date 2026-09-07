#' @title Draw a 2D chemical structure with ggplot2
#'
#' @description
#' Draws a 2D structure from the atom and bond tables returned by
#' \code{\link{ggchemplot1}}. Bonds that end on a labelled atom are
#' shortened so they do not run through the letter. Multiple bonds
#' follow IUPAC GR-1.10 (offset second stroke; ends shortened by
#' local substitution). V2000 stereo is drawn with
#' \code{geom_wedge}. Charges and radicals from the SDF are
#' superscripts on the element symbol. Collapsed hydrogens are
#' placed beside the parent atom when
#' \code{collapse_hydrogens = TRUE} was set in \code{ggchemplot1()}.
#'
#' Arguments that are \code{NULL} are taken from \code{result$params},
#' then from the defaults below. Geometry fixed in
#' \code{\link{ggchemplot1}} (\code{collapse_hydrogens},
#' \code{rotation}, \code{flip_horizontal}, \code{flip_vertical},
#' \code{use_sdf_stereo}) cannot be changed here.
#'
#' @param result List. Object returned by \code{\link{ggchemplot1}}.
#' @param title Character. Optional plot title. Default: value stored
#'   in \code{result$params}.
#' @param show_atom_circles Logical. Draw discs behind atom labels.
#'   Default \code{TRUE}.
#' @param hide_carbon_circles Logical. Omit discs on carbon atoms.
#'   Default \code{TRUE}.
#' @param atom_size Numeric. Disc size in mm. If \code{NULL},
#'   \code{label_size * 0.353 * 2}.
#' @param atom_circle_color Character. Fill of atom discs. Default
#'   \code{"transparent"}.
#' @param circle_stroke Numeric. Disc border width. Default \code{0}
#'   from \code{ggchemplot1()}; if still unset, \code{0}.
#' @param show_atom_labels Logical. Draw element symbols. Default
#'   \code{TRUE}.
#' @param hide_carbon_labels Logical. Omit carbon symbols unless the
#'   atom has a charge or a radical. Default \code{TRUE}.
#' @param label_size Numeric. Atom label size in points. Default
#'   \code{12}.
#' @param label_fontface Character. One of \code{"plain"},
#'   \code{"bold"}, \code{"italic"}, \code{"bold.italic"}. Plotmath
#'   labels (charges) need the face wrapped in the expression; plain
#'   hydrogen labels use this argument directly. Default \code{"plain"}.
#' @param label_family Character. Font family for labels. Default
#'   \code{"sans"}. On Windows the built-in families are
#'   \code{"sans"}, \code{"serif"}, and \code{"mono"}.
#' @param label_padding Numeric. Extra data-unit margin so labels are
#'   not clipped. Default \code{0.6}.
#' @param H_offset Numeric vector of length 2: data-unit gap from the
#'   heteroatom to a collapsed \code{H} (first) or \code{H2}/\code{H3}
#'   (second). Used only when \code{collapse_hydrogens = TRUE} was set
#'   in \code{\link{ggchemplot1}}. If \code{NULL}, a gap is computed
#'   from \code{label_size} and \code{target_bond_length}.
#' @param bond_width Numeric. Bond width in mm. If \code{NULL},
#'   \code{label_size * 0.353 * 0.1}.
#' @param bond_color Character. Colour of bonds. Default \code{"black"}.
#' @param double_bond_offset Numeric. Perpendicular offset of the
#'   second (and third) stroke of a multiple bond, in data units.
#'   Default \code{0.15}.
#' @param paint_it_black Logical. If \code{TRUE}, atom colours are
#'   ignored and everything is drawn black. Default \code{FALSE}.
#' @param show_ids Logical. If \code{TRUE}, draw atom ids (red) and
#'   bond ids (blue) for editing. Default \code{FALSE}.
#' @param normalize Logical. If \code{TRUE}, scale coordinates so the
#'   median bond length equals \code{target_bond_length}. Default
#'   \code{TRUE}.
#' @param target_bond_length Numeric. Target median bond length in data
#'   units when \code{normalize = TRUE}. Default \code{1}.
#'
#' @return A \code{ggplot} object.
#'
#' @details
#' Arguments that are \code{NULL} are taken from \code{result$params},
#' then from the defaults above. Geometry fixed in
#' \code{\link{ggchemplot1}} (\code{collapse_hydrogens},
#' \code{rotation}, \code{flip_horizontal}, \code{flip_vertical},
#' \code{use_sdf_stereo}) cannot be changed here.
#'
#' Bonds that end on a labelled atom are shortened by
#' \code{atom_clearance} so they do not run through the letter.
#' Multiple-bond offset strokes use the same clearance or the IUPAC
#' GR-1.10 fractional shorten, whichever is larger. Charges from the
#' SDF are drawn as superscripts on the element symbol. A doublet
#' radical is drawn as a separate mark. Collapsed hydrogens are
#' placed left or right of the parent when possible.
#'
#' @seealso \code{\link{ggchemplot1}}, \code{\link{save_chemplot}},
#'   \code{\link{add_wedge_bond}}
#'
#' @import ggplot2
#' @import dplyr
#'
#' @importFrom stats setNames
#'
#' @export
ggchemplot2 <- function(result,
                        title = NULL,
                        show_atom_circles = NULL,
                        hide_carbon_circles = NULL,
                        atom_size = NULL,
                        atom_circle_color = NULL,
                        circle_stroke = NULL,
                        show_atom_labels = NULL,
                        hide_carbon_labels = NULL,
                        label_size = NULL,
                        label_fontface = NULL,
                        label_family = NULL,
                        label_padding = NULL,
                        H_offset = NULL,
                        bond_width = NULL,
                        bond_color = "black",
                        double_bond_offset = NULL,
                        paint_it_black = NULL,
                        show_ids = FALSE,
                        normalize = NULL,
                        target_bond_length = NULL) {

  pms <- result$params

  title               <- if (!is.null(title)) title else pms$title
  label_padding       <- if (!is.null(label_padding)) label_padding else pms$label_padding %||% 0.6
  show_atom_circles   <- if (!is.null(show_atom_circles)) show_atom_circles else pms$show_atom_circles %||% TRUE
  hide_carbon_circles <- if (!is.null(hide_carbon_circles)) hide_carbon_circles else pms$hide_carbon_circles %||% TRUE
  circle_stroke       <- if (!is.null(circle_stroke)) circle_stroke else pms$circle_stroke %||% 0
  show_atom_labels    <- if (!is.null(show_atom_labels)) show_atom_labels else pms$show_atom_labels %||% TRUE
  hide_carbon_labels  <- if (!is.null(hide_carbon_labels)) hide_carbon_labels else pms$hide_carbon_labels %||% TRUE
  label_size          <- if (!is.null(label_size)) label_size else pms$label_size %||% 12
  bond_width          <- if (!is.null(bond_width)) bond_width else pms$bond_width %||% (label_size * 0.353 * 0.1)
  atom_size           <- if (!is.null(atom_size)) atom_size else pms$atom_size %||% (label_size * 0.353 * 2)
  double_bond_offset  <- if (!is.null(double_bond_offset)) double_bond_offset else pms$double_bond_offset %||% 0.15
  paint_it_black      <- if (!is.null(paint_it_black)) paint_it_black else pms$paint_it_black %||% FALSE
  bond_color          <- if (!is.null(bond_color)) bond_color else pms$bond_color %||% "black"
  atom_circle_color   <- if (!is.null(atom_circle_color)) atom_circle_color else pms$atom_circle_color %||% "transparent"
  H_offset            <- if (!is.null(H_offset)) H_offset else pms$H_offset %||% c(-0.01250 + 0.02969 * label_size, 0.06375 + 0.03000 * label_size)
  label_fontface      <- if (!is.null(label_fontface)) label_fontface else pms$label_fontface %||% "plain"
  label_family        <- if (!is.null(label_family)) label_family else pms$label_family %||% "sans"
  normalize           <- if (!is.null(normalize)) normalize else pms$normalize %||% TRUE
  target_bond_length  <- if (!is.null(target_bond_length)) target_bond_length else pms$target_bond_length %||% 1

  if (isTRUE(normalize)) {
    result <- normalize_structure(result, target_bond_length = target_bond_length)
  }

  atoms <- result$atoms
  bond_coords <- result$bond_coords

  if (!"charge"     %in% names(atoms)) atoms$charge     <- 0L
  if (!"radical"    %in% names(atoms)) atoms$radical    <- 0L
  if (!"show_label" %in% names(atoms)) atoms$show_label <- FALSE
  atoms$show_label[is.na(atoms$show_label)] <- FALSE

  # ====================== ATOM-DISK CLEARANCE ======================
  atom_size_ref  <- 12 * 0.353 * 1.4
  clearance_ref  <- 0.22
  atom_clearance <- clearance_ref * (atom_size / atom_size_ref) * target_bond_length

  atoms_need_clearance <- atoms %>%
    mutate(
      is_carbon = .data$symbol == "C",
      has_charge_or_rad = .data$charge != 0 | .data$radical != 0,
      force_label = isTRUE(.data$show_label) | .data$show_label %in% TRUE,
      needs_clearance =
        (isTRUE(show_atom_circles) &
           (!(isTRUE(hide_carbon_circles) & .data$is_carbon) | .data$force_label)) |
        (isTRUE(show_atom_labels) &
           (!(isTRUE(hide_carbon_labels) & .data$is_carbon) |
              .data$has_charge_or_rad |
              .data$force_label))
    ) %>%
    select(.data$atom_id, .data$needs_clearance)

  if (nrow(bond_coords) > 0) {
    bond_coords <- bond_coords %>%
      left_join(atoms_need_clearance, by = c("from" = "atom_id")) %>%
      rename(clear_from = .data$needs_clearance) %>%
      left_join(atoms_need_clearance, by = c("to" = "atom_id")) %>%
      rename(clear_to = .data$needs_clearance) %>%
      mutate(
        clear_from = !is.na(.data$clear_from) & .data$clear_from,
        clear_to   = !is.na(.data$clear_to)   & .data$clear_to,
        dx  = .data$x2 - .data$x1,
        dy  = .data$y2 - .data$y1,
        len = sqrt(.data$dx^2 + .data$dy^2),
        ux  = ifelse(.data$len > 1e-8, .data$dx / .data$len, 0),
        uy  = ifelse(.data$len > 1e-8, .data$dy / .data$len, 0),
        s1  = ifelse(.data$clear_from, atom_clearance, 0),
        s2  = ifelse(.data$clear_to,   atom_clearance, 0),
        s1  = pmin(.data$s1, 0.40 * .data$len),
        s2  = pmin(.data$s2, 0.40 * .data$len),
        x1s = .data$x1 + .data$ux * .data$s1,
        y1s = .data$y1 + .data$uy * .data$s1,
        x2s = .data$x2 - .data$ux * .data$s2,
        y2s = .data$y2 - .data$uy * .data$s2
      )
  }

  p <- ggplot()

  # ===================== BONDS (Normal + Wedge + Hashed) ======================
  if (nrow(bond_coords) > 0) {

    normal_bonds <- bond_coords %>% filter(is.na(.data$bond_type))
    if (nrow(normal_bonds) > 0) {
      p <- p + geom_segment(
        data = normal_bonds,
        aes(x = .data$x1s, y = .data$y1s, xend = .data$x2s, yend = .data$y2s),
        linewidth = bond_width, color = bond_color
      )
    }

    # IUPAC GR-1.3: gap ≈ 3 × plain-bond visual thickness
    # bond_width is mm; convert with the same mm-per-data-unit used for circles
    mm_per_unit  <- (3.0 / 0.20) * (1.0 / target_bond_length)
    bond_width_data <- bond_width / mm_per_unit
    hash_gap <- pmax(0.10, 3 * bond_width_data)

    special_bonds <- bond_coords %>%
      filter(!is.na(.data$bond_type), .data$bond_type %in% c("solid", "hashed")) %>%
      mutate(
        draw_len = pmax(0.05, .data$len - .data$s1 - .data$s2),

        # half-width at the FAR end, in data units
        # ~0.13 when bond length ≈ 1  → much slimmer than 0.26
        wedge_w = pmin(0.16, pmax(0.09, 0.13 * .data$draw_len / target_bond_length)),

        # at least 3 hashes; typically 4–8 on a normal bond
        n_hashes_adapt = pmax(3L, pmin(9L, as.integer(round(.data$draw_len / hash_gap)))),

        # thickness of each hash bar (along the bond), data units
        hash_thick = pmax(0.025, 0.55 * bond_width_data)
      )

    if (nrow(special_bonds) > 0) {
      p <- p + geom_wedge(
        data = special_bonds,
        aes(x1 = .data$x1s, y1 = .data$y1s,
            x2 = .data$x2s, y2 = .data$y2s,
            wedgetype = .data$bond_type,
            width = .data$wedge_w,
            n_hashes = .data$n_hashes_adapt,
            wedge_thickness = .data$hash_thick),
        shorten_start = 0,
        shorten_end   = 0,
        colour = special_bonds$colour %||% bond_color,
        fill   = special_bonds$colour %||% bond_color
      )
    }
  }

  # PROTEIN LINKS (Wavy lines)
  if (!is.null(result$protein_links) && length(result$protein_links) > 0) {
    for (link in result$protein_links) {
      p <- p + geom_path(data = link$data,
                         aes(x = .data$x, y = .data$y),
                         linewidth = link$linewidth %||% 1.3,
                         color = bond_color)
    }
  }

  if (!is.null(result$protein_labels) && length(result$protein_labels) > 0) {
    labels_df <- do.call(rbind, result$protein_labels)
    p <- p + geom_text(data = labels_df,
                       aes(x = .data$x, y = .data$y, label = .data$label),
                       size = 4.2,
                       fontface = label_fontface,
                       color = "black")
  }

  # ====================== MULTIPLE BONDS ======================
  multi_bonds <- bond_coords %>% filter(.data$order > 1)

  if (nrow(multi_bonds) > 0) {

    offset_dist <- double_bond_offset

    get_shorten_frac <- function(atom_id, partner_id, bonds, atoms) {
      nb_ids <- bonds %>%
        filter((.data$from == atom_id & .data$to != partner_id) |
                 (.data$to == atom_id & .data$from != partner_id)) %>%
        mutate(nb = ifelse(.data$from == atom_id, .data$to, .data$from)) %>%
        pull(.data$nb) %>% unique()

      n_heavy <- atoms %>%
        filter(.data$atom_id %in% nb_ids, .data$symbol != "H") %>%
        nrow()

      if (n_heavy == 0) return(0.0125)
      return(0.10)
    }

    multi_bonds <- multi_bonds %>%
      rowwise() %>%
      mutate(
        dx  = .data$x2 - .data$x1,
        dy  = .data$y2 - .data$y1,
        len = sqrt(.data$dx * .data$dx + .data$dy * .data$dy),
        ux  = ifelse(.data$len > 1e-8, .data$dx / .data$len, 0),
        uy  = ifelse(.data$len > 1e-8, .data$dy / .data$len, 0),
        px  = -.data$uy,
        py  =  .data$ux,

        side = .mb_choose_double_bond_side(
          from  = .data$from,
          to    = .data$to,
          bonds = bond_coords,
          atoms = atoms
        ),

        off_x = .data$px * offset_dist * .data$side,
        off_y = .data$py * offset_dist * .data$side,

        frac1 = get_shorten_frac(.data$from, .data$to,   bond_coords, atoms),
        frac2 = get_shorten_frac(.data$to,   .data$from, bond_coords, atoms),

        # disk clearance already stored on bond_coords as s1/s2
        # IUPAC cut measured from the ATOM CENTRE
        trim1 = pmax(.data$s1, .data$len * .data$frac1),
        trim2 = pmax(.data$s2, .data$len * .data$frac2),

        ox1 = .data$x1 + .data$ux * .data$trim1 + .data$off_x,
        oy1 = .data$y1 + .data$uy * .data$trim1 + .data$off_y,
        ox2 = .data$x2 - .data$ux * .data$trim2 + .data$off_x,
        oy2 = .data$y2 - .data$uy * .data$trim2 + .data$off_y
      ) %>%
      ungroup()

    p <- p + geom_segment(
      data = multi_bonds,
      aes(x = .data$x1s, y = .data$y1s, xend = .data$x2s, yend = .data$y2s),
      linewidth = bond_width, color = bond_color
    )

    p <- p + geom_segment(
      data = multi_bonds,
      aes(x = .data$ox1, y = .data$oy1, xend = .data$ox2, yend = .data$oy2),
      linewidth = bond_width, color = bond_color
    )

    triples <- multi_bonds %>% filter(.data$order >= 3)
    if (nrow(triples) > 0) {
      triples <- triples %>%
        mutate(
          ox1b = .data$x1s + .data$ux * (.data$len * .data$frac1) - .data$off_x,
          oy1b = .data$y1s + .data$uy * (.data$len * .data$frac1) - .data$off_y,
          ox2b = .data$x2s - .data$ux * (.data$len * .data$frac2) - .data$off_x,
          oy2b = .data$y2s - .data$uy * (.data$len * .data$frac2) - .data$off_y
        )
      p <- p + geom_segment(
        data = triples,
        aes(x = .data$ox1b, y = .data$oy1b, xend = .data$ox2b, yend = .data$oy2b),
        linewidth = bond_width, color = bond_color
      )
    }
  }

  # Atom circles
  atoms_circle <- atoms
  if (hide_carbon_circles) {
    atoms_circle <- atoms %>%
      filter(.data$symbol != "C" | .data$show_label)
  }
  if (show_atom_circles && nrow(atoms_circle) > 0) {
    p <- p + geom_point(data = atoms_circle,
                        aes(x = .data$x, y = .data$y),
                        color = if (circle_stroke <= 0) NA else atoms_circle$color,
                        fill = atom_circle_color,
                        size = atom_size,
                        shape = 21,
                        stroke = circle_stroke)
  }

  # ================= LABELS =================
  if (show_atom_labels) {
    if (!"charge"  %in% names(atoms)) atoms$charge  <- 0L
    if (!"radical" %in% names(atoms)) atoms$radical <- 0L

    atoms_label <- atoms

    if (hide_carbon_labels) {
      atoms_label <- atoms_label %>%
        filter(
          .data$symbol != "C" |
            .data$charge != 0 |
            .data$radical != 0 |
            .data$show_label
        )
    }

    atoms_label <- atoms_label %>%
      mutate(
        charge_txt = dplyr::case_when(
          .data$charge == 0 ~ "",
          abs(.data$charge) == 1 ~ ifelse(.data$charge > 0, "+", "-"),
          TRUE ~ paste0(abs(.data$charge), ifelse(.data$charge > 0, "+", "-"))
        ),
        rad_txt = ifelse(.data$radical %in% c(2, 3), ".", ""),
        super_txt = paste0(.data$rad_txt, .data$charge_txt),
        label_draw = ifelse(
          .data$super_txt == "",
          .data$symbol,
          paste0(.data$symbol, '^{"', .data$super_txt, '"}')
        ),
        label_draw = dplyr::case_when(
          label_fontface == "bold"        ~ paste0("bold(", .data$label_draw, ")"),
          label_fontface == "italic"      ~ paste0("italic(", .data$label_draw, ")"),
          label_fontface == "bold.italic" ~ paste0("bolditalic(", .data$label_draw, ")"),
          TRUE                            ~ .data$label_draw
        )
      )

    p <- p + geom_text(
      data = atoms_label,
      aes(x = .data$x, y = .data$y, label = .data$label_draw),
      parse = TRUE,
      vjust = 0.5,
      size = label_size,
      size.unit = "pt",
      family = label_family,
      fontface = label_fontface,
      color = atoms_label$color
    )

    if (!is.null(result$h_labels) && nrow(result$h_labels) > 0) {
      h_labels <- result$h_labels

      others <- atoms %>%
        filter(.data$symbol != "H") %>%
        select(oid = .data$atom_id, ox = .data$x, oy = .data$y)

      h_off_h  <- H_offset[1]
      h_off_h2 <- if (length(H_offset) >= 2) H_offset[2] else H_offset[1] * 1.35
      h_off_v  <- 0.28 * target_bond_length

      # clash score: smaller = more crowded
      clash <- function(hx, hy, parent_id, others) {
        d <- sqrt((others$ox - hx)^2 + (others$oy - hy)^2)
        d <- d[is.finite(d)]
        if (!length(d)) return(0)
        sum(1 / pmax(d, 0.08)^2)
      }

      placed <- vector("list", nrow(h_labels))
      for (i in seq_len(nrow(h_labels))) {
        row  <- h_labels[i, ]
        base <- if (row$nH == 1) h_off_h else h_off_h2

        cands <- data.frame(
          side  = c("left", "right", "up", "down"),
          hx    = c(row$x - base, row$x + base, row$x, row$x),
          hy    = c(row$y, row$y, row$y + h_off_v, row$y - h_off_v),
          horiz = c(TRUE, TRUE, FALSE, FALSE),
          stringsAsFactors = FALSE
        )

        forced <- if ("side_force" %in% names(h_labels)) row$side_force else NA_character_

        if (!is.na(forced) && nzchar(forced)) {
          best <- cands[cands$side == forced, ]
        } else {
          cands$score <- vapply(seq_len(nrow(cands)), function(j) {
            clash(cands$hx[j], cands$hy[j], row$parent_id, others)
          }, numeric(1))
          cands$score <- cands$score + ifelse(cands$horiz, 0, 1.5)
          best <- cands[which.min(cands$score), ]
        }

        dist_left  <- if (row$nH == 1) base else base * 1.12
        dist_right <- if (row$nH == 1) base else base * 0.82

        hx <- switch(best$side,
                     left  = row$x - dist_left,
                     right = row$x + dist_right,
                     row$x)
        hy <- best$hy

        placed[[i]] <- data.frame(
          parent_id = row$parent_id,
          nH        = row$nH,
          h_text    = row$h_text,
          hx        = hx,
          hy        = hy,
          stringsAsFactors = FALSE
        )
      }
      h_place <- bind_rows(placed)

      parents <- atoms %>% select(parent_id = .data$atom_id, ax = .data$x)

      h_H <- h_place

      h_num <- h_place %>%
        filter(.data$nH > 1) %>%
        mutate(
          lab = as.character(.data$nH),
          hx  = .data$hx + 0.34 * h_off_h2,
          hy  = .data$hy - 0.18 * h_off_h2
        )

      p <- p + geom_text(
        data = h_H,
        aes(x = .data$hx, y = .data$hy, label = "H"),
        parse = FALSE,
        vjust = 0.5,
        hjust = 0.5,
        size = label_size,
        size.unit = "pt",
        family = label_family,
        fontface = label_fontface,
        color = "black"
      )

      if (nrow(h_num) > 0) {
        p <- p + geom_text(
          data = h_num,
          aes(x = .data$hx, y = .data$hy, label = .data$lab),
          parse = FALSE,
          vjust = 1,
          hjust = 0,
          size = label_size * 0.65,
          size.unit = "pt",
          family = label_family,
          fontface = label_fontface,
          color = "black"
        )
      }
    }
  }

  if (show_ids) {
    p <- p + geom_text(data = atoms,
                       aes(x = .data$x + 0.25, y = .data$y + 0.25,
                           label = .data$atom_id),
                       size = 6,
                       size.unit = "pt",
                       color = "red",
                       fontface = "bold")

    if (nrow(bond_coords) > 0) {
      bond_mid <- bond_coords %>%
        mutate(mid_x = (.data$x1 + .data$x2) / 2,
               mid_y = (.data$y1 + .data$y2) / 2,
               bond_id = row_number())

      p <- p + geom_text(data = bond_mid,
                         aes(x = .data$mid_x, y = .data$mid_y,
                             label = .data$bond_id),
                         size = 6,
                         size.unit = "pt",
                         color = "blue",
                         fontface = "bold")
    }
  }

  p <- p +
    coord_fixed(ratio = 1, clip = "off") +
    theme_void() +
    labs(title = title)

  if (nrow(atoms) > 0) {
    x_range <- range(atoms$x, na.rm = TRUE)
    y_range <- range(atoms$y, na.rm = TRUE)

    p <- p +
      xlim(x_range[1] - label_padding, x_range[2] + label_padding) +
      ylim(y_range[1] - label_padding, y_range[2] + label_padding)
  }

  return(p)
}
