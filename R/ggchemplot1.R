#' @title Read a small compound data in SDF format and generate an initial plot
#'
#' @description
#' Parses a Structure Data File (SDF, V2000) and returns a list of atom and
#' bond tables plus a ggplot object. Coordinates can be rotated, flipped,
#' and normalized. Stereochemical bonds, charges and radicals in the file
#' are read when present. Hydrogens on heteroatoms can be collapsed into
#' \code{OH} / \code{NH2} / \code{SH} labels.
#'
#' @param sdf_file Character. Path to an SDF file, or the file contents as
#'   a character vector of lines.
#' @param use_sdf_stereo Logical. If \code{TRUE} (default), V2000 bond stereo
#'   codes are mapped to wedge (\code{stereo = 1}) and hashed
#'   (\code{stereo = 6}) bonds. If \code{FALSE}, all bonds are drawn as
#'   plain lines.
#' @param title Character. Optional plot title.
#' @param collapse_hydrogens Logical. If \code{TRUE}, hydrogens attached to
#'   heteroatoms are removed and drawn as collapsed labels (\code{OH},
#'   \code{NH2}, \code{SH}). Must be set here; it cannot be changed later
#'   in \code{ggchemplot2()}. Default \code{FALSE}.
#' @param rotation Numeric. Rotation in degrees, counterclockwise. Applied
#'   before flipping. Must be set here, not in \code{ggchemplot2()}.
#'   Default \code{0}.
#' @param flip_horizontal Logical. If \code{TRUE}, reflect the structure
#'   through a vertical axis (\code{x -> -x}). A single-axis flip inverts
#'   wedge and hashed bonds so that the absolute configuration is kept.
#'   Flipping both axes is a 180 degrees in-plane rotation and does not invert
#'   stereo. Default \code{FALSE}.
#' @param flip_vertical Logical. If \code{TRUE}, reflect the structure
#'   through a horizontal axis (\code{y -> -y}). See
#'   \code{flip_horizontal}. Default \code{FALSE}.
#' @param label_padding Numeric. Extra data-unit margin around the molecule
#'   so labels are not clipped. Default \code{0.6} when \code{NULL} is
#'   resolved in \code{ggchemplot2()}.
#' @param show_atom_circles Logical. Draw filled discs behind atom labels.
#'   Default \code{TRUE}.
#' @param hide_carbon_circles Logical. If \code{TRUE}, carbon discs are
#'   omitted (skeletal drawing). Carbons that carry a charge or a radical
#'   are still labelled. Default \code{TRUE}.
#' @param circle_stroke Numeric. Border width of atom discs. Default \code{0}.
#' @param show_atom_labels Logical. Draw element symbols. Default \code{TRUE}.
#' @param hide_carbon_labels Logical. If \code{TRUE}, carbon symbols are
#'   omitted unless the atom has a charge or a radical. Default \code{TRUE}.
#' @param bond_width Numeric. Bond line width in mm. If \code{NULL}, derived
#'   from \code{label_size}.
#' @param atom_size Numeric. Atom disc size in mm. If \code{NULL}, derived
#'   from \code{label_size}.
#' @param label_size Numeric. Atom label size in points. If \code{NULL},
#'   \code{12} is used in \code{ggchemplot2()}.
#' @param double_bond_offset Numeric. Perpendicular offset of the second
#'   (and third) stroke of a multiple bond, in data units. If \code{NULL},
#'   \code{0.15} is used.
#' @param custom_atom_colors Named character vector of colours, e.g.
#'   \code{c("N" = "black", "O" = "red")}. Unlisted elements keep the
#'   default palette.
#' @param paint_it_black Logical. If \code{TRUE}, every atom is drawn black.
#'   Default \code{FALSE}.
#' @param bond_color Character. Colour of bonds. Default \code{"black"}.
#' @param atom_circle_color Character. Fill of atom discs. Default
#'   \code{"transparent"}.
#' @param H_offset Numeric. Collapsed-hydrogen gap in data units.
#'   Length 1: same horizontal gap for H and H2; vertical = 1.15 times that.
#'   Length 2: horizontal for H, horizontal for H2; vertical defaults to 1.15 times each.
#'   Length 4: horizontal H, horizontal H2, vertical H, vertical H2.
#'   If \code{NULL}, all four are computed from \code{label_size}.
#' @param label_fontface Character. One of \code{"plain"}, \code{"bold"},
#'   \code{"italic"}, \code{"bold.italic"}. Default \code{"plain"}.
#' @param normalize Logical. If \code{TRUE}, scale coordinates so the median
#'   bond length equals \code{target_bond_length}. Default \code{TRUE}.
#' @param target_bond_length Numeric. Target median bond length in data units
#'   when \code{normalize = TRUE}. Default \code{1}.
#'
#' @return A list containing the processed data and a ggplot object (\code{$plot}).
#'
#' @import ggplot2
#' @import dplyr
#'
#' @importFrom stats setNames
#'
#' @export
ggchemplot1 <- function(sdf_file,
                        use_sdf_stereo = TRUE,
                        title = NULL,
                        collapse_hydrogens = FALSE,
                        rotation = 0,
                        flip_horizontal = FALSE,
                        flip_vertical   = FALSE,
                        label_padding = NULL,
                        show_atom_circles = TRUE,
                        hide_carbon_circles = TRUE,
                        circle_stroke = 0,
                        show_atom_labels = TRUE,
                        hide_carbon_labels = TRUE,
                        bond_width = NULL,
                        atom_size = NULL,
                        label_size = NULL,
                        double_bond_offset = NULL,
                        custom_atom_colors = NULL,
                        paint_it_black = FALSE,
                        bond_color = "black",
                        atom_circle_color = "transparent",
                        H_offset = NULL,
                        label_fontface = "plain",
                        normalize = TRUE,
                        target_bond_length = 1.0) {

  # Read SDF
  lines <-
    if (file.exists(sdf_file)) readLines(sdf_file, warn = FALSE) else
      sdf_file
  counts_line <- grep("^\\s*\\d+\\s+\\d+", lines)[1]
  if (is.na(counts_line)) stop("Invalid SDF format")

  counts <- strsplit(trimws(lines[counts_line]), "\\s+")[[1]]
  n_atoms <- as.integer(counts[1])
  n_bonds <- as.integer(counts[2])

  # Atoms
  atom_start <- counts_line + 1
  atom_lines <- lines[atom_start:(atom_start + n_atoms - 1)]

  atoms <- data.frame(
    atom_id = 1:n_atoms,
    x = as.numeric(substr(atom_lines, 1, 10)),
    y = as.numeric(substr(atom_lines, 11, 20)),
    symbol = trimws(substr(atom_lines, 32, 34)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(default_color = case_when(
      .data$symbol == "C"  ~ "black",
      .data$symbol == "O"  ~ "#e31a1c",
      .data$symbol == "N"  ~ "#1f78b4",
      .data$symbol == "P"  ~ "violet",
      .data$symbol == "S"  ~ "#e6ab02",
      .data$symbol == "H"  ~ "gray40",
      TRUE ~ "black"
    )) %>%
    mutate(default_color = if (paint_it_black) "black" else .data$default_color)

  if (!is.null(custom_atom_colors)) {
    atoms$color <- ifelse(atoms$symbol %in% names(custom_atom_colors),
                          custom_atom_colors[atoms$symbol],
                          atoms$default_color)
  } else {
    atoms$color <- atoms$default_color
  }

  # legacy atom-block charge/radical
  legacy_code <- suppressWarnings(as.integer(trimws(substr(atom_lines, 37, 39))))
  legacy_code[is.na(legacy_code)] <- 0L

  atoms$charge  <- dplyr::case_when(
    legacy_code == 1L ~  3L,
    legacy_code == 2L ~  2L,
    legacy_code == 3L ~  1L,
    legacy_code == 5L ~ -1L,
    legacy_code == 6L ~ -2L,
    legacy_code == 7L ~ -3L,
    TRUE              ~  0L
  )
  atoms$radical <- ifelse(legacy_code == 4L, 2L, 0L)

  # ====================== APPLY ROTATION ======================
  if (rotation != 0) {
    theta <- rotation * pi / 180
    rot_matrix <- matrix(c(cos(theta), -sin(theta),
                           sin(theta),  cos(theta)), nrow = 2)

    coords <- as.matrix(atoms[, c("x", "y")])
    rotated <- coords %*% t(rot_matrix)

    atoms$x <- rotated[, 1]
    atoms$y <- rotated[, 2]
  }

  # ====================== FLIP (keep absolute stereo) ======================
  if (isTRUE(flip_horizontal)) atoms$x <- -atoms$x
  if (isTRUE(flip_vertical))   atoms$y <- -atoms$y

  invert_stereo <- xor(isTRUE(flip_horizontal), isTRUE(flip_vertical))

  # Bonds
  bond_start <- atom_start + n_atoms
  bond_lines <- lines[bond_start:(bond_start + n_bonds - 1)]

  # ----- property block -----
  prop_start <- bond_start + n_bonds
  prop_lines <- if (prop_start <= length(lines)) {
    lines[prop_start:length(lines)]
  } else {
    character(0)
  }
  chg_lines <- grep("^M  CHG", prop_lines, value = TRUE)
  rad_lines <- grep("^M  RAD", prop_lines, value = TRUE)

  parse_m_pairs <- function(plines, tag) {
    if (length(plines) == 0) return(NULL)
    out <- list()
    for (ln in plines) {
      tok <- strsplit(trimws(ln), "\\s+")[[1]]
      if (length(tok) < 5) next
      n <- suppressWarnings(as.integer(tok[3]))
      if (is.na(n) || n < 1) next
      for (i in seq_len(n)) {
        aid <- as.integer(tok[2 + 2 * i])
        val <- as.integer(tok[3 + 2 * i])
        if (!is.na(aid) && !is.na(val)) out[[length(out) + 1]] <- c(aid, val)
      }
    }
    if (!length(out)) return(NULL)
    as.data.frame(do.call(rbind, out), stringsAsFactors = FALSE) |>
      setNames(c("atom_id", tag))
  }

  if (length(chg_lines) || length(rad_lines)) {
    atoms$charge  <- 0L
    atoms$radical <- 0L
  }

  chg_df <- parse_m_pairs(chg_lines, "charge")
  if (!is.null(chg_df)) {
    idx <- match(chg_df$atom_id, atoms$atom_id)
    ok  <- !is.na(idx)
    atoms$charge[idx[ok]] <- chg_df$charge[ok]
  }

  rad_df <- parse_m_pairs(rad_lines, "radical")
  if (!is.null(rad_df)) {
    idx <- match(rad_df$atom_id, atoms$atom_id)
    ok  <- !is.na(idx)
    atoms$radical[idx[ok]] <- rad_df$radical[ok]
  }

  bonds <- data.frame(
    from  = as.integer(substr(bond_lines, 1, 3)),
    to    = as.integer(substr(bond_lines, 4, 6)),
    order = as.integer(substr(bond_lines, 7, 9)),
    stereo = {
      raw <- suppressWarnings(as.integer(substr(bond_lines, 10, 12)))
      ifelse(is.na(raw), 0L, raw)
    },
    stringsAsFactors = FALSE
  )

  # map V2000 stereo
  bonds$bond_type <- dplyr::case_when(
    bonds$order == 1 & bonds$stereo == 1L ~ "solid",   # wedge
    bonds$order == 1 & bonds$stereo == 6L ~ "hashed",  # hash
    TRUE ~ NA_character_
  )

  # defaults used by geom_wedge
  bonds$shorten_start    <- ifelse(!is.na(bonds$bond_type), 0.15, NA_real_)
  bonds$shorten_end      <- ifelse(!is.na(bonds$bond_type), 0.22, NA_real_)
  bonds$width            <- ifelse(bonds$bond_type %in% "hashed", 0.28,
                                   ifelse(bonds$bond_type %in% "solid", 0.26, NA_real_))
  bonds$wedge_thickness  <- ifelse(bonds$bond_type %in% "hashed", 0.10, NA_real_)
  bonds$n_hashes         <- ifelse(bonds$bond_type %in% "hashed", 8L, NA_integer_)
  bonds$colour           <- ifelse(!is.na(bonds$bond_type), "black", NA_character_)

  # Ignore stereo option
  if (!isTRUE(use_sdf_stereo)) {
    bonds$bond_type        <- NA_character_
    bonds$shorten_start    <- NA_real_
    bonds$shorten_end      <- NA_real_
    bonds$width            <- NA_real_
    bonds$wedge_thickness  <- NA_real_
    bonds$n_hashes         <- NA_integer_
    bonds$colour           <- NA_character_
  }

  if (isTRUE(invert_stereo) && isTRUE(use_sdf_stereo)) {
    bonds$bond_type <- dplyr::case_when(
      bonds$bond_type == "solid"  ~ "hashed",
      bonds$bond_type == "hashed" ~ "solid",
      TRUE ~ bonds$bond_type
    )
    bonds$stereo <- dplyr::case_when(
      bonds$stereo == 1L ~ 6L,
      bonds$stereo == 6L ~ 1L,
      TRUE ~ bonds$stereo
    )
    # keep hashed / solid visual defaults consistent after the swap
    bonds$width <- ifelse(bonds$bond_type %in% "hashed", 0.28,
                          ifelse(bonds$bond_type %in% "solid", 0.26, bonds$width))
    bonds$wedge_thickness <- ifelse(bonds$bond_type %in% "hashed", 0.10, NA_real_)
    bonds$n_hashes <- ifelse(bonds$bond_type %in% "hashed", 8L, NA_integer_)
  }

  bond_coords <- bonds %>%
    left_join(atoms %>% select(.data$atom_id, .data$x, .data$y, .data$symbol,
                               .data$color),
              by = c("from" = "atom_id")) %>%
    rename(x1 = .data$x, y1 = .data$y, sym1 = .data$symbol,
           col1 = .data$color) %>%
    left_join(atoms %>% select(.data$atom_id, .data$x, .data$y, .data$symbol),
              by = c("to" = "atom_id")) %>%
    rename(x2 = .data$x, y2 = .data$y, sym2 = .data$symbol)

  # ====================== COLLAPSE HYDROGENS ======================
  result <- list(atoms = atoms,
                 bond_coords = bond_coords,
                 original_atoms = atoms,
                 original_bond_coords = bond_coords
                 )

  if (collapse_hydrogens) {
    result <- collapse_hydrogens_func(result)
  }

  # Save parameters for easy tweaking with ggchemplot2
  result$params <- list(
    title = title,
    collapse_hydrogens = collapse_hydrogens,
    rotation = rotation,
    flip_horizontal = flip_horizontal,
    flip_vertical   = flip_vertical,
    label_padding = label_padding,
    show_atom_circles = show_atom_circles,
    hide_carbon_circles = hide_carbon_circles,
    circle_stroke = circle_stroke,
    show_atom_labels = show_atom_labels,
    hide_carbon_labels = hide_carbon_labels,
    bond_width = bond_width,
    atom_size = atom_size,
    label_size = label_size,
    double_bond_offset = double_bond_offset,
    custom_atom_colors = custom_atom_colors,
    paint_it_black = paint_it_black,
    H_offset = H_offset,
    label_fontface = label_fontface,
    target_bond_length = target_bond_length,
    normalize = normalize
  )

  # Final plot using ggchemplot2
  p <- ggchemplot2(result,
                   title = title,
                   label_padding = label_padding,
                   show_atom_circles = show_atom_circles,
                   hide_carbon_circles = hide_carbon_circles,
                   circle_stroke = circle_stroke,
                   show_atom_labels = show_atom_labels,
                   hide_carbon_labels = hide_carbon_labels,
                   bond_width = bond_width,
                   atom_size = atom_size,
                   label_size = label_size,
                   double_bond_offset = double_bond_offset,
                   paint_it_black = paint_it_black,
                   bond_color = bond_color,
                   atom_circle_color = atom_circle_color,
                   H_offset = H_offset,
                   label_fontface = label_fontface,
                   target_bond_length = target_bond_length,
                   normalize = normalize)


  result$plot <- p
  return(result)
}
