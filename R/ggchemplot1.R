#' @title Read a small compound data in SDF format and generate an initial plot
#'
#' @description
#' This function parses a Structure Data File (SDF) and creates a list of data objects
#' specifying atom connectivity together with parameters for plotting and customization.
#'
#' @param sdf_file String. Path to a SDF file containing the chemical structure.
#' @param title String. Optional title for the plot.
#' @param collapse_hydrogens Logical. Whether to collapse hydrogens attached to
#'        heteroatoms (O, N, S, etc.) into labels like OH, NH2, SH.
#'        Default is \code{FALSE} (shows explicit hydrogens). This option needs
#'        to be specified before \code{ggchemplo2}.
#' @param rotation Numeric. Rotation angle in degrees (counterclockwise). This option needs to be specified before \code{ggchemplo2}.
#' @param label_padding Numeric. Amount of padding around the molecule to prevent
#'        labels from being cut off. Higher values add more white space.
#'        Default is 0.6.
#' @param show_atom_circles Logical. Whether to draw circles around atoms.
#' @param hide_carbon_circles Logical. Whether to hide circles around carbon atoms. Useful for showing organic compounds.
#' @param circle_stroke Numeric. Thickness of the atom circle border.
#' @param show_atom_labels Logical. Whether to show atom labels.
#' @param hide_carbon_labels Logical. Whether to hide labels for carbon atoms.
#' @param bond_width Numeric. Thickness of the bonds.
#' @param atom_size Numeric. Size of the atom circles.
#' @param label_size Numeric. Size of the atom labels.
#' @param double_bond_offset Numeric. Offset distance for drawing double and triple bonds.
#' @param custom_atom_colors Named character vector. Custom colors for specific atoms
#'        (e.g. \code{c("N" = "black", "O" = "red")}).
#' @param paint_it_black Logical. If \code{TRUE}, overrides colors to make all atoms black.
#' @param H_offset Numeric vector of length 2. Controls the distance of collapsed
#'        hydrogen labels (first value for single H, second for H2/H3).
#'        Only used when \code{collapse_hydrogens = TRUE}.
#'
#' @return A list containing the processed data and a ggplot object (\code{$plot}).
#'
#' @import ggplot2
#' @import dplyr
#'
#' @export
ggchemplot1 <- function(sdf_file,
                        title = NULL,
                        collapse_hydrogens = FALSE,
                        rotation = 0,
                        label_padding = 0.6,
                        show_atom_circles = TRUE,
                        hide_carbon_circles = TRUE,
                        circle_stroke = 0,
                        show_atom_labels = TRUE,
                        hide_carbon_labels = TRUE,
                        bond_width = 2,
                        atom_size = 20,
                        label_size = 14,
                        double_bond_offset = 0.15,
                        custom_atom_colors = NULL,
                        paint_it_black = FALSE,
                        H_offset = c(0.35, 0.55)) {

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

  # Bonds
  bond_start <- atom_start + n_atoms
  bond_lines <- lines[bond_start:(bond_start + n_bonds - 1)]

  bonds <- data.frame(
    from = as.integer(substr(bond_lines, 1, 3)),
    to   = as.integer(substr(bond_lines, 4, 6)),
    order = as.integer(substr(bond_lines, 7, 9)),
    stringsAsFactors = FALSE
  )

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
  result <- list(atoms = atoms, bond_coords = bond_coords, original_atoms = atoms)

  if (collapse_hydrogens) {
    result <- collapse_hydrogens_func(result)
  }

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
                   H_offset = H_offset)

  # Save parameters for easy tweaking with ggchemplot2
  result$params <- list(
    title = title,
    collapse_hydrogens = collapse_hydrogens,
    rotation = rotation,
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
    H_offset = H_offset
  )

  result$plot <- p
  return(result)
}
