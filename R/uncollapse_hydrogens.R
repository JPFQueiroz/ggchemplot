#' Uncollapse (restore explicit) hydrogens bound to one or more atoms
#'
#' @param result A list returned by \code{ggchemplot1} (after collapse).
#' @param atom_id Integer vector of atom_id(s) whose attached hydrogens
#'        should be shown explicitly.
#'
#' @return The modified result list (ready for \code{ggchemplot2}).
#' @export
uncollapse_hydrogens <- function(result, atom_id) {
  if (is.null(result$original_atoms) || is.null(result$original_bond_coords)) {
    stop("original_atoms / original_bond_coords missing. ")
  }

  atom_id <- as.integer(atom_id)

  # Bonds that involve a hydrogen attached to the requested parent(s)
  h_bonds <- result$original_bond_coords %>%
    filter(
      (.data$from %in% atom_id & .data$sym2 == "H") |
        (.data$to   %in% atom_id & .data$sym1 == "H")
    )

  if (nrow(h_bonds) == 0) {
    message("No hydrogens attached to atom_id(s): ", paste(atom_id, collapse = ", "))
    return(result)
  }

  # The hydrogen atom_ids themselves
  h_ids <- unique(c(
    h_bonds$from[h_bonds$sym1 == "H"],
    h_bonds$to  [h_bonds$sym2 == "H"]
  ))

  # Restore the H atoms (match the columns that currently exist in result$atoms)
  h_atoms <- result$original_atoms %>%
    filter(.data$atom_id %in% h_ids)

  # Align columns with the (possibly collapsed) atoms table
  extra_cols <- setdiff(names(result$atoms), names(h_atoms))
  for (col in extra_cols) {
    h_atoms[[col]] <- if (col == "display_symbol") h_atoms$symbol else NA
  }
  h_atoms <- h_atoms[, names(result$atoms), drop = FALSE]

  # Avoid duplicates if the user calls the function twice
  already_present <- result$atoms$atom_id %in% h_ids
  result$atoms <- bind_rows(
    result$atoms[!already_present, , drop = FALSE],
    h_atoms
  )

  # Restore the H-bonds (keep any extra columns that may exist)
  common_cols <- intersect(names(result$bond_coords), names(h_bonds))
  h_bonds <- h_bonds[, common_cols, drop = FALSE]

  # Drop any already-restored H-bonds
  key_exist <- paste(result$bond_coords$from, result$bond_coords$to, sep = "-")
  key_new   <- paste(h_bonds$from, h_bonds$to, sep = "-")
  h_bonds   <- h_bonds[!key_new %in% key_exist, , drop = FALSE]

  result$bond_coords <- bind_rows(result$bond_coords, h_bonds)

  # Remove the collapsed H-label for these parents (otherwise you get both
  # the explicit H atoms and the old "OH"/"NH2" text)
  if (!is.null(result$h_labels) && nrow(result$h_labels) > 0) {
    result$h_labels <- result$h_labels %>%
      filter(!.data$parent_id %in% atom_id)
  }

  # Optional but safe: rebuild coordinates from the current atoms table
  # (useful if any later helper moved atoms)
  if (exists("rebuild_bond_coords", mode = "function")) {
    result <- rebuild_bond_coords(result)
  }

  result$plot <- ggchemplot2(result)
  result
}
