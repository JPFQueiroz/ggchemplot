#' @title Add a metal to a sketched structure
#'
#' @param result A ggchemplot result.
#' @param atom_id Atom that binds the metal.
#' @param metal Element symbol, such as \code{"Ni"} or \code{"Fe"}.
#' @param oxidation Oxidation state. Stored as the charge. Default \code{2}.
#' @param bond_length Bond length. Default: median heavy-atom bond.
#' @return The updated result.
#' @export
add_metal <- function(result, atom_id, metal, oxidation = 2L, bond_length = NULL) {
  atoms <- result$atoms
  parent <- atoms[match(atom_id, atoms$atom_id), ]
  if (!nrow(parent)) stop("Atom ", atom_id, " not found.")
  if (.sk_free(result, atom_id) < 1L) {
    dbl <- result$bond_coords$order >= 2 &
      (result$bond_coords$from == atom_id | result$bond_coords$to == atom_id)
    if (any(dbl)) {
      result$bond_coords$order[which(dbl)[1]] <-
        result$bond_coords$order[which(dbl)[1]] - 1L
    }
  }
  len <- bond_length %||% .sk_length(result, "C")
  dir <- .sk_direction(result, atom_id, NULL)
  id <- max(atoms$atom_id) + 1L
  metal_row <- .sk_atom(
    id,
    parent$x + dir[1] * len,
    parent$y + dir[2] * len,
    metal,
    as.integer(oxidation)
  )
  metal_row$show_label <- TRUE
  result$atoms <- rbind(result$atoms, metal_row)
  result$bond_coords <- rbind(
    result$bond_coords[, .sk_bond_cols()],
    .sk_bond(atom_id, id, 1L)
  )
  .sk_finish(result$atoms, result$bond_coords, len, result$params)
}
