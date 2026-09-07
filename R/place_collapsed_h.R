#' @title Place a collapsed hydrogen
#'
#' @param result Object from \code{ggchemplot1}.
#' @param atom_id Parent heteroatom id.
#' @param side One of \code{"left"}, \code{"right"}, \code{"up"}, \code{"down"}.
#'
#' @export
place_collapsed_h <- function(result, atom_id, side = c("left", "right", "up", "down")) {
  side <- match.arg(side)
  if (is.null(result$h_labels) || nrow(result$h_labels) == 0) {
    stop("No collapsed hydrogens. Set collapse_hydrogens = TRUE in ggchemplot1().")
  }
  if (!"side_force" %in% names(result$h_labels)) {
    result$h_labels$side_force <- NA_character_
  }
  hit <- result$h_labels$parent_id == atom_id
  if (!any(hit)) stop("No collapsed H on atom_id ", atom_id)
  result$h_labels$side_force[hit] <- side
  result$plot <- ggchemplot2(result)
  result
}
