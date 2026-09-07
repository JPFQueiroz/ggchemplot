#' @title Show the label of selected atoms
#' @param result Object from \code{ggchemplot1}.
#' @param atom_id Integer vector of atom ids.
#' @export
show_atom_label <- function(result, atom_id) {
  if (!"show_label" %in% names(result$atoms)) {
    result$atoms$show_label <- FALSE
  }
  result$atoms$show_label[result$atoms$atom_id %in% atom_id] <- TRUE
  result$plot <- ggchemplot2(result)
  result
}
