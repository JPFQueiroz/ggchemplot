#' @title Rotate a fragment around an atom
#'
#' @description
#' Rotates every atom on one side of a bond around a pivot atom.
#' The pivot stays fixed. The bond is given by the pivot and the first
#' atom of the fragment that should move.
#'
#' @param result A ggchemplot result.
#' @param pivot_id Atom id used as the centre of rotation.
#' @param through_id Atom id of the neighbour that defines the side to rotate.
#' @param angle Rotation in degrees, counterclockwise.
#' @return The updated result.
#' @export
rotate_fragment <- function(result, pivot_id, through_id, angle) {
  atoms <- result$atoms
  bonds <- result$bond_coords
  pivot <- atoms[match(pivot_id, atoms$atom_id), ]
  if (!nrow(pivot)) stop("Pivot atom ", pivot_id, " not found.")
  if (!any((bonds$from == pivot_id & bonds$to == through_id) |
           (bonds$from == through_id & bonds$to == pivot_id))) {
    stop("No bond between ", pivot_id, " and ", through_id, ".")
  }

  # atoms on the moving side (BFS that does not cross back to the pivot)
  adj <- split(c(bonds$to, bonds$from), c(bonds$from, bonds$to))
  seen <- integer()
  stack <- through_id
  while (length(stack)) {
    id <- stack[[1]]
    stack <- stack[-1]
    if (id %in% seen || id == pivot_id) next
    seen <- c(seen, id)
    nbs <- adj[[as.character(id)]]
    if (!is.null(nbs)) stack <- c(stack, nbs[!nbs %in% c(seen, pivot_id)])
  }
  if (!length(seen)) return(result)

  theta <- angle * pi / 180
  rot <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2)
  hit <- atoms$atom_id %in% seen
  xy <- as.matrix(atoms[hit, c("x", "y")])
  xy <- sweep(xy, 2, c(pivot$x, pivot$y), "-")
  xy <- xy %*% t(rot)
  atoms$x[hit] <- xy[, 1] + pivot$x
  atoms$y[hit] <- xy[, 2] + pivot$y
  result$atoms <- atoms

  # rebuild geometry columns if present
  if (all(c("x1", "y1", "x2", "y2") %in% names(bonds))) {
    bonds$x1 <- atoms$x[match(bonds$from, atoms$atom_id)]
    bonds$y1 <- atoms$y[match(bonds$from, atoms$atom_id)]
    bonds$x2 <- atoms$x[match(bonds$to, atoms$atom_id)]
    bonds$y2 <- atoms$y[match(bonds$to, atoms$atom_id)]
    result$bond_coords <- bonds
  }

  # collapsed-H labels that belong to a moved atom
  if (!is.null(result$h_labels) && nrow(result$h_labels) &&
      "parent_id" %in% names(result$h_labels)) {
    h_hit <- result$h_labels$parent_id %in% seen
    if (any(h_hit)) {
      hxy <- as.matrix(result$h_labels[h_hit, c("x", "y")])
      hxy <- sweep(hxy, 2, c(pivot$x, pivot$y), "-")
      hxy <- hxy %*% t(rot)
      result$h_labels$x[h_hit] <- hxy[, 1] + pivot$x
      result$h_labels$y[h_hit] <- hxy[, 2] + pivot$y
    }
  }

  if (!is.null(result$plot)) result$plot <- ggchemplot2(result)
  result
}
