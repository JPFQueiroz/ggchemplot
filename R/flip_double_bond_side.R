#' @title Swap the offset side of a double bond
#'
#' @param result Object from \code{ggchemplot1}.
#' @param from_id,to_id Atom ids of the bond (order does not matter).
#' @param side \code{+1} or \code{-1}. If \code{NULL}, flip whatever
#'   \code{ggchemplot2()} would choose.
#'
#' @export
flip_double_bond_side <- function(result, from_id, to_id, side = NULL) {
  bc <- result$bond_coords
  hit <- (bc$from == from_id & bc$to == to_id) |
    (bc$from == to_id   & bc$to == from_id)

  if (!any(hit)) stop("No bond between ", from_id, " and ", to_id)
  if (!any(bc$order[hit] == 2)) {
    stop("That bond is not a double bond")
  }

  if (!"mb_side" %in% names(bc)) bc$mb_side <- NA_integer_

  if (is.null(side)) {
    cur <- bc$mb_side[hit]
    bc$mb_side[hit] <- ifelse(is.na(cur) | cur != -1L, -1L, 1L)
  } else {
    side <- as.integer(sign(side))
    if (!side %in% c(-1L, 1L)) stop("side must be +1 or -1")
    bc$mb_side[hit] <- side
  }

  result$bond_coords <- bc
  result$plot <- ggchemplot2(result)
  result
}
