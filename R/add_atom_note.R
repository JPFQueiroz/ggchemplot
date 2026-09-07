#' @title Annotate an atom with subscript text
#'
#' @param result Object from \code{ggchemplot1}.
#' @param atom_id Integer. Atom to hang the note on.
#' @param label Character. Text (e.g. \code{"1"}, \code{"R"}, \code{"S"}).
#' @param distance Numeric. Offset from the atom, in data units.
#'   Default \code{0.28}.
#' @param direction Numeric. Direction of the offset in degrees:
#'   \code{0} = right, \code{90} = up, \code{180} = left, \code{270} = down.
#' @param angle Numeric. Rotation of the letters, in degrees
#'   counterclockwise. Default \code{0}.
#'
#' @export
add_atom_note <- function(result, atom_id, label,
                          distance = 0.28,
                          direction = 0,
                          angle = 0) {
  atoms <- result$atoms
  hit <- atoms$atom_id == atom_id
  if (!any(hit)) stop("No atom_id ", atom_id)

  note <- data.frame(
    atom_id   = as.integer(atom_id),
    label     = as.character(label),
    distance  = distance,
    direction = direction,
    angle     = angle,
    stringsAsFactors = FALSE
  )

  if (is.null(result$atom_notes) || nrow(result$atom_notes) == 0) {
    result$atom_notes <- note
  } else {
    result$atom_notes <- rbind(result$atom_notes, note)
  }

  result$plot <- ggchemplot2(result)
  result
}
