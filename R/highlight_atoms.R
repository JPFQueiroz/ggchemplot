#' Highlight atoms with a translucent disc
#'
#' @param result Object from \code{ggchemplot1()}.
#' @param atom_id Integer vector of atom ids.
#' @param alpha Fill transparency. Default \code{0.35}.
#' @param size Disc size in mm. Default \code{10}.
#' @param fill Fill colour. Default \code{"gold"}.
#' @param colour Stroke colour. Default \code{NA}.
#' @param stroke Stroke width in mm. Default \code{0}.
#'
#' @return The same \code{result} list, with \code{$highlights} updated
#'   and \code{$plot} rebuilt.
#' @export
highlight_atoms <- function(result, atom_id, alpha = 0.35, size = 10,
                            fill = "gold", colour = NA, stroke = 0) {
  atom_id <- as.integer(atom_id)
  new <- data.frame(
    atom_id = atom_id,
    alpha = alpha, size = size,
    fill = as.character(fill),
    colour = if (length(colour)) as.character(colour) else NA_character_,
    stroke = stroke,
    stringsAsFactors = FALSE
  )
  result$highlights <- if (is.null(result$highlights)) new else rbind(result$highlights, new)
  result$plot <- ggchemplot2(result)   # respects result$params$normalize
  result
}
