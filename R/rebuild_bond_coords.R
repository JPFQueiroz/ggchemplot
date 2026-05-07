#' @import ggplot2
#' @import dplyr
#'
#' @keywords internal
# Helper to rebuild coordinates
rebuild_bond_coords <- function(result) {
  if (is.null(result$atoms) || is.null(result$bond_coords)) return(result)

  atoms_temp <- result$atoms %>% select(.data$atom_id, .data$x, .data$y,
                                        .data$symbol)

  result$bond_coords <- result$bond_coords %>%
    select(.data$from, .data$to, .data$order, any_of("bond_type"),
           any_of("shorten")) %>%
    left_join(atoms_temp, by = c("from" = "atom_id")) %>%
    rename(x1 = .data$x, y1 = .data$y, sym1 = .data$symbol) %>%
    left_join(atoms_temp, by = c("to" = "atom_id")) %>%
    rename(x2 = .data$x, y2 = .data$y, sym2 = .data$symbol)

  return(result)
}
