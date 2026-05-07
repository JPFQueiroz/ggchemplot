#' @import ggplot2
#' @import dplyr
#'
#' @keywords internal
collapse_hydrogens_func <- function(res) {
  atoms <- res$atoms
  bonds <- res$bond_coords

  if (!any(atoms$symbol == "H")) return(res)

  hetero <- c("O", "N", "S", "P")

  h_bonds <- bonds %>%
    filter(.data$sym1 == "H" | .data$sym2 == "H") %>%
    mutate(parent_id = ifelse(.data$sym1 == "H", .data$to, .data$from),
           parent_sym = ifelse(.data$sym1 == "H", .data$sym2, .data$sym1))

  # Only heteroatoms get collapsed H
  h_counts <- h_bonds %>%
    filter(.data$parent_sym %in% hetero) %>%
    count(.data$parent_id, name = "nH")

  # Heavy atoms stay exactly at bond attachment point
  heavy_atoms <- atoms %>%
    left_join(h_counts, by = c("atom_id" = "parent_id")) %>%
    mutate(display_symbol = .data$symbol) %>%
    filter(.data$symbol != "H")

  # Separate H labels
  h_labels <- h_counts %>%
    left_join(atoms %>% select(.data$atom_id, .data$x, .data$y),
              by = c("parent_id" = "atom_id")) %>%
    mutate(
      h_text = case_when(
        nH == 1 ~ "H",
        TRUE ~ paste0("H[", nH, "]")
      )
    )

  res$atoms <- heavy_atoms
  res$h_labels <- h_labels
  res$bond_coords <- bonds %>% filter(!(.data$sym1 == "H" | .data$sym2 == "H"))

  return(res)
}
