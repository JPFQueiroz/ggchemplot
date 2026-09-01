#' @importFrom stats median
#'
#' @import dplyr
#'
#' @keywords internal
normalize_structure <- function(result, target_bond_length = target_bond_length) {

  if (is.null(result$atoms) || nrow(result$atoms) == 0) {
    return(result)
  }

  bond_coords <- result$bond_coords

  if (is.null(bond_coords) || nrow(bond_coords) == 0) {
    warning("Normalization skipped: No bonds found.")
    return(result)
  }

  # Calculate bond lengths
  bond_coords <- bond_coords %>%
    mutate(
      length = sqrt((.data$x2 - .data$x1)^2 + (.data$y2 - .data$y1)^2)
    )

  median_length <- median(bond_coords$length, na.rm = TRUE)

  if (is.na(median_length) || median_length < 1e-8) {
    warning("Normalization skipped: Invalid bond lengths.")
    return(result)
  }

  scale_factor <- target_bond_length / median_length

  # Debug message (helpful during development)
  message(sprintf("Normalizing: median bond length = %.3f, target = %.2f (scale = %.3f)",
                  median_length, target_bond_length, scale_factor))

  # Scale atoms
  result$atoms <- result$atoms %>%
    mutate(x = .data$x * scale_factor,
           y = .data$y * scale_factor)

  # Scale bonds
  result$bond_coords <- result$bond_coords %>%
    mutate(x1 = .data$x1 * scale_factor,
           y1 = .data$y1 * scale_factor,
           x2 = .data$x2 * scale_factor,
           y2 = .data$y2 * scale_factor)

  # Scale H labels
  if (!is.null(result$h_labels) && nrow(result$h_labels) > 0) {
    result$h_labels <- result$h_labels %>%
      mutate(x = .data$x * scale_factor,
             y = .data$y * scale_factor)
  }

  # Scale protein links
  if (!is.null(result$protein_links) && length(result$protein_links) > 0) {
    result$protein_links <- lapply(result$protein_links, function(link) {
      link$data <- link$data %>%
        mutate(x = .data$x * scale_factor,
               y = .data$y * scale_factor)
      link
    })
  }

  # Store metadata
  result$params$scale_factor <- scale_factor
  result$params$target_bond_length <- target_bond_length
  result$params$was_normalized <- TRUE
  result$params$original_median_bond <- median_length

  return(result)
}
