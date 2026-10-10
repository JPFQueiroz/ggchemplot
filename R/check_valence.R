#' @title Check valence of a ggchemplot result
#'
#' @description
#' Compares the bond-order sum at each atom with the valence expected from
#' the element, charge and radical. Uses \code{original_atoms} and
#' \code{original_bond_coords} when present, so collapsed carbon hydrogens
#' are not lost. Metals are reported but not flagged.
#'
#' @param result A ggchemplot or sketch result.
#' @param only_problems Logical. If \code{TRUE}, return only atoms that
#'   do not match. Default \code{FALSE}.
#' @return A data frame with one row per atom.
#' @export
check_valence <- function(result, only_problems = FALSE) {
  atoms <- result$original_atoms %||% result$atoms
  bonds <- result$original_bond_coords %||% result$bond_coords
  metals <- c("Fe", "Ni", "Co", "Cu", "Zn", "Mg", "Mn", "Ca", "Na", "K")

  expected_bonds <- function(symbol, charge) {
    base <- c(C = 4L, N = 3L, O = 2L, S = 2L, P = 5L, H = 1L,
              F = 1L, Cl = 1L, Br = 1L, I = 1L)[symbol]
    if (is.na(base)) return(NA_integer_)
    charge <- charge %||% 0L
    if (symbol == "C") as.integer(base - abs(charge)) else as.integer(base + charge)
  }

  heavy <- atoms[atoms$symbol != "H", ]
  out <- lapply(seq_len(nrow(heavy)), function(i) {
    a <- heavy[i, ]
    hit <- bonds$from == a$atom_id | bonds$to == a$atom_id
    used <- if (!any(hit)) 0L else sum(bonds$order[hit], na.rm = TRUE)
    charge <- a$charge %||% 0L
    radical <- a$radical %||% 0L
    target <- expected_bonds(a$symbol, charge)
    if (!is.na(target)) target <- target - radical
    problem <- !is.na(target) && used != target && !(a$symbol %in% metals)
    tag <- if (charge == 0) "" else if (charge > 0) paste0(charge, "+") else paste0(abs(charge), "-")
    data.frame(
      atom_id = a$atom_id,
      symbol = a$symbol,
      charge = charge,
      radical = radical,
      bond_order_sum = used,
      expected = target,
      problem = problem,
      message = if (!problem) "ok" else {
        sprintf("%s%s has %d bonds, expected %d", a$symbol, tag, used, target)
      },
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, out)
  if (only_problems) tab <- tab[tab$problem, , drop = FALSE]
  rownames(tab) <- NULL
  tab
}
