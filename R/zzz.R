#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to ggchemplot - Publication-ready 2D chemical structures with ggplot2")
}

.onLoad <- function(libname, pkgname) {
  # This runs earlier than .onAttach
  invisible()
}
