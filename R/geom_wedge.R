#' @import ggplot2
#' @import grid
#'
#' @keywords internal
GeomWedge <-
  ggplot2::ggproto("GeomWedge", ggplot2::Geom,
                   required_aes = c("x1", "y1", "x2", "y2"),

                   default_aes = aes(
                     width = 0.25,
                     wedgetype = "solid",
                     colour = "black",
                     fill = "black",
                     alpha = 1,
                     wedge_thickness = 0.12,
                     shorten_start = 0.15,
                     shorten_end = 0.22,
                     n_hashes = 8
                   ),

                   draw_panel = function(data, panel_params, coord) {
                     grobs <- lapply(seq_len(nrow(data)), function(i) {
                       row <- data[i, ]

                       x1 <- row$x1; y1 <- row$y1
                       x2 <- row$x2; y2 <- row$y2
                       w  <- row$width

                       dx <- x2 - x1; dy <- y2 - y1
                       len <- sqrt(dx^2 + dy^2)
                       if (len == 0) return(NULL)

                       dx <- dx / len; dy <- dy / len
                       px <- -dy; py <- dx

                       # Shorten both ends independently
                       x1s <- x1 + dx * row$shorten_start
                       y1s <- y1 + dy * row$shorten_start
                       x2s <- x2 - dx * row$shorten_end
                       y2s <- y2 - dy * row$shorten_end

                       if (row$wedgetype == "solid") {
                         df <- data.frame(
                           x = c(x1s, x2s + px * w, x2s - px * w),
                           y = c(y1s, y2s + py * w, y2s - py * w)
                         )
                         coords <- coord$transform(df, panel_params)
                         polygonGrob(coords$x, coords$y,
                                     gp = gpar(fill = row$fill, col = NA, alpha = row$alpha))
                       } else {
                         # Hashed
                         n <- row$n_hashes
                         t_vals <- seq(0.05, 0.95, length.out = n)
                         polys <- lapply(t_vals, function(t) {
                           xt <- x1s + (x2s - x1s) * t
                           yt <- y1s + (y2s - y1s) * t
                           wt <- w * t * 0.9
                           thickness <- w * (row$wedge_thickness * (0.3 + 0.7 * t))
                           wx <- px * wt; wy <- py * wt
                           tx <- dx * thickness; ty <- dy * thickness

                           df <- data.frame(
                             x = c(xt + wx, xt - wx, xt - wx - tx, xt + wx - tx),
                             y = c(yt + wy, yt - wy, yt - wy - ty, yt + wy - ty)
                           )
                           coords <- coord$transform(df, panel_params)
                           polygonGrob(coords$x, coords$y,
                                       gp = gpar(fill = row$colour, col = NA, alpha = row$alpha))
                         })
                         do.call(grobTree, polys)
                       }
                     })
                     do.call(grobTree, grobs)
                   }
  )

geom_wedge <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       ..., na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE) {

  layer(
    geom = GeomWedge,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
