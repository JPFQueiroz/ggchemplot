#' @import ggplot2
#' @import grid
#'
#' @keywords internal
GeomWedge <-
  ggplot2::ggproto("GeomWedge", ggplot2::Geom,
          required_aes = c("x1", "y1", "x2", "y2"),

          default_aes = ggplot2::aes(
            width = 0.25,
            wedgetype = "solid",
            colour = "black",
            fill = "black",
            alpha = 1,
            wedge_thickness = 0.12,
            shorten = 0.18   # shortens the WIDE end
          ),

          draw_panel = function(data, panel_params, coord) {

            grobs <- lapply(seq_len(nrow(data)), function(i) {
              row <- data[i, ]

              x1 <- row$x1; y1 <- row$y1   # narrow end (origin)
              x2 <- row$x2; y2 <- row$y2   # wide end
              w  <- row$width
              shorten <- row$shorten

              dx <- x2 - x1
              dy <- y2 - y1
              len <- sqrt(dx^2 + dy^2)
              if (len == 0) return(NULL)

              dx <- dx / len
              dy <- dy / len
              px <- -dy
              py <- dx

              # Shorten the WIDE end (near the receiving atom)
              x2s <- x2 - dx * shorten
              y2s <- y2 - dy * shorten

              if (row$wedgetype == "solid") {
                df <- data.frame(
                  x = c(x1, x2s + px * w, x2s - px * w),
                  y = c(y1, y2s + py * w, y2s - py * w)
                )

                coords <- coord$transform(df, panel_params)

                polygonGrob(coords$x, coords$y,
                            gp = gpar(fill = row$fill, col = NA,
                                      alpha = row$alpha))

              } else {
                # Hashed wedge - shortened wide end
                n <- 10
                t_vals <- seq(0, 0.95, length.out = n)   # avoid going to the very end

                polys <- lapply(t_vals, function(t) {
                  xt <- x1 + (x2s - x1) * t
                  yt <- y1 + (y2s - y1) * t

                  wt <- w * t * 0.9
                  thickness <- w * (row$wedge_thickness * (0.3 + 0.7 * t))

                  wx <- px * wt
                  wy <- py * wt
                  tx <- dx * thickness
                  ty <- dy * thickness

                  df <- data.frame(
                    x = c(xt + wx, xt - wx, xt - wx - tx, xt + wx - tx),
                    y = c(yt + wy, yt - wy, yt - wy - ty, yt + wy - ty)
                  )

                  coords <- coord$transform(df, panel_params)
                  polygonGrob(coords$x, coords$y,
                              gp = gpar(fill = row$colour, col = NA,
                                        alpha = row$alpha))
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
