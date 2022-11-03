#' Converts an interactive ggplot2 object into an HTML widget
#'
#' @param   ggobj   An interactive \code{ggplot2} object.
#' @param   width_svg,height_svg    Width and height for the output SVG.
#' @param   sizing_width    Sets the width option in \code{ggiraph::opts_sizing}.
#'
#' @return  A \code{girafe} / \code{htmlwidget} object for use in a web browser.

create_widget <- function(ggobj,
                          width_svg,
                          height_svg,
                          sizing_width = 0.8) {
  tooltip_css <- paste0(
    "background-color:black;",
    "color:grey;",
    "padding:14px;",
    "border-radius:8px;",
    "font-family:\"Courier New\",monospace;"
  )
  ggiraph::girafe(
    ggobj = ggobj,
    width_svg = width_svg,
    height_svg = height_svg,
    options = list(
      ggiraph::opts_sizing(width = sizing_width),
      ggiraph::opts_tooltip(
        css = tooltip_css,
        use_fill = FALSE
      )
    )
  )
}
