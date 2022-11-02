#' Create a sina plot against lineage
#'
#' Side-effect: this function creates a plot and stores it to a file.
#'
#' @param pldf  Data-frame. The data element of a tree plot.
#' @param output_dir  String. The directory within which the plot should be stored.
#' @param varx  String. Which variable in \code{pldf} should be plotted on the y-axis.
#' @param mut_regexp  String. Regular-expression(s) for restricting the plot to a subset of the rows
#'   in \code{pldf}. Only rows with an "allmuts" entry that matches one of these regular expressions
#'   will be presented. If NULL, all rows are presented.
#' @param lineage_regexp  String. Regular-expression(s) for restricting the lineages that are
#'   presented in the plot. Only those rows of \code{pldf} with a "lineage" entry that matches one
#'   of these regular expressions will be presented. If NULL, all rows are presented.

plot_cluster_sina <- function(pldf,
                              output_dir,
                              varx = "logistic_growth_rate",
                              mut_regexp = "S:A222V",
                              lineage_regexp = NULL) {
  sc1 <- format_cluster_sina_data(
    pldf,
    varx = varx,
    mut_regexp = mut_regexp,
    lineage_regexp = lineage_regexp
  )

  p1 <- ggplot2::ggplot(sc1) +
    ggiraph::geom_point_interactive(
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        colour = .data$mutation_lineage,
        size = .data$cluster_size + 1
      ),
      tooltip = sc1$mouseover,
      data_id = sc1$node,
      alpha = .5,
      data = sc1
    ) +
    ggplot2::xlab("Lineage and/or mutation") +
    ggplot2::ylab(varx) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::theme(
      axis.text.x =  ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  g1 <- create_widget(
      ggobj = p1,
      width_svg = 8,
      height_svg = 8
    )

  htmlwidgets::saveWidget(g1,
    file = as.character(glue::glue("{output_dir}/sina-{varx}.html"))
  )
}

#' Reformats data for a tree-plot for presentation in a sina plot
#'
#' @inheritParams plot_cluster_sina

format_cluster_sina_data <- function(pldf,
                                     varx,
                                     mut_regexp = NULL,
                                     lineage_regexp = NULL) {
  sc1 <- pldf[pldf$isTip, ]
  sc1$varx <- sc1[[varx]]
  sc1$colour_var <- ""

  if (!is.null(mut_regexp)) {
    y <- do.call(cbind, lapply(mut_regexp, function(xre) {
      ifelse(grepl(sc1$allmuts, pattern = xre), xre, "")
    }))
    ymut <- apply(y, 1, function(x) paste(x, collapse = "."))
  } else {
    ymut <- rep("", nrow(sc1))
  }

  if (!is.null(lineage_regexp)) {
    y <- do.call(cbind, lapply(lineage_regexp, function(xre) {
      ifelse(grepl(sc1$lineage, pattern = xre), xre, "")
    }))
    ylin <- apply(y, 1, function(x) paste(x, collapse = "."))
  } else {
    ylin <- rep("", nrow(sc1))
  }

  sc1$mutation_lineage <- paste(ymut, ylin, sep = "_")
  fx <- as.factor(sc1$mutation_lineage)
  sc1$x <- as.numeric(fx)
  sc1$x <- stats::rnorm(nrow(sc1), sc1$x, .15) ## seed value???
  sc1$y <- sc1[[varx]]

  sc1
}
