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

  p1 <- create_cluster_sina_ggplot(sc1, y_lab = varx)

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

  sc1$mutation_lineage <- get_mutation_and_lineage(
    sc1,
    mut_regexp = mut_regexp,
    lineage_regexp = lineage_regexp
  )

  fx <- as.factor(sc1$mutation_lineage)
  sc1$x <- as.numeric(fx)
  sc1$x <- stats::rnorm(nrow(sc1), sc1$x, .15) ## seed value???
  sc1$y <- sc1[[varx]]

  sc1
}

#' Creates an interactive sina plot using ggplot2 and ggiraph
#'
#' @param   sc1   Data-frame. As created using \code{format_cluster_sina_data}.
#' @param   x_lab,y_lab   Scalar string. Labels for the x- and y-axis.
#'
#' @return  an interactive ggplot2 object

create_cluster_sina_ggplot <- function(sc1,
                                       x_lab = "Lineage and/or mutation",
                                       y_lab = "logistic_growth_rate") {
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
    ggplot2::xlab(x_lab) +
    ggplot2::ylab(y_lab) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
    ggplot2::theme(
      axis.text.x =  ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  p1
}

#' Returns a string vector containing mutation and lineage information
#'
#' The strings are "_"-separated, with mutations on the left and lineages on the right.
#' In a given string, there is an entry for each regexp that matches.
#' The regex matches are "."-separated (there is n-1 dots on the LHS if there are n `mut_regexps`).
#'
#' @inheritParams   plot_cluster_sina
#' @param   x   A data-frame, it must contain a "allmuts" and a "lineage" column.
#'
#' @return   Vector of strings, one for each row of \code{x}.

get_mutation_and_lineage <- function(x,
                                     mut_regexp = NULL,
                                     lineage_regexp = NULL) {
  if (!is.null(mut_regexp)) {
    y <- do.call(cbind, lapply(mut_regexp, function(xre) {
      ifelse(grepl(x[["allmuts"]], pattern = xre), xre, "")
    }))
    ymut <- apply(y, 1, function(z) paste(z, collapse = "."))
  } else {
    ymut <- rep("", nrow(x))
  }

  if (!is.null(lineage_regexp)) {
    y <- do.call(cbind, lapply(lineage_regexp, function(xre) {
      ifelse(grepl(x[["lineage"]], pattern = xre), xre, "")
    }))
    ylin <- apply(y, 1, function(z) paste(z, collapse = "."))
  } else {
    ylin <- rep("", nrow(x))
  }

  paste(ymut, ylin, sep = "_")
}
