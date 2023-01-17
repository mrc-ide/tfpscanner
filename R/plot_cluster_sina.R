#' Create a sina plot against lineage
#'
#' @param pldf  Data-frame. The data element of a tree plot.
#' @param varx  String. Which variable in \code{pldf} should be plotted on the y-axis.
#' @param mut_regexp  String. Regular-expression(s) for restricting the plot to a subset of the rows
#'   in \code{pldf}. Only rows with an "allmuts" entry that matches one of these regular expressions
#'   will be presented. If NULL, all rows are presented.
#' @param lineage_regexp  String. Regular-expression(s) for restricting the lineages that are
#'   presented in the plot. Only those rows of \code{pldf} with a "lineage" entry that matches one
#'   of these regular expressions will be presented. If NULL, all rows are presented.
#'
#' @return   A \code{ggplot2} object storing the sina-cluster plot.

plot_cluster_sina <- function(pldf,
                              varx = "logistic_growth_rate",
                              mut_regexp = "S:A222V",
                              lineage_regexp = NULL) {
  sc1 <- format_cluster_sina_data(
    pldf,
    varx = varx,
    mut_regexp = mut_regexp,
    lineage_regexp = lineage_regexp
  )

  create_cluster_sina_ggplot(sc1, y_lab = varx)
}

#' Save a sina-cluster-plot to a file
#'
#' Saves as either an htmlwidget (in a \code{html} file) or a ggplot object (in a \code{rds} file).
#'
#' @param   ggobj   \code{ggplot2} object. Contains the plot that is to be saved.
#' @param   varx   Scalar string. Which variable is depicted in the plot?
#' @param   output_dir   File path. The directory where the plot will be stored.
#' @param   output_format   String (either \code{rds}, \code{html} or both). Default: both. In
#'   which file format(s) should the plots be saved?
#' @param   width_svg,height_svg   The width and height of the plot (only used when
#'   \code{output_format == "html"}).
#'
#' @return   Invisibly returns the file paths (one for each output format) where the plots were
#'   saved.

save_sina_plot <- function(ggobj,
                           varx,
                           output_dir,
                           output_format = c("rds", "html"),
                           width_svg = 8,
                           height_svg = 8) {
  save_widget <- function(file_path) {
    g1 <- create_widget(
      ggobj = ggobj,
      width_svg = width_svg,
      height_svg = height_svg
    )

    htmlwidgets::saveWidget(g1, file = file_path)
  }

  output_format <- match.arg(output_format, several.ok = TRUE)

  file_paths <- stats::setNames(
    file.path(output_dir, glue::glue("sina-{varx}.{output_format}")),
    output_format
  )

  for (fmt in output_format) {
    file_path <- file_paths[fmt]

    if (fmt == "rds") {
      saveRDS(ggobj, file = file_path)
    } else {
      save_widget(file_path)
    }
  }

  invisible(file_paths)
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
