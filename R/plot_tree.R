#' Create a non-interactive \code{ggtree} object for presenting mutation and lineage data
#'
#' @param   ggtree_data    Tree data for passing to \code{ggtree}.
#' @param   branch_col    Scalar string. The name of a column within \code{ggtree_data} for which
#'   the \code{ggtree} object should be produced.
#' @param   lins    String. Vector of lineages that are under study.
#' @param   lin_nodes   Integer. Vector of node numbers. The nodes are defined in
#'   \code{ggtree_data$node}. The order of entries matches that for \code{lins}.
#' @param   lin_node_names    String. Name of the lineage. The order of entries matches that for
#'   \code{lins}.
#' @param   shapes    Shapes for the branches and leaves in the tree.
#' @param   colours    Vector of colours.
#' @param   colour_limits   Min and max values for the colours.
#'
#' @return  A \code{ggtree} object.

create_noninteractive_ggtree <- function(ggtree_data,
                                         branch_col,
                                         lins,
                                         lin_nodes,
                                         lin_node_names,
                                         shapes,
                                         colours,
                                         colour_limits) {
  gtr1 <- ggtree::ggtree(
    ggtree_data,
    ggplot2::aes_string(colour = branch_col),
    ladderize = TRUE,
    right = TRUE,
    continuous = TRUE
  )

  gtr1.1 <- gtr1 +
    ggplot2::scale_color_gradientn(
      name = gsub(branch_col, pattern = "_", replacement = " "),
      colours = colours,
      limits = colour_limits,
      oob = scales::squish
    ) +
    ggplot2::geom_point(
      ggplot2::aes_string(
        color = branch_col,
        size = "cluster_size",
        shape = "as.factor(internal)"
      ),
      data = gtr1$data
    ) +
    ggplot2::scale_shape_manual(
      name = NULL,
      labels = NULL,
      values = shapes
    ) +
    ggplot2::scale_size(
      name = "Cluster size",
      range = c(2, 16)
    ) +
    ggplot2::ggtitle(glue::glue("{Sys.Date()}, colour: {branch_col}")) +
    ggplot2::theme(legend.position = "top")

  for (i in seq_along(lins)) {
    if (!is.na(lin_nodes[i])) {
      gtr1.1 <- gtr1.1 +
        ggtree::geom_cladelabel(
          node = lin_nodes[i],
          label = lin_node_names[i],
          offset = .00001,
          colour = "black"
        )
    }
  }

  gtr1.1
}
