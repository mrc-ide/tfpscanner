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

#' Sorts a vector of mutations
#'
#' @param   muts    String. Vector of mutations. Each string must be a separate mutation
#'   (e.g., "S:A243del"). The mutations have a prefix ("S:", "N:") and a positional description of
#'   the protein-level mutation ("T205I" for Thr to Ile mutation at position 205).
#'
#' @return  String. A vector of the same length as \code{muts}. The mutations are sorted by prefix
#'   and then by the location of the mutation.

sort_mutations <- function(muts) {
  if (length(muts) == 0) {
    return("")
  }
  pre <- sapply(strsplit(muts, split = ":"), "[", 1)
  upres <- sort(unique(pre))
  sorted_mutations <- do.call(c, lapply(upres, function(.pre) {
    .muts <- muts[pre == .pre]
    .muts1 <- sapply(strsplit(.muts,
      split = ":"
    ), "[", 2)
    sites <- regmatches(
      .muts1,
      regexpr(.muts1,
        pattern = "[0-9]+"
      )
    )
    o <- order(as.numeric(sites))
    .muts[o]
  }))

  sorted_mutations
}
