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

#' Adds data to a ggtree object to allow mouse-over tooltips etc when presented interactively
#'
#' @param   ggobj   A \code{ggtree} object.
#' @param   branch_col    Scalar string. The name of a column within \code{ggobj$data} defining the
#'   statistic under study here (`logistic_growth_rate`, `clock_outlier`).
#' @param   sc0,cmuts   Data-frames.
#' @param   mut_regex   Regular expression. Defines the mutations under study here.
#'
#' @return  A ggtree object. The \code{$data} entry has additional entries (\code{mouseover},
#'   \code{colour_var}, \code{defmuts}, \code{allmuts}) that are used when presented interactively
#'   by \code{ggiraph}.

append_interactivity_data <- function(ggobj,
                                      branch_col,
                                      sc0,
                                      cmuts,
                                      mut_regex) {
  # make mouseover info
  ## standard meta data
  ttdfs <- apply(ggobj$data, 1, FUN = function(x) {
    z <- as.list(x)
    lgr <- as.numeric(z$logistic_growth_rate)
    # TODO: replace with() with explicit z$cluster_id etc
    y <- with(
      z,
      data.frame(
        `Cluster ID` = glue::glue("#{cluster_id}"),
        `Cluster size` = cluster_size,
        `Date range` = date_range,
        `Example sequence` = label,
        `Logistic growth` = paste0(
          ifelse(lgr > 0, "+", ""),
          round(lgr * 100), "%"
        ),
        `Mol clock outlier` = clock_outlier,
        `Lineages` = lineages
      )
    )
    y <- t(y)
    colnames(y) <- ""
    tryCatch(
      paste(knitr::kable(y, "simple"), collapse = "\n"),
      error = function(e) paste(knitr::kable(y, "markdown"), collapse = "\n")
    )
  })

  ## table with geo composition
  ttregtabs <- ggobj$data$region_summary #
  ## cocirc
  ttcocirc <- ggobj$data$cocirc_summary #

  ## defining muts
  ttdefmuts <- sapply(match(ggobj$data$cluster_id, sc0$cluster_id), function(isc0) {
    if (is.na(isc0)) {
      return("")
    }
    paste(
      sep = "\n",
      "Cluster branch mutations:",
      gsub(
        x = tryCatch(
          stringr::str_wrap(
            paste(
              collapse = " ",
              sort_mutations(cmuts[[as.character(sc0$node_number[isc0])]]$defining)
            ),
            width = 60
          ),
          error = function(e) browser()
        ),
        pattern = " ",
        replacement = ", "
      ),
      "\n"
    )
  }) # end of sapply

  ttallmuts <- sapply(match(ggobj$data$cluster_id, sc0$cluster_id), function(isc0) {
    if (is.na(isc0)) {
      return("")
    }
    paste(
      sep = "\n",
      "All mutations:",
      gsub(
        x = stringr::str_wrap(
          paste(
            collapse = " ",
            sort_mutations(cmuts[[as.character(sc0$node_number[isc0])]]$all)
          ),
          width = 60
        ),
        pattern = " ",
        replacement = ", "
      ),
      "\n"
    )
  }) # end of sapply

  ggobj$data$defmuts <- ttdefmuts
  ggobj$data$allmuts <- ttallmuts
  if (!is.null(mut_regex)) {
    for (mre in mut_regex) {
      i <- which(grepl(ggobj$data$allmuts, pattern = mre))
      ggobj$data[[mre]] <- grepl(ggobj$data$allmuts, pattern = mre)
    }
  }

  # make html widget
  ggobj$data$mouseover <- sapply(seq_along(ttdfs), function(i) {
    paste0(
      "Statistics:\n", ttdfs[i],
      "\n\nGeography:\n", ttregtabs[i],
      "\n\nCo-circulating with:\n", ttcocirc[i],
      "\n\n", ttdefmuts[i],
      "\n", ttallmuts[i],
      "\n",
      collapse = "\n"
    )
  })
  ggobj$data$colour_var <- ggobj$data[[branch_col]]

  ggobj
}

#' Extract the data about viral genotypes from a \code{ggtree} object
#'
#' @param   ggobj   A \code{ggtree} object, as generated by \code{append_interactivity_data}. The
#'   \code{data} entry for this object should contain the columns "node", "label" and a column for
#'   each of the \code{mut_regex} values.
#' @param   n_leaves   Scalar integer. The number of leaves in the tree.
#' @param   mut_regex   String. Regular expression defining the mutations under study here. This
#'   should be a subset of the column-names in \code{ggobj$data}.
#'
#' @return   Data-frame.

extract_genotype_data <- function(ggobj,
                                  n_leaves,
                                  mut_regex) {
  genotype <- as.data.frame(
    ggobj$data[
      ggobj$data$node <= n_leaves,
      c("label", mut_regex)
    ]
  )
  rownames(genotype) <- genotype$label
  genotype <- genotype[, -1, drop = FALSE]

  genotype
}

#' Adds a heatmap to the right of a ggtree object
#'
#' @param   ggobj   A ggtree object.
#' @param   genotype    The heatmap data.
#' @param   heatmap_width,heatmap_lab_offset    Parameters for positioning of the heatmap.
#'
#' @return  A \code{ggtree} / \code{gg} / \code{ggplot} object with an appended heatmap.

append_heatmap <- function(ggobj,
                           genotype,
                           heatmap_width = 1,
                           heatmap_lab_offset = 0) {
  ggtree::gheatmap(
    p = ggobj,
    data = genotype,
    width = heatmap_width,
    offset = 0.0005,
    colnames_angle = -90,
    colnames_position = "top",
    colnames_offset_y = heatmap_lab_offset,
    legend_title = "Genotype"
  )
}

#' Converts a \code{ggtree} object into a \code{ggiraph} object with interactive potential
#'
#' @param   ggobj   A \code{ggtree} object.
#' @param   branch_col    Scalar string. Name of the column in \code{ggobj$data} that we are
#'   creating the plot for.
#' @param   cluster_size_range    Numeric (length-2). min and max values for cluster sizes on the
#'   chart.
#' @inheritParams   create_noninteractive_ggtree
#'
#' @return  A  \code{ggtree} object with interactive data for presentation by
#'   \code{ggiraph::girafe}.

create_interactive_ggtree <- function(ggobj,
                                      branch_col,
                                      cluster_size_range,
                                      shapes,
                                      colours,
                                      colour_limits) {
  ggobj +
    ggiraph::geom_point_interactive(
      ggplot2::aes(
        x = .data$x,
        y = .data$y,
        color = .data$colour_var,
        tooltip = .data$mouseover,
        data_id = .data$node,
        size = .data$cluster_size + 1,
        shape = as.factor(.data$internal)
      )
    ) +
    ggplot2::scale_shape_manual(
      name = NULL,
      labels = NULL,
      values = shapes
    ) +
    ggplot2::scale_size(
      name = "Cluster size",
      range = cluster_size_range
    ) +
    ggplot2::scale_color_gradientn(
      name = stringr::str_to_title(
        gsub(branch_col, pattern = "_", replacement = " ")
      ),
      colours = colours,
      limits = colour_limits,
      oob = scales::squish
    ) +
    ggplot2::theme(legend.position = "top")
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
