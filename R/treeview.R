#' Generate interactive tree visualisations and scatter plots for illustrating scannint statistics.
#'
#' This will produce a set of html widgets which will highlight by colour and tooltips statistics
#' such as growth rate and molecular clock outliers.
#'
#' @param e0 Path to the scanner environment produced by \code{tfpscan}. Alternatively can pass the
#'   environment directly.
#' @param branch_cols A character vector of statistics for which outputs should be produced. The
#'   logistic growth rate plot will always be produced.
#' @param mutations A character vector of mutations which will be illustrated in a heatmap
#' @param lineages A set of lineage names which will be used to subdivide outputs in scatter plots.
#' @param output_dir Outputs will be saved in this directory. Will create the directory if it does
#'   not exist.
#' @param sina_output_format String (either \code{rds}, \code{html} or both). In which formats
#'   should the sina-cluster plots be saved?
#' @param heatmap_width,heatmap_lab_offset Width and label-offset parameters for the constructed
#'   heatmap.
#'
#' @importFrom rlang .data
#'
#' @return A ggtree plot
#'
#' @export

treeview <- function(e0,
                     branch_cols = c("logistic_growth_rate", "clock_outlier"),
                     mutations = c("S:A222V", "S:Y145H", "N:Q9L", "S:E484K"),
                     lineages = c("AY\\.9", "AY\\.43", "AY\\.4\\.2"),
                     output_dir = "treeview",
                     sina_output_format = c("rds", "html"),
                     heatmap_width = .075,
                     heatmap_lab_offset = -6) {
  sina_output_format <- match.arg(sina_output_format, several.ok = TRUE)

  # require logistic growth rate, prevent non-empty
  branch_cols <- unique(c(
    "logistic_growth_rate",
    branch_cols
  ))

  dir.create(output_dir,
    showWarnings = FALSE
  )

  # load env
  if (is.character(e0)) {
    e0 <- readRDS(e0)
  }
  sc0 <- e0$Y
  cmuts <- lapply(seq_len(nrow(sc0)), function(i) {
    list(
      defining = strsplit(sc0$defining_mutations[i], split = "\\|")[[1]],
      all = strsplit(sc0$all_mutations[i], split = "\\|")[[1]]
    )
  })
  names(cmuts) <- sc0$cluster_id
  tr1 <- e0$tre

  stopifnot(all(branch_cols %in% colnames(e0$Y)))

  # pick one tip from each cluster
  sc0$representative <- NA
  for (i in seq_len(nrow(sc0))) {
    if (sc0$external_cluster[i]) {
      tu <- strsplit(sc0$tip[i], split = "\\|")[[1]]
      .sts <- e0$sts[tu]
      sc0$representative[i] <- tu[which.max(.sts)] # most recent sample in cluster
    }
  }
  tr2 <- ape::keep.tip(
    tr1,
    stats::na.omit(sc0$representative)
  )

  # subtrees
  message("Computing sub-trees")
  stres2 <- ape::subtrees(tr2, wait = TRUE)
  Nstres2 <- sapply(stres2, ape::Ntip)

  # for each cluster find the node in tr2 representing it
  ## internal nodes first
  for (i in seq_along(stres2)) {
    inode <- i + ape::Ntip(tr2)
    uv <- stats::na.omit(sc0$node_number[match(
      stres2[[i]]$tip.label,
      sc0$representative
    )])
    shared_anc <- Reduce(
      intersect,
      e0$ancestors[uv]
    )
    shared_anc2 <- setdiff(
      intersect(
        shared_anc,
        sc0$node_number
      ),
      uv
    )
    if (length(shared_anc2) > 0) {
      a <- shared_anc2[which.min(e0$ndesc[shared_anc2])]
      sc0$tr2mrca[sc0$node_number == a] <- inode
    }
  }

  ## tips (takes precedence if overlap in tr2mrca )
  i <- which(sc0$representative %in% tr2$tip.label)
  sc0$tr2mrca[i] <- match(
    sc0$representative[i],
    tr2$tip.label
  )

  # tree data frames
  ## tips
  sc2 <- sc0[!is.na(sc0$tr2mrca), ]
  sc2$date_range <- sapply(
    seq_len(nrow(sc2)),
    function(i) glue::glue("{sc2$least_recent_tip[i]} -> {sc2$most_recent_tip[i]}")
  )
  tdvars <- unique(c(
    branch_cols,
    "logistic_growth_rate",
    "clock_outlier",
    "cluster_size",
    "date_range",
    "cluster_id",
    "region_summary",
    "cocirc_lineage_summary",
    "lineage",
    "tr2mrca"
  ))

  td0 <- sc2[sc2$tr2mrca <= ape::Ntip(tr2), tdvars]
  td0$lineages <- td0$lineage
  td0$cocirc_summary <- td0$cocirc_lineage_summary
  td0$node <- td0$tr2mrca
  td0$internal <- "N"
  ## internal
  td1 <- sc2[sc2$tr2mrca > ape::Ntip(tr2), tdvars]
  if (nrow(td1) > 0) {
    td1$lineages <- td1$lineage
    td1$cocirc_summary <- td1$cocirc_lineage_summary
    td1$node <- td1$tr2mrca
    td1$internal <- "Y"
    td1$cluster_size <- 0
    x <- setdiff(
      (ape::Ntip(tr2) + 1):(ape::Ntip(tr2) + ape::Nnode(tr2)),
      td1$node
    ) # make sure every node represented
    td1 <- merge(td1,
      data.frame(node = x),
      all = TRUE
    )
    td <- rbind(td0, td1)
  } else {
    td <- td0
  }
  td <- td[order(td$node), ] # important

  # rescale clock ?
  td$clock_outlier <- scale(td$clock_outlier) / 2

  # interpolate missing values  &  repair cluster sizes
  td$logistic_growth_rate[(td$node <= ape::Ntip(tr2)) & (is.na(td$logistic_growth_rate))] <- 0
  td$clock_outlier[(td$node <= ape::Ntip(tr2)) & (is.na(td$clock_outlier))] <- 0
  for (ie in ape::postorder(tr2)) {
    a <- tr2$edge[ie, 1]
    u <- tr2$edge[ie, 2]
    td$cluster_size[a] <- td$cluster_size[a] + td$cluster_size[u]
    for (vn in branch_cols) {
      if (is.na(td[[vn]][a])) {
        td[[vn]][a] <- td[[vn]][u]
      }
    }
  }

  # cols for continuous stats
  cols <- rev(c("red", "orange", "green", "cyan", "blue"))

  # lineages for clade labels
  td$lineages1 <- sapply(strsplit(td$lineages, split = "\\|"), "[", 1)
  sc0$lineage1 <- sapply(strsplit(sc0$lineage, split = "\\|"), "[", 1)

  tablin <- table(td$lineages1[!(sc0$lineage1 %in% c("None", "B.1"))])
  tablin <- tablin[order(tablin)]

  # find a good internal node to represent the mrca of each lineage
  lins <- list(
    lineages = names(tablin),
    nodes = c(),
    node_names = c()
  )

  for (lin in lins$lineages) {
    whichrep <- stats::na.omit(sc0$representative[sc0$lineage1 == lin])
    res <- NULL
    if (length(whichrep) == 1) {
      res <- which(tr2$tip.label == whichrep)
    } else {
      res <- ape::getMRCA(tr2, whichrep)
      # ~ 			get desc tips
      itr2 <- which(sc0$tr2mrca == res)
      if (length(itr2) == 1) {
        tips <- strsplit(sc0[itr2, "tips"], split = "\\|")[[1]]
        # ~ 			check lineage
        tipslins <- e0$amd$lineage[match(tips, e0$amd$sequence_name)]
        # ~ 			check major lineage is lin
        txtl <- table(tipslins)
        if (names(txtl[which.max(txtl)]) != lin) {
          res <- NULL
        }
      } else {
        res <- NULL
      }
    }
    if (!is.null(res)) {
      lins$nodes <- c(lins[["nodes"]], res)
      lins$node_names <- c(lins[["node_names"]], lin)
    }

    res
  }

  message("Generating figures")

  ggtree_data <- dplyr::full_join(tr2, td, by = "node")
  n_leaves <- ape::Ntip(tr2)

  create_trees_curried <- function(...) {
    create_trees(
      ...,
      ggtree_data = ggtree_data,
      n_leaves = n_leaves,
      lins = lins,
      sc0 = sc0,
      cmuts = cmuts,
      mut_regex = mutations,
      colours = cols,
      heatmap_width = heatmap_width,
      heatmap_lab_offset = heatmap_lab_offset
    )
  }

  suppressWarnings({
    lgr_trees <- create_trees_curried(
      branch_col = "logistic_growth_rate",
      colour_limits = c(-.5, .5)
    )
    save_trees(
      lgr_trees,
      branch_col = "logistic_growth_rate",
      n_leaves = n_leaves,
      output_dir = output_dir
    )

    for (branch_col in setdiff(branch_cols, c("logistic_growth_rate"))) {
      tree_list <- create_trees_curried(
        branch_col = branch_col,
        colour_limits = range(td[[branch_col]])
      )
      save_trees(
        tree_list,
        branch_col = branch_col,
        n_leaves = n_leaves,
        output_dir = output_dir
      )
    }
  })

  pl <- lgr_trees[["with_interactivity_data"]]
  pldf <- pl[["data"]]

  suppressWarnings({
    for (vn in unique(c("logistic_growth_rate", branch_cols))) {
      sina_plot <- plot_cluster_sina(
        pldf,
        varx = vn,
        mut_regexp = mutations,
        lineage_regexp = lineages
      )

      save_sina_plot(
        sina_plot,
        varx = vn,
        output_dir = output_dir,
        output_format = sina_output_format
      )
    }
  })

  invisible(pl)
}

# ~ treeview( e0 = 'tfpscan-2021-11-25/scanner-env-2021-11-25.rds'
# ~ , branch_cols = c('logistic_growth_rate')
# ~ , mutations = c( 'S:Y145H', 'N:Q9L')
# ~ , lineages = c( 'AY\\.9' , 'AY\\.43' )
# ~ , output_dir = 'treeview'
# ~ )
