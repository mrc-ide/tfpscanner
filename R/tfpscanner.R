#' Compute scanner statistics for nodes in tree.
#'
#' Takes standard inputs in the form of a rooted phylogeny and data frame with required metadata (see below).
#' Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to tip regression).
#' Comparison sample is matched by time and region.
#' If using parallel computation (ncpu>1), the code should be called via mpirun, e.g.
#' mpirun -n <ncpu+1> R --slave -f <script>
#'
#' @param tre A phylogeny in ape::phylo or treeio::treedata form. If not rooted, must provide an outgroup (see below)
#' @param amd A data frame containing required metadata for each tip in tree: sequence_name, sample_date, region. Optional metadata includes: sample_time(numeric), lineage, mutations.
#' @param min_descendants Clade must have at least this many tips
#' @param max_descendants Clade can have at most this many tips
#' @param min_date Only include samples after this data
#' @param max_date Only include samples before and including this date
#' @param min_cluster_age_yrs Only include clades that have sample tips that span at least this value
#' @param min_blen Only compute statistics for nodes descended from branches of at least this length
#' @param ncpu number cpu for multicore ops
#' @param output_dir Path to directory where results will be saved
#' @param num_ancestor_comparison When finding comparison sample for molecular clock stat, make sure sister clade has at least this many tips
#' @param factor_geo_comparison  When finding comparison sample based on geography, make sure sample has this factor times the number within clade of interest
#' @param Tg Approximate generation time in years. Growth rate is reported in these units.
#' @param report_freq Print progress for every n'th node
#' @param mutation_cluster_frequency_threshold If mutation is detected with more than this frequency within a cluster it may be called as a defining mutation
#' @param test_cluster_odds A character vector of variable names in \code{amd}. The odds of a sample belonging to each cluster given tis variable will be estimated using conditional logistic regression and adjusting for time.
#' @param test_cluster_odds_values Vector of same length as \code{test_cluster_odds}. This variable will be dichotomised by testing for equality of the variable with this value (e.g. vaccine_breakthrough == 'yes'). If NULL, the variable is assumed to be continuous (e.g. patient_age).
#' @param root_on_tip If input tree is not rooted, will root on this tip
#' @param root_on_tip_sample_time Numeric time that tip was sampled
#' @param detailed_output If TRUE will provide detailed figures for each cluster
#' @import ape lubridate glue mgcv
#' @import ggplot2 ggtree phangorn
#' @import foreach doMPI
#' @return Invisibly returns a data frame with cluster statistics.
#' @export
#'

tfpscan <- function(tre,
                    amd,
                    min_descendants = 100,
                    max_descendants = 20e3,
                    min_cluster_age_yrs = 1 / 12,
                    min_date = NULL,
                    max_date = NULL,
                    min_blen = 1 / 30e3 / 2,
                    ncpu = 1,
                    output_dir = paste0("tfpscan-", Sys.Date()),
                    num_ancestor_comparison = 500,
                    factor_geo_comparison = 5,
                    Tg = 7 / 365,
                    report_freq = 50,
                    mutation_cluster_frequency_threshold = 0.75,
                    test_cluster_odds = c(),
                    test_cluster_odds_value = c(),
                    root_on_tip = "Wuhan/WH04/2020",
                    root_on_tip_sample_time = 2020,
                    detailed_output = FALSE,
                    compute_gam = TRUE,
                    compute_cluster_muts = TRUE) {
  message(paste("Starting scan", Sys.time()), "\n")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  max_time <- Inf
  if (!is.null(max_date)) {
    max_time <- decimal_date(max_date)
  } else {
    max_date <- Sys.Date()
  }

  min_time <- -Inf
  if (!is.null(min_date)) {
    min_time <- decimal_date(min_date)
  }

  # tree data
  # load algn md
  amd <- amd[!is.na(amd$sequence_name), ]
  amd$sample_date <- as.Date(amd$sample_date)
  stopifnot(all(tre$tip.label %in% amd$sequence_name))
  if (!("sample_time" %in% colnames(amd))) {
    amd$sample_time <- decimal_date(amd$sample_date)
  }
  amd$sts <- amd$sample_time
  if (!("lineage" %in% colnames(amd))) {
    amd$lineage <- "lineage_not_provided"
  }
  # exclude missing dates
  amd <- amd[!is.na(amd$sample_time), ]

  # filter by sample time
  amd <- amd[(amd$sample_time >= min_time) & (amd$sts <= max_time), ]
  sts <- setNames(amd$sample_time, amd$sequence_name)

  # mutations var
  if (!("mutations" %in% colnames(amd))) {
    amd$mutations <- "-"
  }

  # retain only required variables
  # amd <- amd[ , unique( c('sequence_name', 'sample_time', 'sample_date', 'region', 'lineage', 'mutations', test_cluster_odds) ) ]
  amd <- amd[amd$sequence_name %in% tre$tip.label, ]

  # prune tree
  if (!is.rooted(tre)) {
    if (!(root_on_tip %in% amd$sequence_name)) {
      stopifnot(root_on_tip %in% tre$tip.label)
      amd <- rbind(amd, data.frame(
        sequence_name = root_on_tip,
        sample_time = root_on_tip_sample_time,
        sample_date = as.Date(lubridate::date_decimal(root_on_tip_sample_time)),
        region = NA,
        lineage = NA
      ))
    }
  }
  tre <- keep.tip(tre, intersect(tre$tip.label, amd$sequence_name))
  tr2 <- tre
  if (!is.rooted(tre)) {
    stopifnot(root_on_tip %in% tre$tip.label)
    if (!(root_on_tip %in% tre$tip.label)) {
      stop("Outgroup sequence missing from input tree.")
    }
    tr2 <- root(tre, outgroup = root_on_tip, resolve.root = TRUE)
    tre <- tr2
  }

  # var's for fast lookup time and region
  itr2 <- match(tr2$tip.label, amd$sequence_name)
  sts <- setNames(amd$sample_time[itr2], tr2$tip.label)
  sequence_name2region <- setNames(amd$region[itr2], tr2$tip.label)

  # root to tip
  ndel <- node.depth.edgelength(tre)

  message(paste("Loaded tree & filtered by inclusion criteria", Sys.time()))

  # data structures to quickly look up tree data
  # copied/adapted from treestructure
  {
    n <- Ntip(tre)
    nnode <- Nnode(tre)
    poedges <- tre$edge[ape::postorder(tre), ]
    preedges <- tre$edge[rev(ape::postorder(tre)), ]


    # PRECOMPUTE for each node
    # num descendents
    ndesc <- rep(0, n + nnode) # n tips descending
    ndesc[1:n] <- 1
    Ndesc <- rep(1, n + nnode) # include internal nodes
    # max and min sample times of descendants :
    tlsts <- sts[tre$tip.label]
    max_desc_time <- rep(-Inf, n + nnode)
    max_desc_time[1:n] <- tlsts
    min_desc_time <- rep(Inf, n + nnode)
    min_desc_time[1:n] <- tlsts
    # postorder traverse
    for (ie in 1:nrow(poedges)) {
      a <- poedges[ie, 1]
      u <- poedges[ie, 2]
      ndesc[a] <- ndesc[a] + ndesc[u]
      Ndesc[a] <- Ndesc[a] + Ndesc[u]
      max_desc_time[a] <- max(max_desc_time[a], max_desc_time[u])
      min_desc_time[a] <- min(min_desc_time[a], min_desc_time[u])
    }
    clade_age <- max_desc_time - min_desc_time # time span of descendant tips


    # descendants; note also counts mrca node
    descendants <- lapply(1:(n + nnode), function(u) integer(Ndesc[u])) # pre allocate
    for (u in 1:(n + nnode)) {
      descendants[[u]][1] <- u
    }
    Ndesc_index <- rep(2, n + nnode)
    for (ie in 1:nrow(poedges)) {
      a <- poedges[ie, 1]
      u <- poedges[ie, 2]
      ## mem efficient way to fill in values for a
      i0 <- Ndesc_index[a]
      i1 <- Ndesc[u] + i0 - 1
      descendants[[a]][i0:i1] <- descendants[[u]]
      Ndesc_index[a] <- i1 + 1
    }

    # descendant tips; index of tip
    descendantTips <- descendants
    for (a in (n + 1):(n + nnode)) {
      descendantTips[[a]] <- descendantTips[[a]][descendantTips[[a]] <= n]
    }
    ## tip labels, only including names in amd which meet inclusion criteria (excl NA )
    descendantSids <- lapply(1:(n + nnode), function(u) na.omit(tre$tip.label[descendantTips[[u]]]))

    # ancestors
    st0 <- Sys.time()
    ancestors <- lapply(1:(n + nnode), function(u) integer())
    for (ie in 1:nrow(preedges)) {
      a <- preedges[ie, 1]
      u <- preedges[ie, 2]
      ancestors[[u]] <- c(ancestors[[a]], a)
    }
    st1 <- Sys.time()

    ## vector samp time desc
    descsts <- NULL # deprecate
  }

  message(paste("Derived lookup variables", Sys.time()))

  .get_comparator_ancestor <- function(u, num_comparison = num_ancestor_comparison) {
    nu <- ndesc[u]
    asu <- ancestors[[u]]
    na <- -1
    for (a in asu) {
      na <- ndesc[a]
      nc <- na - nu
      if (nc >= num_comparison) {
        break
      }
    }
    if (na < nu) {
      if (na > 0) {
        message("Failed to find comparator ancestor with more tips. Returning NA. Node: ", u)
      }
      return(NA)
    }
    a
  }

  # matched by time and in proportion to region prevalence
  .get_comparator_sample <- function(u, nX = factor_geo_comparison) {
    nu <- ndesc[u]
    tu <- descendantSids[[u]]
    stu <- sts[tu]
    # ~ 		 stu = descsts [[ u ]]
    minstu <- min(na.omit(stu))
    maxstu <- max(na.omit(stu))

    rtu <- sequence_name2region[tu]
    ta_r <- names(sequence_name2region[sequence_name2region %in% unique(rtu)])
    ta_t <- names(sts)[sts >= minstu & sts <= maxstu]
    ta0 <- setdiff(intersect(ta_r, ta_t), tu)
    ta0r <- sequence_name2region[ta0]

    # weight for region
    w <- table(rtu)
    w <- w / (table(ta0r)[names(w)])
    w[is.na(w)] <- 0
    wa <- w[ta0r]
    wa[is.na(wa)] <- 0

    na <- min(length(ta0), nu * nX)
    if (na < nu) {
      return(NULL)
    }
    ta <- sample(ta0, replace = FALSE, size = na, prob = wa)

    ta
  }

  .freq_figure <- function(.s1, variable = "type", value = "clade") {
    estdf <- data.frame(time = .s1$time, estimated_logodds = .s1$estimated)
    estdf <- estdf[!duplicated(estdf$time), ]
    estdf <- estdf[order(estdf$time), ]
    estdf <- cbind(estdf, t(sapply(1:nrow(estdf), function(k) {
      .s2 <- .s1[.s1$time == estdf$time[k], ]
      n <- nrow(.s2)
      lo <- estdf$estimated_logodds[k]
      f <- exp(lo) / (1 + exp(lo))
      ub <- qbinom(.975, size = n, prob = f) / n
      lb <- qbinom(.025, size = n, prob = f) / n
      n1 <- sum(.s2[[variable]] == value)
      n2 <- n - n1
      ef <- n1 / n
      c(
        lb = log(lb / (1 - lb)), ub = log(ub / (1 - ub)), n = n,
        weights = 1 / sqrt(f * (1 - f) / n),
        logodds = log(ef / (1 - ef))
      )
    })))

    p <- ggplot() +
      geom_point(aes(x = as.Date(date_decimal(time)), y = logodds, size = n), data = estdf) +
      theme_minimal() +
      theme(legend.pos = "") +
      ylab("Relative cluster frequency (log odds)") +
      xlab("") +
      # ggtitle( paste0('Frequency of ', variable, '=', value) ) +
      geom_path(aes(x = as.Date(date_decimal(time)), y = estimated_logodds), data = estdf, color = "blue", size = 1) +
      geom_ribbon(aes(x = as.Date(date_decimal(time)), ymin = lb, ymax = ub), data = estdf, alpha = .25)

    p
  }

  .logistic_growth_stat <- function(u, ta = NULL, generation_time_scale = Tg) {
    if (is.null(ta)) {
      ta <- .get_comparator_sample(u)
    }
    if (is.null(ta)) { # if still null, cant get good comparison
      return(
        list(
          lgr = NA, lgrp = .5, gam_r = NA,
          AIC = NA,
          AIC_gam = NA,
          growthrates = setNames(c(NA, NA), c("Logistic", "GAM")),
          relative_model_support = NA,
          plot = NULL
        )
      )
    }
    tu <- descendantSids[[u]]
    sta <- sts[ta]
    stu <- sts[tu]
    X <- data.frame(time = c(sta, stu), type = c(rep("control", length(ta)), rep("clade", length(tu))))
    X <- na.omit(X)
    m <- glm(type == "clade" ~ time, data = X, family = binomial(link = "logit"))
    s <- summary(m)
    rv <- unname(coef(m)[2] * generation_time_scale)
    p <- NA
    if (is.na(rv)) {
      message("NA growth stat, node: ", u)
    } else {
      p <- s$coefficients[2, 4]
    }
    ## time dep growth ; needs a larger sample size
    X$estimated <- predict(m)
    if (compute_gam & (length(tu) > 50)) {
      m1 <- mgcv::gam(type == "clade" ~ s(time, bs = "bs", k = 4, m = 1), family = binomial(link = "logit"), data = X)
      X$estimated <- predict(m1)

      tout <- seq(min(X$time), max(X$time), length = 5)
      tout1 <- tout[4] + diff(tout)[1] / 2
      dlo <- diff(predict(m1, newdata = data.frame(type = NA, time = c(tout1, max(tout)))))
      r <- dlo * Tg / ((max(tout) - tout1))
      aic <- c(AIC(m), AIC(m1))
    } else {
      r <- NA
      aic <- c(AIC(m), Inf)
    }

    list(
      lgr = rv, lgrp = p, gam_r = r,
      AIC = aic[1],
      AIC_gam = aic[2],
      growthrates = setNames(c(rv, r), c("Logistic", "GAM")),
      relative_model_support = setNames(exp(min(aic) - aic) / 2, c("Logistic", "GAM")),
      plot = .freq_figure(X)
    )
  }

  # log median p-value of rtt predicted divergence of tips under u
  .clock_outlier_stat <- function(u, a = NULL) {
    if (is.null(a)) {
      a <- .get_comparator_ancestor(u)
    }
    if (is.na(a)) {
      return(NA)
    }

    tu <- descendantSids[[u]]
    ta <- setdiff(descendantSids[[a]], tu)

    sta <- sts[ta]
    stu <- sts[tu]
    # ~ 		stu = descsts [[ u ]]

    iu <- match(tu, tre$tip.label)
    ia <- match(ta, tre$tip.label)
    ndelu <- ndel[iu]
    ndela <- ndel[ia]
    m <- lm(ndela ~ sta)
    r2 <- summary(m)$r.squared
    oosp <- predict(m, newdata = data.frame(sta = unname(stu))) - ndelu
    sqrt(mean(oosp^2)) * sign(mean(oosp))
  }

  # Compute 'proportionality' statistics for node (u) and given variable (var)
  # e.g. if sample is from vaccine breakthrough, is there higher odds that sample is in clade?
  .var_proportionality_stat <- function(u, ta,
                                        var = "breakthrough2", value = TRUE,
                                        f = clade ~ time + var,
                                        form_index = 3) {
    if (is.null(ta)) {
      return(setNames(c(NA, 1), paste(sep = "_", var, c("logodds", "p"))))
    }
    tu <- descendantSids[[u]]
    sta <- sts[ta]
    stu <- sts[tu]
    X <- data.frame(
      tip = c(ta, tu),
      time = c(sta, stu),
      type = c(rep("control", length(ta)), rep("clade", length(tu)))
    )
    X$var <- (amd[[var]][match(X$tip, amd$sequence_name)])
    if (!is.null(value)) {
      X$var <- (X$var == value)
    }
    X$clade <- (X$type == "clade")

    m <- glm(f, data = X, family = binomial(link = "logit"))
    s <- summary(m)
    rv <- unname(coef(m)[form_index])
    p <- NA
    if (is.na(rv)) {
      message("NA proportionality stat, node: ", u)
    } else {
      p <- s$coefficients[form_index, 4]
    }
    setNames(c(rv, p), paste(sep = "_", var, c("logodds", "p")))
  }

  .cluster_muts <- function(u, a = NULL, mut_variable = "mutations") {
    tu <- descendantSids[[u]]
    mdf.u <- amd[amd$sequence_name %in% tu, ]
    ## find comparator ancestor
    if (is.null(a)) {
      a <- .get_comparator_ancestor(u)
    }
    asids <- setdiff(descendantSids[[a]], descendantSids[[u]])
    mdf.a <- amd[amd$sequence_name %in% asids, ]
    if (nrow(mdf.a) == 0 | nrow(mdf.u) == 0) {
      return(list(defining = NA, all = NA))
    }
    vtabu <- sort(table(do.call(c, strsplit(mdf.u[[mut_variable]], split = "\\|"))) / nrow(mdf.u))
    vtaba <- sort(table(do.call(c, strsplit(mdf.a[[mut_variable]], split = "\\|"))) / nrow(mdf.a))

    umuts <- names(vtabu[vtabu > mutation_cluster_frequency_threshold])
    defining_muts <- setdiff(
      names(vtabu[vtabu > mutation_cluster_frequency_threshold]),
      names(vtaba[vtaba > mutation_cluster_frequency_threshold])
    )
    list(defining = defining_muts, all = umuts)
  }

  .lineage_summary <- function(tips, maxrows = 4) {
    if (is.null(tips)) {
      return("")
    }
    lins <- amd$lineage[match(tips, amd$sequence_name)]
    tx <- sort(table(lins), decreasing = TRUE) / length(lins)
    if (length(tx) > 1) {
      y <- as.data.frame(tx)
      colnames(y) <- c("Lineage", "Frequency")
    } else {
      y <- data.frame(Lineage = lins[1], Frequency = 1)
    }
    y <- y[1:min(nrow(y), maxrows), ]
    y$Frequency <- paste0(round(y$Frequency * 100), "%")
    paste(knitr::kable(y, "simple"), collapse = "\n") # convert to string
  }

  .region_summary <- function(tips, maxrows = 5) {
    if (is.null(tips)) {
      return("")
    }
    regs <- amd$region[match(tips, amd$sequence_name)]
    tx <- sort(table(regs), decreasing = TRUE) / length(regs)
    if (length(tx) > 1) {
      y <- as.data.frame(tx)
      colnames(y) <- c("Region", "Frequency")
    } else {
      y <- data.frame(Region = regs[1], Frequency = 1)
    }
    y <- y[1:min(nrow(y), maxrows), ]
    y$Frequency <- paste0(round(y$Frequency * 100), "%")
    paste(knitr::kable(y, "simple"), collapse = "\n") # convert to string
  }

  .cluster_tree <- function(tips) # tr2, amd
  {
    tr <- keep.tip(tre, tips)
    gtr <- ggtree(tr)
    mutlist <- strsplit(amd[match(tips, amd$sequence_name), ]$mutations, split = "\\|")
    # ~ sharedmut = Reduce( intersect,  mutlist )
    mutthresh <- .75
    tx <- table(do.call(c, mutlist))
    tx <- tx / Ntip(tr)
    sharedmut <- names(tx)[tx >= mutthresh] # TODO need a better way to get sharedmut
    segregatingmut <- lapply(mutlist, function(x) setdiff(x, sharedmut))
    allsegregating <- Reduce(union, segregatingmut)

    # remove stops from allseg
    allsegregating <- allsegregating[!grepl(allsegregating, patt = "[*]$")]
    allsegregating <- allsegregating[!grepl(allsegregating, patt = ":[*]")]
    annots <- rep("", Ntip(tr) + Nnode(tr))
    if (length(allsegregating) > 0) {
      allseg1 <- substr(regmatches(allsegregating, regexpr(allsegregating, patt = ":[A-Z]")), 2, 2)
      allseg2 <- regmatches(allsegregating, regexpr(allsegregating, patt = "[A-Z*]$"))
      sites_post <- regmatches(allsegregating, regexpr(allsegregating, patt = ":.*$"))
      sites_post <- substr(sites_post, 3, nchar(sites_post) - 1)
      sites_pre <- regmatches(allsegregating, regexpr(allsegregating, patt = "^.*:"))
      sites <- paste0(sites_pre, sites_post)

      aas <- c()
      aaslist <- lapply(seq_along(allsegregating), function(i) {
        do.call(rbind, lapply(
          mutlist,
          function(x) ifelse(allsegregating[i] %in% x, allseg1[i], allseg2[i])
        ))
      })
      aas <- do.call(cbind, aaslist)
      colnames(aas) <- allsegregating # sites
      rownames(aas) <- tr$tip.label
      aas <- aas[, order(sites)]
      aas[is.na(aas)] <- "X"
      sites <- sort(sites)
      # aas1 = as.AAbin( aas )
      aas2 <- as.phyDat(aas, type = "AA")
      ap <- ancestral.pars(tr, aas2, return = "phyDat")
      ap1 <- as.character(ap)
      for (ie in postorder(tr)) {
        a <- tr$edge[ie, 1]
        u <- tr$edge[ie, 2]
        j <- which(ap1[a, ] != ap1[u, ])
        # keep only S and N annots
        j <- j[grepl(sites[j], patt = "^[SN]:")]
        annots[u] <- paste(paste0(sites[j], ap1[u, j]), collapse = ",")
        if (nchar(annots[u]) == 0) annots[u] <- NA
      }
    }
    nodedf <- data.frame(node = 1:(Ntip(tr) + Nnode(tr)), annot = annots, stringsAsFactors = FALSE)
    gtr1 <- gtr %<+% nodedf
    gtr1 <- gtr1 + geom_label(aes(x = branch, label = annot, size = 6))

    gtr1 <- gtr1 + geom_tiplab(align = FALSE)
    gtr2 <- gtr1

    if (length(allsegregating) < 100) {
      gtr2 <- gheatmap(gtr1, as.data.frame(aas), width = .66, offset = 0.0005, colnames = FALSE, colnames_angle = -90, colnames_position = "top", colnames_offset_y = -2) + theme(legend.position = "none")
    }
    gtr2
  }


  # all analyses for a particular node
  .process.node <- function(u) {
    tu <- descendantSids[[u]]
    ta <- .get_comparator_sample(u)
    ulins <- amd$lineage[match(tu, amd$sequence_name)]
    alins <- amd$lineage[match(ta, amd$sequence_name)]
    lgs <- .logistic_growth_stat(u, ta)
    best_gr <- NA
    if (!is.na(lgs$lgr)) {
      best_gr <- lgs$growthrates[which.max(lgs$relative_model_support)]
    }


    reg_summary <- tryCatch(.region_summary(tu), error = function(e) as.character(e))
    cocirc_summary <- tryCatch(.lineage_summary(ta), error = function(e) as.character(e))
    lineage_summary <- tryCatch(.lineage_summary(tu), error = function(e) as.character(e))

    a <- .get_comparator_ancestor(u)
    if (compute_cluster_muts) {
      cmut <- .cluster_muts(u, a)
    } else {
      cmut <- list(defining = NA, all = NA)
    }
    X <- data.frame(
      cluster_id = as.character(u),
      node_number = u,
      parent_number = ifelse(is.null(ancestors[[u]]), NA, tail(ancestors[[u]], 1)),
      most_recent_tip = as.Date(date_decimal(max(na.omit(sts[tu])))),
      least_recent_tip = as.Date(date_decimal(min(na.omit(sts[tu])))),
      cluster_size = length(tu),
      logistic_growth_rate = best_gr,
      logistic_growth_rate_p = lgs$lgrp,
      simple_logistic_growth_rate = lgs$lgr,
      gam_logistic_growth_rate = lgs$gam_r,
      simple_logistic_model_support = lgs$relative_model_support["Logistic"],
      clock_outlier = .clock_outlier_stat(u, a),
      lineage = paste(names(sort(table(ulins), decreasing = TRUE)), collapse = "|"),
      lineage_summary = lineage_summary,
      cocirc_lineage_summary = cocirc_summary,
      region_summary = reg_summary,
      external_cluster = !(u %in% node_ancestors),
      tips = paste(tu, collapse = "|"),
      defining_mutations = paste(cmut$defining, collapse = "|"),
      all_mutations = paste(cmut$all, collapse = "|"),
      stringsAsFactors = FALSE
    )
    for (i in seq_along(test_cluster_odds)) {
      vn <- test_cluster_odds[i]
      val <- test_cluster_odds_value[i]
      clodds <- .var_proportionality_stat(u, ta, var = vn, value = val)
      X <- cbind(X, t(clodds))
    }
    rownames(X) <- as.character(u)
    if (u %in% report_nodes) { # print progress
      i <- which(nodes == u)
      message(paste("Progress", round(100 * i / length(nodes)), "%"))
    }
    if (detailed_output) {
      cldir <- glue("{output_dir}/{as.character(u)}")
      dir.create(cldir, showWarnings = FALSE)
      # summary stat data
      write.csv(data.frame(statistic = t(X[1, c("logistic_growth_rate", "simple_logistic_growth_rate", "logistic_growth_rate_p", "gam_logistic_growth_rate", "simple_logistic_model_support", "clock_outlier")])), file = glue("{cldir}/summary.csv"))
      # freq plot
      if (!is.null(lgs$plot)) {
        suppressMessages(ggsave(lgs$plot, file = glue("{cldir}/frequency.pdf")))
        suppressMessages(ggsave(lgs$plot, file = glue("{cldir}/frequency.png"), bg = "white"))
      }
      # tree plot
      if (length(tu) < 1e3) {
        gtr <- .cluster_tree(tu)
        suppressMessages(
          ggsave(gtr,
            file = glue("{cldir}/clustertree.pdf"),
            height = max(6, floor(length(tu) / 5)),
            width = min(64, max(36, sqrt(length(tu)))),
            limitsize = FALSE
          )
        )
      }
      # clock figure TODO
      # tip table
      write.csv(amd[amd$sequence_name %in% tu, ], file = glue("{cldir}/sequences.csv"))
      # reg summary
      write.csv(reg_summary, file = glue("{cldir}/regional_composition.csv"))
      # lineage summary
      write.csv(lineage_summary, file = glue("{cldir}/lineage_composition.csv"))
      # cocirc lineage summary
      write.csv(cocirc_summary, file = glue("{cldir}/cocirculating_lineages.csv"))
    }
    X
  }

  # main analysis thread
  ## compute stats for subset of nodes based on size and age
  nodes <- which((ndesc >= min_descendants) & (ndesc <= max_descendants) & (clade_age >= min_cluster_age_yrs))
  nodes_blenConstraintSatisfied <- tre$edge[tre$edge.length > min_blen, 2]
  nodes <- intersect(nodes, nodes_blenConstraintSatisfied)
  report_nodes <- nodes[seq(1, length(nodes), by = report_freq)] # progress reporting
  node_ancestors <- do.call(c, lapply(nodes, function(u) ancestors[[u]]))
  if (ncpu > 1) {
    message("Initiating MPI cluster")
    mpiclust <- startMPIcluster(count = ncpu)
    registerDoMPI(mpiclust)
    foreach(
      u = nodes,
      .combine = rbind,
      .packages = c("lubridate", "glue", "mgcv", "ggplot2", "ggtree", "phangorn"),
      .export = c("output_dir"),
      .errorhandling = "remove",
      .verbose = TRUE
    ) %dopar% {
      tryCatch(.process.node(u),
        error = function(e) {
          saveRDS(e, file = glue("{u}-err.rds"))
          return(e)
        }
      )
    } -> Y
  } else {
    Y <- c()
    for (u in nodes) {
      X <- tryCatch(.process.node(u), error = function(e) {
        saveRDS(e, file = glue("{u}-err.rds"))
        return(e)
      })
      Y <- rbind(Y, X)
    }
  }
  if (!all(nodes %in% Y$node_number)) {
    stop("Statistics not computed for all nodes. Possible that memory exceeded with ncpu > 1. Try with ncpu = 1.")
  }
  Y <- Y[order(Y$logistic_growth_rate, decreasing = TRUE), ]

  ofn1 <- glue("{output_dir}/scanner-{max_date}.rds")
  ofn3 <- glue("{output_dir}/scanner-env-{max_date}.rds")
  saveRDS(Y, file = ofn1)
  message("saving image ... ")
  # save internal variables and functions
  e0 <- environment()
  saveRDS(e0, file = ofn3)
  message(glue("Data written to {ofn1} and {ofn3}. Returning data frame invisibly."))

  if (ncpu > 1) {
    closeCluster(mpiclust)
    mpi.finalize()
  }
  invisible(Y)
}




#' Run a fast treedater/mlesky analysis for a given node in the scanner output
#' @import mlesky
#' @import treedater
#' @export
get_clusternode_mlesky <- function(u = 406318, scanner_env = readRDS("scanner-env-2021-03-03.rds")) {
  e1 <- as.environment(scanner_env)
  attach(e1)

  mr <- 5.9158E-4
  utre <- keep.tip(tre, descendantSids[[u]])
  sample_times <- sts[utre$tip.label]

  tr <- di2multi(utre, tol = 1e-05)
  tr <- unroot(multi2di(tr))
  tr$edge.length <- pmax(1 / 29000 / 5, tr$edge.length)
  tr3 <- dater(unroot(tr), sts[tr$tip.label],
    s = 29000, omega0 = mr,
    numStartConditions = 0, meanRateLimits = c(mr, mr + 1e-6), ncpu = 6
  )

  msg <- mlskygrid(tr3,
    tau = NULL, tau_lower = .001, tau_upper = 10, sampleTimes = sts[tr3$tip.label],
    res = 10, ncpu = 3
  )

  list(mlesky = msg, timetree = tr3, tree = tr)
}
