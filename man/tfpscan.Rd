% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfpscanner.R
\name{tfpscan}
\alias{tfpscan}
\title{Compute scanner statistics for nodes in tree.}
\usage{
tfpscan(
  tre,
  amd,
  min_descendants = 100,
  max_descendants = 20000,
  min_cluster_age_yrs = 1/12,
  min_date = NULL,
  max_date = NULL,
  min_blen = 1/30000/2,
  ncpu = 1,
  output_dir = paste0("tfpscan-", Sys.Date()),
  num_ancestor_comparison = 500,
  factor_geo_comparison = 5,
  Tg = 7/365,
  report_freq = 50,
  mutation_cluster_frequency_threshold = 0.75,
  test_cluster_odds = c(),
  test_cluster_odds_value = c(),
  root_on_tip = "Wuhan/WH04/2020",
  root_on_tip_sample_time = 2020,
  detailed_output = FALSE,
  compute_gam = TRUE,
  compute_cluster_muts = TRUE
)
}
\arguments{
\item{tre}{A phylogeny in ape::phylo or treeio::treedata form. If not rooted, must provide an
outgroup (see below)}

\item{amd}{A data frame containing required metadata for each tip in tree: sequence_name,
sample_date, region. Optional metadata includes: sample_time(numeric), lineage, mutations.}

\item{min_descendants}{Clade must have at least this many tips}

\item{max_descendants}{Clade can have at most this many tips}

\item{min_cluster_age_yrs}{Only include clades that have sample tips that span at least this
value}

\item{min_date}{Only include samples after this data}

\item{max_date}{Only include samples before and including this date}

\item{min_blen}{Only compute statistics for nodes descended from branches of at least this length}

\item{ncpu}{number cpu for multicore ops}

\item{output_dir}{Path to directory where results will be saved}

\item{num_ancestor_comparison}{When finding comparison sample for molecular clock stat, make sure
sister clade has at least this many tips}

\item{factor_geo_comparison}{When finding comparison sample based on geography, make sure sample
has this factor times the number within clade of interest}

\item{Tg}{Approximate generation time in years. Growth rate is reported in these units.}

\item{report_freq}{Print progress for every n'th node}

\item{mutation_cluster_frequency_threshold}{If mutation is detected with more than this frequency
within a cluster it may be called as a defining mutation}

\item{test_cluster_odds}{A character vector of variable names in \code{amd}. The odds of a sample
belonging to each cluster given tis variable will be estimated using conditional logistic
regression and adjusting for time.}

\item{test_cluster_odds_value}{Vector of same length as \code{test_cluster_odds}. This variable
will be dichotomised by testing for equality of the variable with this value (e.g.
vaccine_breakthrough == 'yes'). If NULL, the variable is assumed to be continuous (e.g.
patient_age).}

\item{root_on_tip}{If input tree is not rooted, will root on this tip}

\item{root_on_tip_sample_time}{Numeric time that tip was sampled}

\item{detailed_output}{If TRUE will provide detailed figures for each cluster}

\item{compute_gam}{Boolean. Should growth rates be summarised using a generalised additive model
(in addition to the logistic model)?}

\item{compute_cluster_muts}{Boolean. Should clusters of mutations (defining mutations etc) be
identified?}
}
\value{
Invisibly returns a data frame with cluster statistics.
}
\description{
Takes standard inputs in the form of a rooted phylogeny and data frame with required metadata
(see below).
Computes a logistic growth rate and a statistic for outlying values in a molecular clock (root to
tip regression).
Comparison sample is matched by time and region.
If using parallel computation (ncpu>1), the code should be called via mpirun, e.g.
mpirun -n <ncpu+1> R --slave -f <script>
}
