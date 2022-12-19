#' Extract lists of 'defining' and 'all' mutations for each cluster in the data.frame `sc0`
#'
#' @param   sc0   data.frame. Must contains columns \code{cluster_id}, \code{defining_mutations} and
#'   \code{all_mutations}. The "mutation" columns contain a "|"-separated string of the mutations
#'   present in any given cluster. Each row represents a cluster in the phylogeny.
#'
#' @return   List of lists. Each entry in the named list corresponds to a cluster in the phylogeny.
#'   The inner lists all have entries "defining" and "all", character vectors defining the mutations
#'   that are present in the cluster.

get_mutation_list <- function(sc0) {
  required_columns <- c("defining_mutations", "all_mutations", "cluster_id")
  stopifnot(all(required_columns %in% colnames(sc0)))

  cmuts <- lapply(seq_len(nrow(sc0)), function(i) {
    list(
      defining = strsplit(sc0$defining_mutations[i], split = "\\|")[[1]],
      all = strsplit(sc0$all_mutations[i], split = "\\|")[[1]]
    )
  })
  names(cmuts) <- sc0$cluster_id

  cmuts
}
