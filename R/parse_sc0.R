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

#' Extract data.frames containing the 'defining' and 'all' mutations for each cluster
#'
#' @inheritParams   get_mutation_list
#'
#' @return   List of two data.frames with names "all" and "defining". These contain all mutations-
#'   and just the defining mutations for each cluster in the phylogeny. The data.frames have
#'   identical structure with column names \code{cluster_id} and \code{mutation}.

get_mutation_tables <- function(sc0) {
  mutation_list <- get_mutation_list(sc0)
  defining_mutations <- purrr::map_df(
    mutation_list, ~ tibble::tibble(mutation = .x[["defining"]]),
    .id = "cluster_id"
  )
  all_mutations <- purrr::map_df(
    mutation_list, ~ tibble::tibble(mutation = .x[["all"]]),
    .id = "cluster_id"
  )

  list(
    defining = defining_mutations,
    all = all_mutations
  )
}
