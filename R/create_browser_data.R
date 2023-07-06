#' Generate treeview illustrations and all metadata for presentation by tfpbrowser
#'
#' @param   e0   Path to the scanner environment produced by \code{tfpscan}. Alternatively can pass
#'   the environment directly.
#' @param   output_dir   Where should the output files be saved to? This function will create
#'   subdirectories `<output_dir>/treeview`, `<output_dir>/mutations`, `<output_dir>/sequences` and
#'   assumes that `<output_dir>/scanner_output` already exists.
#' @param   ...   Further arguments for `treeview()`.
#'
#' @export

create_browser_data <- function(e0,
                                output_dir,
                                ...) {
  dirs <- list(
    treeview = file.path(output_dir, "treeview"),
    mutations = file.path(output_dir, "mutations"),
    scanner = file.path(output_dir, "scanner_output"),
    sequences = file.path(output_dir, "sequences")
  )

  treeview_results <- treeview(e0, output_dir = dirs[["treeview"]], ...)

  # Create blank treeview images
  empty_treeview(
    treeview = "tree-logistic_growth_rate.rds",
    treeview_dir = dirs[["treeview"]]
  )

  # Save .csv containing node-lookups
  create_all_node_lookups(output_dir)

  # Save .csv containing cluster IDs and associated sequences
  create_sequences_lookup(output_dir)

  # Save .csvs containing defining- and all-mutations associated with cluster IDs
  create_mutation_files(e0, mutations_dir = dirs[["mutations"]])
}

#' Function to create a treeview with all grey nodes to display
#'
#' @param   treeview   RDS file containing an existing treeview plot
#' @param   treeview_dir   The directory where data should be read from / written to.
#' @param   types   Character vector of new variables to colour by
#'
#' @export

empty_treeview <- function(treeview,
                           treeview_dir,
                           types = c("mutations", "sequences")) {
  filename <- file.path(treeview_dir, treeview)
  stopifnot(file.exists(filename))
  g <- readRDS(filename)

  make_treeview_type <- function(type) {
    new_g <- g +
      ggplot2::scale_colour_gradient(low = "grey", high = "grey") +
      ggplot2::guides(
        colour = "none",
        fill = "none",
        shape = "none"
      ) +
      ggplot2::labs(title = glue::glue("Colour: {type}"))
    new_filename <- file.path(treeview_dir, glue::glue("tree-{type}.rds"))
    saveRDS(new_g, file = new_filename)
  }

  purrr::walk(.x = types, .f = ~ make_treeview_type(.x))
}

#' Function to create lookups for nodes for all .rds files in the treeview directory
#'
#' @param   data_dir   The directory where the data should be read from / written to. This must
#'   contain a `treeview` subdirectory.
#'
#' @export

create_all_node_lookups <- function(data_dir) {
  # get list of all widgets
  all_widgets <- available_treeview(data_dir)
  purrr::walk(.x = all_widgets, .f = ~ create_node_lookup(.x, data_dir = data_dir))
}

#' Function to create lookup for a single treeview
#'
#' @param   widgetChoice   The .rds filename for selected treeview output from radio button.
#' @param   data_dir   The directory where the data should be read from / written to.

create_node_lookup <- function(widgetChoice, data_dir) {
  dirs <- list(
    data = data_dir,
    treeview = file.path(data_dir, "treeview"),
    node_lookup = file.path(data_dir, "treeview", "node_lookup")
  )

  stopifnot(dir.exists(dirs[["treeview"]]))
  create_dir_if_missing(dirs[["node_lookup"]])

  output_basename <- stringr::str_replace(widgetChoice, ".rds", ".csv")
  files <- list(
    input = file.path(dirs[["treeview"]], widgetChoice),
    output = file.path(dirs[["node_lookup"]], output_basename)
  )

  g <- readRDS(files[["input"]])

  built <- suppressWarnings(ggplot2::ggplot_build(g))
  if (widgetChoice %in% c(
    "sina-logistic_growth_rate.rds",
    "sina-simple_logistic_growth_rate.rds",
    "sina-clock_outlier.rds"
  )) {
    ids <- built$data[1][[1]]["data_id"]
    tooltips <- built$data[1][[1]]$tooltip
    tooltip_ids <- get_cluster_ID(tooltips)
  } else {
    n_layers <- length(built$data)
    ids <- built$data[n_layers][[1]]["data_id"]
    tooltips <- built$data[n_layers][[1]]$tooltip
    tooltip_ids <- get_cluster_ID(tooltips)
  }
  ids$cluster_ids <- tooltip_ids

  readr::write_csv(ids, file = files[["output"]])
}

#' Function to save a CSV file of all sequences for all clusterIDs
#'
#' The files `<data_dir>/scanner_output/*/sequences.csv` will be combined together to create the
#' output file `<data_dir>/sequences/all_sequences.csv`.
#'
#' @param   data_dir   The data directory for the application. Must have a `scanner_output`
#'   subdirectory. Within `<data_dir>/scanner_output/` every subdirectory must contain a
#'   `sequences.csv` file.
#'
#' @export

create_sequences_lookup <- function(data_dir) {
  dirs <- list(
    input = file.path(data_dir, "scanner_output"),
    output = file.path(data_dir, "sequences")
  )
  cluster_ids <- list.dirs(
    dirs[["input"]],
    recursive = FALSE,
    full.names = FALSE
  )
  output_filepath <- file.path(dirs[["output"]], "all_sequences.csv")

  create_dir_if_missing(dirs[["output"]])

  lookup_table <- purrr::map_dfr(.x = cluster_ids, .f = ~ process_seq_table(.x, data_dir))
  readr::write_csv(lookup_table, file = output_filepath)
}

#' Function to get lookup table of clusterID and sequence
#'
#' @param   selected_folder   Folder name relating to a single clusterID
#' @param   data_dir   The data directory for the application. Must have a `scanner_output`
#'   subdirectory.

process_seq_table <- function(selected_folder, data_dir) {
  sequences <- file.path(data_dir, "scanner_output", selected_folder, "sequences.csv")
  sequences <- suppressMessages(readr::read_csv(sequences))
  if (nrow(sequences) > 0) {
    seq_names <- unique(sequences$sequence_name)
    output <- tibble::tibble(
      cluster_id = rep(selected_folder, length(seq_names)),
      sequence = seq_names
    )
    return(output)
  }
}

#' Create .csvs containing defining- and all-mutations
#'
#' @param   e0   Path to the scanner environment produced by \code{tfpscan}. Alternatively can pass
#'   the environment directly.
#' @param   mutations_dir   The directory where the mutations .csv files will be placed. The
#'   filenames are `<mutations_dir>/all_mutations.csv` and `<mutations_dir>/defining_mutations.csv`.

create_mutation_files <- function(e0,
                                  mutations_dir) {
  create_dir_if_missing(mutations_dir)

  files <- list(
    defining = file.path(mutations_dir, "defining_mutations.csv"),
    all = file.path(mutations_dir, "all_mutations.csv")
  )

  # load env
  if (is.character(e0)) {
    e0 <- readRDS(e0)
  }
  sc0 <- e0$Y

  cmut_tables <- get_mutation_tables(sc0)
  readr::write_csv(cmut_tables[["defining"]], files[["defining"]])
  readr::write_csv(cmut_tables[["all"]], files[["all"]])
}

#' function to return treeview options
#'
#' @param   data_dir   The directory containing the data for the application.

available_treeview <- function(data_dir) {
  all_trees <- list.files(
    file.path(data_dir, "treeview"),
    pattern = "\\.rds$"
  )
  all_trees <- factor(
    all_trees,
    c(
      stringr::str_subset(all_trees, "tree"),
      stringr::str_subset(all_trees, "sina")
    )
  )
  all_trees <- as.character(sort(all_trees))
  names(all_trees) <- all_trees %>%
    stringr::str_replace_all("_|-|\\.rds", " ") %>%
    stringr::str_trim() %>%
    stringr::str_to_title()
  return(all_trees)
}

#' Create a directory if it doesn't yet exist
#'
#' @param   path   The path for the required directory

create_dir_if_missing <- function(path) {
  stopifnot(length(path) == 1)

  if (!dir.exists(path)) {
    dir.create(path)
  }
}

#' Function to get node id from data_id column of ggplot
#' TO BE REMOVED AFTER TOOLTIPS TFPSCANNER PR IS MERGED
#'
#' @param   tooltip_input   Character vector of tooltip content
#'
#' @export

get_cluster_ID <- function(tooltip_input) {
  # start searching the string after the "Cluster.ID" text
  # until the next new line
  match_matrix <- stringr::str_match(tooltip_input, pattern = r"(Cluster.ID\s+#(\d+))")
  cluster_ids <- as.numeric(match_matrix[, 2])
  return(cluster_ids)
}
