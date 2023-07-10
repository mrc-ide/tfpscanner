get_default_filepaths <- function(dir_path, base_name) {
  filepaths <- c(
    noninteractive = file.path(dir_path, glue::glue("{base_name}.svg")),
    interactive_html = file.path(dir_path, glue::glue("{base_name}.html")),
    interactive_rds = file.path(dir_path, glue::glue("{base_name}.rds"))
  )

  filepaths
}

describe("save_trees", {
  tree_list <- list(
    noninteractive = ggplot2::ggplot(),
    with_interactivity_data = ggplot2::ggplot(),
    with_heatmap = ggplot2::ggplot(),
    interactive = ggplot2::ggplot()
  )
  branch_col <- "logistic_growth_rate"
  default_base_name <- glue::glue("tree-{branch_col}")
  dated_base_name <- glue::glue("tree-{branch_col}-{Sys.Date()}")

  it("creates .svg / .rds / .html files - default: without the current date", {
    td <- withr::local_tempdir(pattern = "no-date-svg")

    filetypes <- list(
      expected = c("noninteractive", "interactive_html", "interactive_rds"),
      unexpected = character(0)
    )

    default_files <- get_default_filepaths(td, default_base_name)
    expected_files <- default_files[filetypes$expected]
    unexpected_files <- default_files[filetypes$unexpected]

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100
    )

    expect_mapequal(created_files, expected_files)
    expect_true(all(file.exists(expected_files)))
  })

  it("creates interactive html plot 'only' based on output_format argument", {
    td <- withr::local_tempdir(pattern = "no-date-svg")

    filetypes <- list(
      expected = c("noninteractive", "interactive_html"),
      unexpected = "interactive_rds"
    )

    default_files <- get_default_filepaths(td, default_base_name)
    expected_files <- default_files[filetypes$expected]
    unexpected_files <- default_files[filetypes$unexpected]

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100,
      output_format = "html"
    )

    expect_mapequal(created_files, expected_files)
    expect_true(all(file.exists(expected_files)))
    expect_false(any(file.exists(unexpected_files)))
  })

  it("creates interactive rds plot 'only' based on output_format argument", {
    td <- withr::local_tempdir(pattern = "no-date-svg")

    filetypes <- list(
      expected = c("noninteractive", "interactive_rds"),
      unexpected = "interactive_html"
    )

    default_files <- get_default_filepaths(td, default_base_name)
    expected_files <- default_files[filetypes$expected]
    unexpected_files <- default_files[filetypes$unexpected]

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100,
      output_format = "rds"
    )

    expect_mapequal(created_files, expected_files)
    expect_true(all(file.exists(expected_files)))
    expect_false(any(file.exists(unexpected_files)))
  })

  it("creates .svg / .rds / .html files - optionally including the current date", {
    td <- withr::local_tempdir(pattern = "noninteractive-svg")

    filetypes <- list(
      expected = c("noninteractive", "interactive_html", "interactive_rds"),
      unexpected = character(0)
    )

    default_files <- get_default_filepaths(td, dated_base_name)
    expected_files <- default_files[filetypes$expected]
    unexpected_files <- default_files[filetypes$unexpected]

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100, include_date = TRUE
    )

    expect_mapequal(created_files, expected_files)
    expect_true(all(file.exists(expected_files)))
    expect_false(any(file.exists(unexpected_files)))
  })
})
