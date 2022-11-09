describe("save_trees", {
  tree_list <- list(
    noninteractive = ggplot2::ggplot(),
    with_interactivity_data = ggplot2::ggplot(),
    with_heatmap = ggplot2::ggplot(),
    interactive = ggplot2::ggplot()
  )
  branch_col <- "logistic_growth_rate"

  it("creates an .svg from the 'noninteractive' entry - default: without the current date", {
    td <- withr::local_tempdir(pattern = "no-date-svg")
    output_path <- file.path(td, glue::glue("tree-{branch_col}.svg"))

    save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100
    )

    expect_true(file.exists(output_path))
  })

  it("creates an .rds from the 'interactive' entry - default: without the current date", {
    td <- withr::local_tempdir(pattern = "no-date-rds")
    output_path <- file.path(td, glue::glue("tree-{branch_col}.rds"))

    save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100
    )

    expect_true(file.exists(output_path))
  })

  it("creates an .svg from the 'noninteractive' entry - including the current date", {
    td <- withr::local_tempdir(pattern = "noninteractive-svg")
    output_path <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.svg"))

    save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100, include_date = TRUE
    )

    expect_true(file.exists(output_path))
  })

  it("creates an .rds from the 'interactive' entry - including the current date", {
    td <- withr::local_tempdir(pattern = "interactive-rds")
    output_path <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.rds"))

    save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100, include_date = TRUE
    )

    expect_true(file.exists(output_path))
  })
})
