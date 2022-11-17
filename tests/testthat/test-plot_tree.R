describe("save_trees", {
  tree_list <- list(
    noninteractive = ggplot2::ggplot(),
    with_interactivity_data = ggplot2::ggplot(),
    with_heatmap = ggplot2::ggplot(),
    interactive = ggplot2::ggplot()
  )
  branch_col <- "logistic_growth_rate"

  it("creates .svg / .rds files - default: without the current date", {
    td <- withr::local_tempdir(pattern = "no-date-svg")
    output_rds <- file.path(td, glue::glue("tree-{branch_col}.rds"))
    output_svg <- file.path(td, glue::glue("tree-{branch_col}.svg"))

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100
    )

    expect_equal(
      created_files,
      c("noninteractive" = output_svg, "interactive" = output_rds)
    )
    expect_true(file.exists(output_svg))
    expect_true(file.exists(output_rds))
  })

  it("creates .svg / .rds files - optionally including the current date", {
    td <- withr::local_tempdir(pattern = "noninteractive-svg")
    output_rds <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.rds"))
    output_svg <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.svg"))

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100, include_date = TRUE
    )

    expect_equal(
      created_files,
      c("noninteractive" = output_svg, "interactive" = output_rds)
    )
    expect_true(file.exists(output_svg))
    expect_true(file.exists(output_rds))
  })
})
