describe("save_trees", {
  tree_list <- list(
    noninteractive = ggplot2::ggplot(),
    with_interactivity_data = ggplot2::ggplot(),
    with_heatmap = ggplot2::ggplot(),
    interactive = ggplot2::ggplot()
  )

  it("creates an .svg from the 'noninteractive' entry", {
    td <- withr::local_tempdir(pattern = "noninteractive-svg")
    branch_col <- "logistic_growth_rate"
    output_path <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.svg"))

    save_trees(tree_list, branch_col = branch_col, output_dir = td, n_leaves = 100)

    expect_true(file.exists(output_path))
  })
})
