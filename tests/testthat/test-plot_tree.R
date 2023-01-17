describe("save_trees", {
  tree_list <- list(
    noninteractive = ggplot2::ggplot(),
    with_interactivity_data = ggplot2::ggplot(),
    with_heatmap = ggplot2::ggplot(),
    interactive = ggplot2::ggplot()
  )
  branch_col <- "logistic_growth_rate"

  it("creates .svg / .rds / .html files - default: without the current date", {
    td <- withr::local_tempdir(pattern = "no-date-svg")
    output_rds <- file.path(td, glue::glue("tree-{branch_col}.rds"))
    output_html <- file.path(td, glue::glue("tree-{branch_col}.html"))
    output_svg <- file.path(td, glue::glue("tree-{branch_col}.svg"))

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100
    )

    expect_mapequal(
      created_files,
      c(
        "noninteractive" = output_svg,
        "interactive_rds" = output_rds,
        "interactive_html" = output_html
      )
    )
    expect_true(file.exists(output_svg))
    expect_true(file.exists(output_html))
    expect_true(file.exists(output_rds))
  })

  it("creates interactive html plot 'only' based on output_format argument", {
    td <- withr::local_tempdir(pattern = "no-date-svg")
    output_html <- file.path(td, glue::glue("tree-{branch_col}.html"))
    output_svg <- file.path(td, glue::glue("tree-{branch_col}.svg"))
    unexpected_rds <- file.path(td, glue::glue("tree-{branch_col}.rds"))

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100,
      output_format = "html"
    )

    expect_mapequal(
      created_files,
      c("noninteractive" = output_svg, "interactive_html" = output_html)
    )
    expect_true(file.exists(output_svg))
    expect_true(file.exists(output_html))
    expect_false(file.exists(unexpected_rds))
  })

  it("creates interactive rds plot 'only' based on output_format argument", {
    td <- withr::local_tempdir(pattern = "no-date-svg")
    output_rds <- file.path(td, glue::glue("tree-{branch_col}.rds"))
    output_svg <- file.path(td, glue::glue("tree-{branch_col}.svg"))
    unexpected_html <- file.path(td, glue::glue("tree-{branch_col}.html"))

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100,
      output_format = "rds"
    )

    expect_mapequal(
      created_files,
      c("noninteractive" = output_svg, "interactive_rds" = output_rds)
    )
    expect_true(file.exists(output_svg))
    expect_true(file.exists(output_rds))
    expect_false(file.exists(unexpected_html))
  })

  it("creates .svg / .rds / .html files - optionally including the current date", {
    td <- withr::local_tempdir(pattern = "noninteractive-svg")
    output_rds <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.rds"))
    output_html <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.html"))
    output_svg <- file.path(td, glue::glue("tree-{branch_col}-{Sys.Date()}.svg"))

    created_files <- save_trees(
      tree_list,
      branch_col = branch_col, output_dir = td, n_leaves = 100, include_date = TRUE
    )

    expect_equal(
      created_files,
      c(
        "noninteractive" = output_svg,
        "interactive_rds" = output_rds,
        "interactive_html" = output_html
      )
    )
    expect_true(file.exists(output_svg))
    expect_true(file.exists(output_rds))
    expect_true(file.exists(output_html))
  })
})
