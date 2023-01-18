# tfpscanner 0.2.2 _2023-01-18_

- Interactive tree view plots can be saved in either `html` (as an htmlwidget) or `rds` (as a 
  ggtree object) files
- Cluster SINA plots and interactive tree view plots are saved as _both_ `html` and `rds` by
  default

# tfpscanner 0.2.1 _2022-12-06_

- "Cluster sina plot"s can be saved in either `html` (as an htmlwidget) or `rds` (as a ggplot2
  object) files

# tfpscanner 0.2.0 _2022-11-09_

- Breaking change: `htmlwidget`s are no longer saved to file by `treeview()`
- Breaking change: date of creation is absent from the "tree-..." files created by `treeview()`
- Breaking change: `create_trees()` no longer creates `htmlwidget` tree-views
- `treeview()` now saves `ggtree` objects to `tree-xyz.rds` files, these can be converted to
  `htmlwidgets` using `ggiraph::girafe`
- `save_trees()` returns the file paths for any files that it creates

# tfpscanner 0.1.7 _2022-11-09_

- Separate the function for creating and saving tree plots from the `treeview` function

# tfpscanner 0.1.6 _2022-11-09_

- Split up the `.plot_tree` function into smaller components

# tfpscanner 0.1.5 _2022-11-07_

- Separate the function for plotting `sina` cluster data from the `treeview` function

# tfpscanner 0.1.4 _2022-11-07_

- Add `testthat` and `pre-commit` infrastructure

# tfpscanner 0.1.3 _2022-11-01_

- Integrate changes from mrc-ide repo with those in Jumping Rivers repo
- Styling with styler
- Fix some undocumented parameters

# tfpscanner 0.1.2 _2022-10-03_

- Revise sina plot
- Output sina for more variables
- Move doMPI to a suggested dependency

# tfpscanner 0.1.1 _2022-09-26_

- Add linting
- Add namespacing
- Pass R CMD Check
- Add system requirements

# tfpscanner 0.1.0 _2022-07-03_

- Fix description file
- Remove`library()` calls
