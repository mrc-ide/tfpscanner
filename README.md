# {tfpscanner}

## Installation

```
remotes::install_github("mrc-ide/tfpscanner")
```

## Use with {tfpbrowser}

{tfpscanner} is able to produce data that can be used in the {tfpbrowser} application. That app
must be pointed to a data-directory that contains the tfpscanner-generated files (see the notes in
{tfpbrowser} for details).

The structure of the data-directory used by tfpbrowser is as follows:

```bash
<data-root>
  |- mutations
  |  |- all_mutations.csv
  |  |- defining_mutations.csv
  |- scanner_output
  |  |- <clusterID1>
  |  |  |- cocirculating_lineages.csv
  |  |  |- lineage_composition.csv
  |  |  |- regional_composition.csv
  |  |  |- sequences.csv
  |  |  |- summary.csv
  |  |- <clusterID2>
  |  |  |- <as for clusterID1>
  |  |- ...
  |- sequences
  |  |- all_sequences.csv
  |- treeview
    |- node_lookup
    |  |- sina-logistic_growth_rate.csv
    |  |- tree-logistic_growth_rate.csv
    |  |- (a .csv for each .rds file in ./treeview)
    |- sina-logistic_growth_rate.rds
    |- tree-logistic_growth_rate.rds
    |- tree-mutations.rds
    |- tree-sequences.rds
    |- (various other .rds files containing trees that are presented by tfpbrowser)
```

The files in that data directory can be populated using `tfpscanner::create_browser_data()`.
To run that function requires that a few files are in place:

- A directory into which the tfpbrowser files are to be added has been added to the file-system
- That data directory contains a `scanner_output` directory (this is typically generated by
  `tfpscanner::tfpscan()`)
- A scanner environment file (typically having the name `scanner_output/scanner-env-<date>.rds`), as
  generated by `tfpscanner::tfpscan()` is present in the data-directory

Suppose:

- your data-directory was `./tfpbrowser_files/`;
- the files generated by `tfpscan()` had been added to `./tfpbrowser_files/scanner_output/`; and
- the scanner-environment file was `./tfpbrowser_files/scanner_output/scanner-env-2023-07-05.rds`.

Then, to set up all the remaining files required by tfpbrowser, you would make the following
function call:

```
tfpscanner::create_browser_files(
  e0 = "tfpbrowser_files/scanner_output/scanner-env-2023-07-05.rds",
  output_dir = "tfpbrowser_files",
  [any additional arguments to be passed on to `tfpscanner::treeview()`]
)
```
