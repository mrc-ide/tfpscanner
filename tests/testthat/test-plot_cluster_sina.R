get_pldf <- function() {
  data.frame(
    isTip = c(TRUE, FALSE, TRUE, FALSE),
    logistic_growth_rate = c(-0.1, 0, 0.5, -10)
  )
}

describe("format_cluster_sina_data", {
  pldf <- get_pldf()
  varx <- "logistic_growth_rate"

  sina_data <- format_cluster_sina_data(pldf, varx = varx)

  it("returns an entry for each row that isTip", {
    expect_true(all(sina_data$isTip))
    expect_equal(
      sina_data[[varx]],
      pldf$logistic_growth_rate[pldf$isTip]
    )
  })
  it("adds a 'varx' and 'y' column that matches a named column", {
    expect_equal(
      sina_data$varx,
      sina_data[[varx]]
    )
    expect_equal(
      sina_data$y,
      sina_data[[varx]]
    )
  })
  it("adds a 'mutation_lineage' column", {
    expect_true("mutation_lineage" %in% colnames(sina_data))
  })
})

describe("get_mutation_and_lineage", {
  df <- data.frame(
    lineage = c(
      "AY.43", "AY.43|AY.9", "B.1.351|None"
    ),
    allmuts = paste(
      "all mutations:\n",
      c("E:P71L", "S:A222V, N:T205I", "S:E484K"),
      "\n\n",
      sep = ""
    )
  )

  it("returns a '<mutation>_<lineage>' string", {
    ml <- get_mutation_and_lineage(df, mut_regexp = "S:A222V", lineage_regexp = "AY\\.43")
    expect_equal(
      ml,
      c("_AY\\.43", "S:A222V_AY\\.43", "_")
    )
  })
  it("is just an underscore if neither regexp pattern is specified", {
    ml <- get_mutation_and_lineage(df, mut_regexp = NULL, lineage_regexp = NULL)
    expect_equal(
      ml,
      c("_", "_", "_")
    )
  })
  it("includes a dot-separated mutation string", {
    ml <- get_mutation_and_lineage(
      df,
      mut_regexp = c("E:P71L", "S:A222V", "N:T205I"),
      lineage_regexp = NULL
    )
    expect_equal(
      ml,
      c("E:P71L.._", ".S:A222V.N:T205I_", ".._")
    )
  })
  it("includes a dot-separated lineage string", {
    ml <- get_mutation_and_lineage(
      df,
      mut_regexp = NULL,
      lineage_regexp = c("AY\\.9", "AY\\.43")
    )
    expect_equal(
      ml,
      c("_.AY\\.43", "_AY\\.9.AY\\.43", "_.")
    )
  })
})
