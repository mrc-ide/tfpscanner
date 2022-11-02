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
