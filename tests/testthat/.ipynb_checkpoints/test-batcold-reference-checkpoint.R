read_reference_csv <- function(path, use_row_names = TRUE) {
  if (use_row_names) {
    read.csv(path, row.names = 1, check.names = FALSE)
  } else {
    read.csv(path, check.names = FALSE)
  }
}

run_batcold_analysis <- function(trans_model) {
  pkg_root <- normalizePath(
    file.path(testthat::test_path("..", "..")),
    winslash = "/",
    mustWork = TRUE
  )

  counts_path <- file.path(pkg_root, "inst", "extdata", "BATcold_ballinger_counts.csv")
  metadata_path <- file.path(pkg_root, "inst", "extdata", "BATcold_ballinger_metadata.csv")

  counts <- as.matrix(read.csv(counts_path, row.names = 1, check.names = FALSE))
  metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

  fit_obj <- new(
    "fitObject",
    counts = counts,
    metadata = metadata,
    trans_model = trans_model,
    fields_to_test = NULL
  )
  fit_obj <- fit_edgeR(fit_obj)

  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::normLibSizes(y)
  y <- edgeR::estimateDisp(y)

  assignments <- getAssignmentsAndPlot(fit_obj, plot = FALSE)

  list(
    weights = as.data.frame(fit_obj@weights, check.names = FALSE),
    tagwise_dispersion = data.frame(
      Genes = rownames(counts),
      tagwise_dispersion = y$tagwise.dispersion,
      row.names = rownames(counts),
      check.names = FALSE
    ),
    raw_pvals = fit_obj@raw_pvals,
    fdrs = fit_obj@BH_FDRs,
    proportion_cis = assignments$df[, c("gene", "cis_prop")]
  )
}

test_that("BATcold example remains stable across trans-model fits", {
  pkg_root <- normalizePath(
    file.path(testthat::test_path("..", "..")),
    winslash = "/",
    mustWork = TRUE
  )
  fixture_dir <- file.path(pkg_root, "tests", "fixtures", "batcold")

  for (trans_model in c("log_additive", "dominant")) {
    actual <- run_batcold_analysis(trans_model)

    expected_weights <- read_reference_csv(
      file.path(fixture_dir, paste0(trans_model, "_weights.csv"))
    )
    expected_dispersion <- read_reference_csv(
      file.path(fixture_dir, paste0(trans_model, "_tagwise_dispersion.csv"))
    )
    expected_raw_pvals <- read_reference_csv(
      file.path(fixture_dir, paste0(trans_model, "_raw_pvals.csv"))
    )
    expected_fdrs <- read_reference_csv(
      file.path(fixture_dir, paste0(trans_model, "_fdrs.csv"))
    )
    expected_prop_cis <- read_reference_csv(
      file.path(fixture_dir, paste0(trans_model, "_proportion_cis.csv")),
      use_row_names = FALSE
    )

    expect_equal(actual$weights, expected_weights, tolerance = 1e-12)
    expect_equal(actual$tagwise_dispersion, expected_dispersion, tolerance = 1e-12)
    expect_equal(actual$raw_pvals, expected_raw_pvals, tolerance = 1e-12)
    expect_equal(actual$fdrs, expected_fdrs, tolerance = 1e-12)
    expect_equal(actual$proportion_cis, expected_prop_cis, tolerance = 1e-12)
  }
})
