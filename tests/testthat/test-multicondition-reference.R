read_reference_csv <- function(path) {
  read.csv(path, row.names = 1, check.names = FALSE)
}

run_multicondition_analysis <- function(trans_model) {
  pkg_root <- normalizePath(
    file.path(testthat::test_path("..", "..")),
    winslash = "/",
    mustWork = TRUE
  )

  counts_path <- file.path(pkg_root, "inst", "extdata", "ballinger_counts.csv")
  metadata_path <- file.path(pkg_root, "inst", "extdata", "ballinger_metadata.csv")

  counts <- as.matrix(read.csv(counts_path, row.names = 1, check.names = FALSE))
  metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)
  if ("Sample" %in% colnames(metadata)) {
    rownames(metadata) <- metadata$Sample
  }

  fit_obj <- new(
    "fitObject",
    counts = counts,
    metadata = metadata,
    trans_model = trans_model,
    fields_to_test = c("Tissue", "Temperature")
  )
  fit_obj <- fit_edgeR(fit_obj)

  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::normLibSizes(y)
  y <- edgeR::estimateDisp(y)

  list(
    design_matrix = as.data.frame(fit_obj@design_matrix_full, check.names = FALSE),
    weights = as.data.frame(fit_obj@weights, check.names = FALSE),
    tagwise_dispersion = data.frame(
      Genes = rownames(counts),
      tagwise_dispersion = y$tagwise.dispersion,
      row.names = rownames(counts),
      check.names = FALSE
    ),
    raw_pvals = fit_obj@raw_pvals,
    fdrs = fit_obj@BH_FDRs
  )
}

test_that("Multicondition example remains stable across trans-model fits", {
  pkg_root <- normalizePath(
    file.path(testthat::test_path("..", "..")),
    winslash = "/",
    mustWork = TRUE
  )
  fixture_dir <- file.path(pkg_root, "tests", "fixtures", "multicondition")

  for (trans_model in c("log_additive", "dominant")) {
    actual <- run_multicondition_analysis(trans_model)

    expected_design <- read_reference_csv(
      file.path(fixture_dir, paste0(trans_model, "_design_matrix.csv"))
    )
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

    expect_true(any(grepl("^beta_cis\\*tissue-", colnames(actual$design_matrix))))
    expect_true(any(grepl("^beta_trans\\*temperature-", colnames(actual$design_matrix))))
    expect_true("null: no Tissue cis" %in% colnames(actual$raw_pvals))
    expect_true("null: no Tissue trans" %in% colnames(actual$raw_pvals))
    expect_true("null: no Temperature cis" %in% colnames(actual$raw_pvals))
    expect_true("null: no Temperature trans" %in% colnames(actual$raw_pvals))

    expect_equal(actual$design_matrix, expected_design, tolerance = 1e-12)
    expect_equal(actual$weights, expected_weights, tolerance = 1e-12)
    expect_equal(actual$tagwise_dispersion, expected_dispersion, tolerance = 1e-12)
    expect_equal(actual$raw_pvals, expected_raw_pvals, tolerance = 1e-12)
    expect_equal(actual$fdrs, expected_fdrs, tolerance = 1e-12)
  }
})
