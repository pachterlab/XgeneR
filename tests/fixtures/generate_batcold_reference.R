script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
pkg_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

sys.source(file.path(pkg_root, "R", "XgeneR_class.R"), envir = globalenv())
sys.source(file.path(pkg_root, "R", "XgeneR_plotting.R"), envir = globalenv())

counts_path <- file.path(pkg_root, "inst", "extdata", "BATcold_ballinger_counts.csv")
metadata_path <- file.path(pkg_root, "inst", "extdata", "BATcold_ballinger_metadata.csv")
output_dir <- file.path(pkg_root, "tests", "fixtures", "batcold")

counts <- as.matrix(read.csv(counts_path, row.names = 1, check.names = FALSE))
metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)

collect_outputs <- function(trans_model) {
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

write_outputs <- function(trans_model, outputs) {
  prefix <- file.path(output_dir, trans_model)

  write.csv(outputs$weights, paste0(prefix, "_weights.csv"), quote = FALSE)
  write.csv(
    outputs$tagwise_dispersion,
    paste0(prefix, "_tagwise_dispersion.csv"),
    quote = FALSE
  )
  write.csv(outputs$raw_pvals, paste0(prefix, "_raw_pvals.csv"), quote = FALSE)
  write.csv(outputs$fdrs, paste0(prefix, "_fdrs.csv"), quote = FALSE)
  write.csv(
    outputs$proportion_cis,
    paste0(prefix, "_proportion_cis.csv"),
    quote = FALSE,
    row.names = FALSE
  )
}

for (trans_model in c("log_additive", "dominant")) {
  outputs <- collect_outputs(trans_model)
  write_outputs(trans_model, outputs)
  message("Wrote BATcold fixtures for ", trans_model)
}
