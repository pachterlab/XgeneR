script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
pkg_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

sys.source(file.path(pkg_root, "R", "XgeneR_class.R"), envir = globalenv())

counts_path <- file.path(pkg_root, "inst", "extdata", "ballinger_counts.csv")
metadata_path <- file.path(pkg_root, "inst", "extdata", "ballinger_metadata.csv")
output_dir <- file.path(pkg_root, "tests", "fixtures", "multicondition")

counts <- as.matrix(read.csv(counts_path, row.names = 1, check.names = FALSE))
metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)
if ("Sample" %in% colnames(metadata)) {
  rownames(metadata) <- metadata$Sample
}

collect_outputs <- function(trans_model) {
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

write_outputs <- function(trans_model, outputs) {
  prefix <- file.path(output_dir, trans_model)
  write.csv(outputs$design_matrix, paste0(prefix, "_design_matrix.csv"), quote = FALSE)
  write.csv(outputs$weights, paste0(prefix, "_weights.csv"), quote = FALSE)
  write.csv(
    outputs$tagwise_dispersion,
    paste0(prefix, "_tagwise_dispersion.csv"),
    quote = FALSE
  )
  write.csv(outputs$raw_pvals, paste0(prefix, "_raw_pvals.csv"), quote = FALSE)
  write.csv(outputs$fdrs, paste0(prefix, "_fdrs.csv"), quote = FALSE)
}

for (trans_model in c("log_additive", "dominant")) {
  outputs <- collect_outputs(trans_model)
  write_outputs(trans_model, outputs)
  message("Wrote multicondition fixtures for ", trans_model)
}
