#' @title fitObject S4 Class
#' @description
#' An S4 class to store input data and fit results for allele-specific expression analysis using edgeR.
#' 
#' @slot counts A matrix of raw counts with genes as rows and samples as columns.
#' @slot metadata A data.frame containing metadata for the samples. Must contain a column named 'Allele'.
#' @slot fields_to_test Character vector of column names in metadata to test in the model.
#' @slot ref Named list specifying the reference levels for each metadata field.
#' @slot higher_order_interactions Optional list of higher-order interaction terms to include in the design.
#' @slot design The design vector for model building.
#' @slot metadata_reformatted Metadata reformatted into binary indicator variables.
#' @slot design_matrix_full Full numeric design matrix for model fitting.
#' @slot contrast_vectors List of contrast vectors for hypothesis testing.
#' @slot raw_pvals Raw p-values from hypothesis tests.
#' @slot BH_FDRs FDR-adjusted p-values using Benjamini-Hochberg.
#' @slot weights Coefficients (weights) estimated from the model.
#' 
#' @export
setClass(
  "fitObject",
  slots = list(
    counts = "matrix",
    metadata = "data.frame",
    fields_to_test = "character",        # vector of fields to test that MUST be columns in metadata
    ref = "list",                        # named list of reference levels for each field in field_to_test
    higher_order_interactions = "listOrNULL",  # extra higher order interactions to test
    design = "ANY",                      # defined upont initialization given fields_to_test and higher_order_interactions
    metadata_reformatted = "ANY",        # defined upon initialization
    design_matrix_full = "ANY",          # defined upon initialization given design
    contrast_vectors = "listOrNULL",     # defined upon initialization given design and fields to test
    raw_pvals = "ANY",                   # defined post fit
    BH_FDRs = "ANY",                     # defined post fit
    weights = "ANY",                     # defined post fit
    
  ),
  prototype = list(
    ref = NULL
    raw_pvals = NULL,
    weights = NULL,
    metadata_reformatted = NULL,
    design_matrix_full = NULL,
    contrast_vectors = NULL,
  ),
  validity = function(object) {
    if (is.null(colnames(object@counts))) {
      return("Counts matrix must have column names (sample IDs).")
    }
    if (is.null(rownames(object@counts))) {
      return("Counts matrix must have row names (gene IDs).")
    }
    if (!"Allele" %in% colnames(object@metadata)) {
      return("Metadata must include a column named 'Allele'.")
    }
    if (!all(colnames(object@counts) %in% rownames(object@metadata))) {
      return("All sample IDs in counts must be present in metadata row names.")
    }
    TRUE
  }
)

# Helper functions
convertMetadata <- function(object, ref) {
  metadata <- object@metadata
  stopifnot("Allele" %in% colnames(metadata))

  n_samples <- nrow(metadata)
  conditions <- setdiff(colnames(metadata), "Allele")
  
  metadata_reformatted <- data.frame(
    Samples = rownames(metadata),
    Allele = metadata$Allele,
    stringsAsFactors = FALSE
  )
  
  for (cond in conditions) {
    ref_cond <- ref[[cond]]
    categories <- setdiff(unique(metadata[[cond]]), ref_cond)
    
    for (cat in categories) {
      entries <- ifelse(metadata[[cond]] == cat, 1, 0)
      metadata_reformatted[[cat]] <- entries
    }
  }
  
  rownames(metadata_reformatted) <- metadata_reformatted$Samples
  metadata_reformatted$Samples <- NULL
  
  return(metadata_reformatted)
}

getDesignVector <- function(metadata,
                            fields,
                            ref = NULL,
                            higher_order_interactions = NULL,
                            interaction_designator = "*") {
  
  # Initialize design vector and reference categories
  design <- list("Reg")
  ref_categories <- character(length(fields))
  names(ref_categories) <- fields
  
  for (f in fields) {
    categories <- unique(metadata[[f]])
    
    if (!is.null(ref) && f %in% names(ref)) {
      category_to_remove <- ref[[f]]
    } else {
      category_to_remove <- categories[1]
    }
    
    # Remove reference from the category list
    categories <- setdiff(categories, category_to_remove)
    ref_categories[[f]] <- category_to_remove
    
    for (c in categories) {
      design <- append(design, list(paste0("Reg", interaction_designator, c)))
      design <- append(design, list(c))
    }
  }
  
  # Append higher order interactions if provided
  if (!is.null(higher_order_interactions)) {
    design <- append(design, higher_order_interactions)
  }
  
  return(list(design = design, ref = ref_categories))
}

getContrastVectors <- function(metadata, fields_to_test, weight_names, interaction_designator = "*") {
  contrast_vector_list <- list()
  
  # Step 1: get unique categories for each field
  unique_categories <- lapply(fields_to_test, function(f) unique(metadata[[f]]))
  names(unique_categories) <- fields_to_test
  
  # Step 2: build all combinations
  combos <- expand.grid(unique_categories, stringsAsFactors = FALSE)
  combo_labels <- apply(combos, 1, function(row) paste(row, collapse = interaction_designator))
  
  # Step 3: collect all levels
  all_fields <- unlist(unique_categories, use.names = FALSE)
  
  # Step 4: error checking
  forbidden <- c("beta_cis", "beta_trans1", "beta_trans2")
  if (any(all_fields %in% forbidden)) {
    stop("Field name in beta_cis, beta_trans1, or beta_trans2. Please change.")
  }
  
  n_weights <- length(weight_names)
  
  # Step 5: build contrast vectors
  for (i in seq_along(combo_labels)) {
    combo <- combo_labels[i]
    combo_split <- unlist(strsplit(combo, interaction_designator))
    bad_fields <- setdiff(all_fields, combo_split)
    
    # cis contrast
    cis_contrast <- numeric(n_weights)
    for (j in seq_along(weight_names)) {
      weight <- weight_names[j]
      keep <- !any(sapply(bad_fields, function(f) grepl(f, weight)))
      if (keep && grepl("cis", weight)) {
        cis_contrast[j] <- 1
      }
    }
    
    # trans contrast
    trans_contrast <- numeric(n_weights)
    for (j in seq_along(weight_names)) {
      weight <- weight_names[j]
      keep <- !any(sapply(bad_fields, function(f) grepl(f, weight)))
      if (keep && grepl("trans1", weight)) {
        trans_contrast[j] <- 1
      }
      if (keep && grepl("trans2", weight)) {
        trans_contrast[j] <- -1
      }
    }
    
    contrast_vector_list[[paste0(combo, " null: no cis")]] <- cis_contrast
    contrast_vector_list[[paste0(combo, " null: no trans")]] <- trans_contrast
  }
  
  return(contrast_vector_list)
}


.buildDesignMatrix <- function(metadata_reformatted, design,
                               reg_designator = "Reg",
                               interact_designator = "*") {
  stopifnot(design[1] == reg_designator)
  stopifnot("Allele" %in% colnames(metadata_reformatted))
  
  n_samples <- nrow(metadata_reformatted)
  design_full <- matrix(0, nrow = n_samples, ncol = 0)
  weight_names <- c()

  # Add intercept
  design_full <- cbind(design_full, rep(1, n_samples))
  weight_names <- c(weight_names, "Intercept")

  P1 <- metadata_reformatted$Allele == "P1"
  P2 <- metadata_reformatted$Allele == "P2"
  H1 <- metadata_reformatted$Allele == "H1"
  H2 <- metadata_reformatted$Allele == "H2"

  beta_cis <- as.integer(P2 | H2)
  beta_trans1 <- as.integer(P1 | H1 | H2)
  beta_trans2 <- as.integer(P2 | H1 | H2)

  design_full <- cbind(design_full, beta_cis, beta_trans1, beta_trans2)
  weight_names <- c(weight_names, "beta_cis", "beta_trans1", "beta_trans2")

  for (des in design[-1]) {
    if (grepl(paste0(reg_designator, interact_designator), des, fixed = TRUE)) {
      conds <- strsplit(des, split = interact_designator, fixed = TRUE)[[1]][-1]
      cond_filt <- rep(TRUE, n_samples)
      for (cond in conds) {
        cond_filt <- cond_filt & (metadata_reformatted[[cond]] == 1)
      }
      name <- paste(conds, collapse = interact_designator)
      
      design_full <- cbind(
        design_full,
        as.integer((P2 | H2) & cond_filt),
        as.integer((P1 | H1 | H2) & cond_filt),
        as.integer((P2 | H1 | H2) & cond_filt)
      )
      
      weight_names <- c(weight_names,
                        paste0("beta_cis*", name),
                        paste0("beta_trans1*", name),
                        paste0("beta_trans2*", name))
    } else {
      conds <- strsplit(des, split = interact_designator, fixed = TRUE)[[1]]
      cond_filt <- rep(TRUE, n_samples)
      for (cond in conds) {
        cond_filt <- cond_filt & (metadata_reformatted[[cond]] == 1)
      }
      name <- paste(conds, collapse = interact_designator)
      
      design_full <- cbind(design_full, as.integer(cond_filt))
      weight_names <- c(weight_names, paste0("beta_", name))
    }
  }

  rownames(design_full) <- rownames(metadata_reformatted)
  colnames(design_full) <- weight_names

  return(design_full)
}

# object methods                          
setMethod("initialize", "fitObject", function(.Object, counts, metadata, fields_to_test,
                                              ref = NULL, higher_order_interactions = NULL) {
  .Object@counts <- counts
  .Object@metadata <- metadata
  .Object@fields_to_test <- fields_to_test  
  .Object@ref <- ref
  .Object@higher_order_interactions <- higher_order_interactions
  
  # Compute reformatted metadata and design matrix
  metadata_ref <- convertMetadata(.Object, ref)
  result_designVector <- getDesignVector(.Object, metadata, fields_to_test, ref, higher_order_interactions)
  design <- result_designVector[["design"]]
  ref <- result_designVector[["ref"]]
  design_matrix <- .buildDesignMatrix(metadata_ref, design)
  weight_names <- colnames(design_matrix)
  contrast_vectors <- getContrastVectors(metadata, fields_to_test, weight_names)
  
  .Object@contrast_vectors <- contrast_vectors
  .Object@design <- design
  .Object@metadata_reformatted <- metadata_ref
  .Object@design_matrix_full <- design_matrix

  validObject(.Object)
  .Object
})


setGeneric("fit_edgeR", function(object, ...) standardGeneric("fit_edgeR"))
#' @title Fit GLM on parent and hybrid crosses using edgeR.
#' @description
#' Fits a negative binomial GLM using edgeR for a given `fitObject`, testing for cis and trans effects between homozygous parents and their heterozygous crosses.
#'
#' @param object A `fitObject` instance containing counts, metadata, and design.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `fitObject` with updated slots for raw p-values, adjusted FDRs, and model weights.
#' @export
#' @import edgeR           
setMethod("fit_edgeR", "fitObject", function(object, ...) {
  counts <- object@counts
  design_matrix <- object@design_matrix_full
  contrast_vectors <- object@contrast_vectors
  gene_names <- rownames(counts)

  raw_pval_list <- list(Genes = gene_names)
  corrected_fdr_list <- list(Genes = gene_names)

  # Run edgeR pipeline
  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmFit(y, design)

  # Test various hypotheses
  for (contrast_name in names(contrast_vectors)) {
      contrast_vector <- contrast_vectors[[contrast_name]]

      lrt <- edgeR::glmLRT(fit, contrast = contrast_vector)
      raw_pvals <- lrt$table$PValue
      raw_pval_list[[contrast_name]] <- raw_pvals
      corrected_fdr_list[[contrast_name]] <- p.adjust(raw_pvals, method = "BH")
    }

  object@raw_pvals <- as.data.frame(raw_pval_list, row.names = gene_names)
  object@BH_FDRs <- as.data.frame(corrected_fdr_list, row.names = gene_names)
  object@weights <- coef(fit)  # ADD WEIGHT NAMES AND GENE NAMES 

  return(object)
})


    
