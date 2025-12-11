setClassUnion("listOrNULL", c("list", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
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
    covariate_cols = "characterOrNULL",  # vector of covariates that MUST be columns in metadata
    fields_to_test = "characterOrNULL",        # vector of fields to test that MUST be columns in metadata
    ref = "listOrNULL",                  # named list of reference levels for each field in field_to_test
    higher_order_interactions = "listOrNULL",  # extra higher order interactions to test
    design = "ANY",                      # defined upont initialization given fields_to_test and higher_order_interactions
    metadata_reformatted = "ANY",        # defined upon initialization
    design_matrix_full = "ANY",          # defined upon initialization given design
    contrast_vectors = "listOrNULL",     # defined upon initialization given design and fields to test
    raw_pvals = "ANY",                   # defined post fit
    BH_FDRs = "ANY",                     # defined post fit
    weights = "ANY"                     # defined post fit
    
  ),
  prototype = list(
    fields_to_test = NULL,
    ref = NULL,
    raw_pvals = NULL,
    weights = NULL,
    metadata_reformatted = NULL,
    design_matrix_full = NULL,
    contrast_vectors = NULL,
    BH_FDRs = NULL
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
convertMetadata <- function(metadata, ref) {
  stopifnot("Allele" %in% colnames(metadata))

  n_samples <- nrow(metadata)
  conditions <- setdiff(colnames(metadata), "Allele")
  
  metadata_reformatted <- data.frame(
    Samples = rownames(metadata),
    Allele = metadata$Allele,
    stringsAsFactors = FALSE
  )
  
  for (cond in conditions) {
      categories <- unique(metadata[[cond]])

      # Drop the first category to make design matrix identifiable
      categories_to_encode <- categories[-1]

      for (cat in categories_to_encode) {
        entries <- ifelse(metadata[[cond]] == cat, 1, 0)
        metadata_reformatted[[cat]] <- entries
      }
    }
  
  rownames(metadata_reformatted) <- metadata_reformatted$Samples
  metadata_reformatted$Samples <- NULL
  
  return(metadata_reformatted)
}

getDesignVector <- function(metadata,
                            covariate_cols,
                            fields_to_test,
                            ref = NULL,
                            higher_order_interactions = NULL,
                            interaction_designator = "*") {
  
  # Initialize design vector and reference categories
  design <- list("Reg")
  if (!is.null(fields_to_test)) {
      ref_categories <- character(length(fields_to_test))
      names(ref_categories) <- fields_to_test

      for (f in fields_to_test) {
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
  } else {ref_categories<-NULL}
  
  # Append higher order interactions if provided
  if (!is.null(higher_order_interactions)) {
    design <- append(design, higher_order_interactions)
  }
  
  return(list(design = design, ref = ref_categories))
}

getContrastVectors <- function(metadata, covariate_cols, fields_to_test, weight_names, interaction_designator = "*") {
  contrast_vector_list <- list()
  
  if (!is.null(fields_to_test)) {
      # Get unique categories for each field
      unique_categories <- lapply(fields_to_test, function(f) unique(metadata[[f]]))
      names(unique_categories) <- fields_to_test

      # Build all combinations
      combos <- expand.grid(unique_categories, stringsAsFactors = FALSE)
      combo_labels <- apply(combos, 1, function(row) paste(row, collapse = interaction_designator))

      # Collect all levels
      all_fields <- unlist(unique_categories, use.names = FALSE)

      # Error checking
      forbidden <- c("beta_cis", "beta_trans1", "beta_trans2")
      if (any(all_fields %in% forbidden)) {
        stop("Field name in beta_cis, beta_trans1, or beta_trans2. Please change.")
      }

      n_weights <- length(weight_names)

      # Build contrast vectors
      for (i in seq_along(combo_labels)) {
        combo <- combo_labels[i]
        combo_split <- unlist(strsplit(combo, interaction_designator,fixed=TRUE))
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
  } else  {
        n_covariate_weights <- 0

        if (!is.null(covariate_cols)) {
          for (col in covariate_cols) {
            values <- metadata[[col]]
            n_covariate_weights <- n_covariate_weights + (length(unique(values)) - 1)
          }
        }

        n_extra <- n_covariate_weights
        extra_zeros <- rep(0, n_extra)

        # core contrasts + zeros for covariates
        contrast_vector_list[[paste0("null: no cis")]]   <- c(0, 1, 0, extra_zeros)
        contrast_vector_list[[paste0("null: no trans")]] <- c(0, 0, 1, extra_zeros)
      }

      return(contrast_vector_list)
    }


.buildDesignMatrix <- function(metadata_reformatted,
                               covariate_cols,
                               fields_to_test,
                               design,
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
#   beta_trans1 <- as.integer(P1 | H1 | H2)
#   beta_trans1 <- as.integer(P2 | H1 | H2)
#   beta_trans1 <- ifelse(P1, 1,
#                       ifelse(H1 | H2, 0.5, 0))
  beta_trans <- ifelse(P1, 1,
                      ifelse(H1 | H2, 0.5, 0))

  design_full <- cbind(design_full, beta_cis, beta_trans)
  weight_names <- c(weight_names, "beta_cis", "beta_trans")
    
    
  if (is.null(fields_to_test) && !is.null(covariate_cols)) {
          covariate_matrix <- metadata_reformatted[, setdiff(colnames(metadata_reformatted), "Allele"), drop=FALSE]
          design_full <- cbind(design_full, covariate_matrix)
          weight_names <- c(weight_names, colnames(covariate_matrix))
        }

  if (!is.null(fields_to_test)) { 
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
  }

  rownames(design_full) <- rownames(metadata_reformatted)
  colnames(design_full) <- weight_names

  return(design_full)
}

# object methods                          
setMethod("initialize", "fitObject", function(.Object, counts, metadata, covariate_cols = NULL, fields_to_test = NULL,
                                              ref = NULL, higher_order_interactions = NULL) {
  .Object@counts <- counts
  .Object@covariate_cols <- covariate_cols 
  .Object@fields_to_test <- fields_to_test  
  .Object@ref <- ref
  .Object@higher_order_interactions <- higher_order_interactions
  
  cols_to_include <- c("Allele")
  if (!is.null(fields_to_test)) cols_to_include <- c(cols_to_include, fields_to_test)
  if (!is.null(covariate_cols)) cols_to_include <- c(cols_to_include, covariate_cols)
  
  # Subset raw metadata before converting
  metadata <- metadata[, intersect(cols_to_include, colnames(metadata)), drop = FALSE]
    
  .Object@metadata <- metadata

  # Compute reformatted metadata and design matrix
  metadata_ref <- convertMetadata(metadata, ref)
  result_designVector <- getDesignVector(metadata, covariate_cols, fields_to_test, ref, higher_order_interactions)
  design <- result_designVector[["design"]]
  ref <- result_designVector[["ref"]]
  design_matrix <- .buildDesignMatrix(metadata_ref, covariate_cols, fields_to_test, design)
  weight_names <- colnames(design_matrix)
  contrast_vectors <- getContrastVectors(metadata, covariate_cols, fields_to_test, weight_names)

  .Object@contrast_vectors <- contrast_vectors
  .Object@design <- design
  .Object@metadata_reformatted <- metadata_ref
  .Object@design_matrix_full <- design_matrix

  validObject(.Object)
  .Object
})

                              
get_fdrs <- function(pvals) {
  num_test <- length(pvals)
  # Get order of p-values
  order_idx <- order(pvals)
  # Sort p-values
  sorted_p <- pvals[order_idx]
  # BH adjustment formula
  fdr_sorted <- (seq_along(sorted_p) / num_test) * sorted_p
  # Ensure monotonicity of adjusted p-values
  fdr_sorted <- cummin(rev(fdr_sorted))
  fdr_sorted <- rev(fdr_sorted)
  # Put back in original order
  fdr <- numeric(num_test)
  fdr[order_idx] <- fdr_sorted
  return(fdr)
}

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
#   y <- edgeR::calcNormFactors(y)
  y <- edgeR::normLibSizes(y)
  y <- edgeR::estimateDisp(y)
  fit <- edgeR::glmFit(y, design_matrix)

  # Test various hypotheses
  i = 0 
  for (contrast_name in names(contrast_vectors)) {
      i = i+1
      contrast_vector <- contrast_vectors[[contrast_name]]
      lrt <- edgeR::glmLRT(fit, contrast = contrast_vector)
      raw_pvals <- lrt$table$PValue
      raw_pval_list[[contrast_name]] <- raw_pvals
#       corrected_fdr_list[[contrast_name]] <- p.adjust(raw_pvals, method = "BH")  perhaps implement later
      corrected_fdr_list[[contrast_name]] <- get_fdrs(raw_pvals)
    }
  print(i)

  object@raw_pvals <- as.data.frame(raw_pval_list, row.names = gene_names, check.names = FALSE)
  object@BH_FDRs <- as.data.frame(corrected_fdr_list, row.names = gene_names, check.names = FALSE)
  object@weights <- coef(fit)  # ADD WEIGHT NAMES AND GENE NAMES 

  return(object)
})


    
