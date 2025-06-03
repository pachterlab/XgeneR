#' @title Plot p-value histograms nd false discovery rates for certain condition combination. 
#' @description
#' Plots histograms of raw p-values for a given experimental combination.
#'
#' @param fitObject A `fitObject` that has been fitted using `fit_edgeR`.
#' @param combo A string specifying the experimental combination (e.g., "Treatment*Time").
#'
#' @return A ggplot object with the p-value histograms.
#' @export
#' @import ggplot2
#' @import patchwork
plotPvalHistograms <- function(fitObject, combo) {
  # Construct keys for named list lookup
  cis_label <- paste0(combo, " null: no cis")
  trans_label <- paste0(combo, " null: no trans")
  
  # Extract values
  raw_cis <- fitObject@raw_pvals[[cis_label]]
  raw_trans <- fitObject@raw_pvals[[trans_label]]
  fdr_cis <- fitObject@BH_FDRs[[cis_label]]
  fdr_trans <- fitObject@BH_FDRs[[trans_label]]
  
  # Plot objects
  p1 <- ggplot2::ggplot(data.frame(pval = raw_cis), aes(x = pval)) +
    geom_histogram(fill = "steelblue", bins = 30) +
    ggtitle("Raw P Values") +
    xlab(cis_label) + ylab("Count")

  p2 <- ggplot2::ggplot(data.frame(pval = raw_trans), aes(x = pval)) +
    geom_histogram(fill = "steelblue", bins = 30) +
    ggtitle("") +
    xlab(trans_label) + ylab("Count")

  p3 <- ggplot2::ggplot(data.frame(fdr = fdr_cis), aes(x = fdr)) +
    geom_histogram(fill = "darkred", bins = 30) +
    ggtitle("Benjamini-Hochberg corrected FDRs") +
    xlab(cis_label) + ylab("Count")

  p4 <- ggplot2::ggplot(data.frame(fdr = fdr_trans), aes(x = fdr)) +
    geom_histogram(fill = "darkred", bins = 30) +
    ggtitle("") +
    xlab(trans_label) + ylab("Count")

  # Arrange in 2x2 layout
  (p1 | p2) / (p3 | p4)
}



#' @title Classify genes and optionally plot predicted transformed data with cis/trans regulatory assignments.
#' @description
#' Classifies genes into regulatory categories based on hypothesis test results and optionally plots the transformed assignment.
#'
#' @param fitObject A `fitObject` that has been fitted using `fit_edgeR`.
#' @param combo A string specifying the experimental combination to classify.
#' @param plot Logical, whether to plot the assignments (default: TRUE).
#' @param cell_size Number of genes per cell in the hex plot (default: 10000).
#' @param interaction_designator Character used to separate field names (default: "*").
#'
#' @return A data.frame with gene assignments and optionally a ggplot object if `plot=TRUE`.
#' @export
#' @import ggplot2
#' @import gridExtra
getAssignmentsAndPlot <- function(fitObject, combo, plot = TRUE, 
                                  cell_size = 10000, alpha = 0.05,
                                  interaction_designator = "*") {
  
  weight_names <- colnames(fitObject@design_matrix)
  fdrs_no_cis <- fitObject@BH_FDRs[[paste0(combo, " null: no cis")]]
  fdrs_no_trans <- fitObject@BH_FDRs[[paste0(combo, " null: no trans")]]
  
  # Construct unique categories
  unique_categories <- lapply(fitObject@fields_to_test, function(f) {
    unique(fitObject@metadata[[f]])
  })
  names(unique_categories) <- fitObject@fields_to_test
  
  all_fields <- unlist(unique_categories)
  combo_split <- strsplit(combo, interaction_designator)[[1]]
  bad_fields <- all_fields[!all_fields %in% combo_split]
  
  combo_X <- matrix(0, nrow = 4, ncol = length(weight_names))
  
  for (i in seq_along(weight_names)) {
    weight <- weight_names[i]
    keep <- !any(sapply(bad_fields, function(b) grepl(b, weight)))
    cis <- keep && grepl("cis", weight)
    trans1 <- keep && grepl("trans1", weight)
    trans2 <- keep && grepl("trans2", weight)
    
    if (cis) {
      combo_X[2, i] <- 1
      combo_X[4, i] <- 1
    }
    if (trans1) {
      combo_X[1, i] <- 1
      combo_X[3, i] <- 1
      combo_X[4, i] <- 1
    }
    if (trans2) {
      combo_X[2, i] <- 1
      combo_X[3, i] <- 1
      combo_X[4, i] <- 1
    }
    if (keep && !cis && !trans1 && !trans2) {
      combo_X[, i] <- 1
    }
  }
  
  num_test <- nrow(fitObject@coefficients)
  pred_counts <- matrix(NA, nrow = num_test, ncol = 4)
  for (i in 1:num_test) {
    pred_log <- fitObject@coefficients[i, ] %*% t(combo_X) + log(cell_size)
    pred_counts[i, ] <- exp(pred_log)
  }
  
  P1 <- pred_counts[, 1]
  P2 <- pred_counts[, 2]
  H1 <- pred_counts[, 3]
  H2 <- pred_counts[, 4]
  
  df <- data.frame(
    gene = fitObject@genes,
    P1 = P1, P2 = P2, H1 = H1, H2 = H2,
    Parlog2FC = log2(P1 / P2),
    Hyblog2FC = log2(H1 / H2),
    fdr_cis = fdrs_no_cis,
    fdr_trans = fdrs_no_trans
  )
  
  P <- df$Parlog2FC
  H <- df$Hyblog2FC
  delta <- P - H
  
  cis_index <- df$fdr_cis < alpha & df$fdr_trans > alpha
  trans_index <- df$fdr_cis > alpha & df$fdr_trans < alpha
  cis_plus_trans_index <- df$fdr_cis < alpha & df$fdr_trans < alpha & ((delta > 0 & H > 0) | (delta < 0 & H < 0))
  cis_x_trans_index <- df$fdr_cis < alpha & df$fdr_trans < alpha & ((delta <= 0 & H >= 0) | (delta >= 0 & H <= 0))
  
  colors <- rep("lightgray", nrow(df))
  colors[cis_index] <- "orangered"
  colors[trans_index] <- "royalblue"
  colors[cis_plus_trans_index] <- "skyblue"
  colors[cis_x_trans_index] <- "green"
  df$colors <- colors
  
  reg_assignment <- rep("conserved", nrow(df))
  reg_assignment[cis_index] <- "cis"
  reg_assignment[trans_index] <- "trans"
  reg_assignment[cis_plus_trans_index] <- "cis+trans"
  reg_assignment[cis_x_trans_index] <- "cisxtrans"
  df$reg_assignment <- reg_assignment
  
  # Polar transformation
  theta_scaled <- (2 / pi) * atan2(H, delta)
  R <- sqrt(delta^2 + H^2)
  cis_prop_reordered <- theta_scaled - 0.5
  cis_prop_reordered[cis_prop_reordered <= -1] <- cis_prop_reordered[cis_prop_reordered <= -1] + 2.0
  
  df$R <- R
  df$theta_scaled <- theta_scaled
  df$cis_prop_reordered <- cis_prop_reordered
  
  if (plot) {
    p1 <- ggplot2::ggplot(df, aes(x = Parlog2FC, y = Hyblog2FC, color = reg_assignment)) +
      geom_point(size = 0.8, alpha = 0.7) +
      geom_abline(intercept = 0, slope = 1, color = "orangered") +
      geom_abline(intercept = 0, slope = 2, color = "skyblue") +
      geom_abline(intercept = 0, slope = -2, color = "forestgreen") +
      geom_vline(xintercept = 0) + geom_hline(yintercept = 0, color = "darkblue") +
      theme_bw() + ggtitle("Parental vs Hybrid Fold Change")
    
    p2 <- ggplot2::ggplot(df, aes(x = delta, y = H, color = reg_assignment)) +
      geom_point(size = 0.8, alpha = 0.7) +
      geom_vline(xintercept = 0, color = "orangered") +
      geom_hline(yintercept = 0, color = "darkblue") +
      geom_abline(intercept = 0, slope = 1, color = "skyblue") +
      geom_abline(intercept = 0, slope = -1, color = "forestgreen") +
      theme_bw() + ggtitle("Delta vs Hybrid Fold Change")
    
    p3 <- ggplot2::ggplot(df, aes(x = cis_prop_reordered, y = P, color = reg_assignment)) +
      geom_point(size = 0.8, alpha = 0.7) +
      geom_vline(xintercept = c(-1, 0, 0.5, 1), color = c("forestgreen", "skyblue", "orangered", "forestgreen")) +
      geom_hline(yintercept = 0, color = "black") +
      theme_bw() + ggtitle("Polar Representation of Cis-Trans Balance")
    
    grid.arrange(p1, p2, p3, nrow = 3)
  }
  
  return(df)
}
