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
plotPvalHistograms <- function(fitObject, combo = NULL) {
  # Construct keys for named list lookup
  if (!is.null(combo)) {
      cis_label <- paste0(combo, " null: no cis")
      trans_label <- paste0(combo, " null: no trans")} else {
      cis_label <- "null: no cis"
      trans_label <- "null: no trans"
  }
  
  
  # Extract values
  raw_cis <- fitObject@raw_pvals[[cis_label]]
  raw_trans <- fitObject@raw_pvals[[trans_label]]
  fdr_cis <- fitObject@BH_FDRs[[cis_label]]
  fdr_trans <- fitObject@BH_FDRs[[trans_label]]
  
  # Plot objects
    p1 <- ggplot2::ggplot(data.frame(pval = raw_cis), ggplot2::aes(x = pval)) +
      ggplot2::geom_histogram(fill = "steelblue", bins = 30) +
      ggplot2::ggtitle("Raw P Values") +
      ggplot2::xlab(cis_label) + ggplot2::ylab("Count")

    p2 <- ggplot2::ggplot(data.frame(pval = raw_trans), ggplot2::aes(x = pval)) +
      ggplot2::geom_histogram(fill = "steelblue", bins = 30) +
      ggplot2::ggtitle("") +
      ggplot2::xlab(trans_label) + ggplot2::ylab("Count")

    p3 <- ggplot2::ggplot(data.frame(fdr = fdr_cis), ggplot2::aes(x = fdr)) +
      ggplot2::geom_histogram(fill = "darkred", bins = 30) +
      ggplot2::ggtitle("Benjamini-Hochberg corrected FDRs") +
      ggplot2::xlab(cis_label) + ggplot2::ylab("Count")

    p4 <- ggplot2::ggplot(data.frame(fdr = fdr_trans), ggplot2::aes(x = fdr)) +
      ggplot2::geom_histogram(fill = "darkred", bins = 30) +
      ggplot2::ggtitle("") +
      ggplot2::xlab(trans_label) + ggplot2::ylab("Count")


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
getAssignmentsAndPlot <- function(fitObject, combo = NULL, plot = TRUE, 
                                  cell_size = 10000, alpha = 0.05,
                                  interaction_designator = "*") {
  
  weight_names <- colnames(fitObject@design_matrix_full)
    
  if (!is.null(combo)) {
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
                        
  } else {
      fdrs_no_cis <- fitObject@BH_FDRs[["null: no cis"]]
      fdrs_no_trans <- fitObject@BH_FDRs[["null: no trans"]]
      
      combo_X <- matrix(0, nrow = 4, ncol = length(weight_names))
      # intercept
      combo_X[1,1] <- 1
      combo_X[2,1] <- 1
      combo_X[3,1] <- 1
      combo_X[4,1] <- 1
      # cis
      combo_X[2, 2] <- 1
      combo_X[4, 2] <- 1
      # trans1
      combo_X[1, 3] <- 1
      combo_X[3, 3] <- 1
      combo_X[4, 3] <- 1
      # trans2
      combo_X[2, 4] <- 1
      combo_X[3, 4] <- 1
      combo_X[4, 4] <- 1
  }
                        
  num_test <- nrow(fitObject@weights)
  pred_counts <- matrix(NA, nrow = num_test, ncol = 4)
  for (i in 1:num_test) {
    pred_log <- fitObject@weights[i, ] %*% t(combo_X) + log(cell_size)
    pred_counts[i, ] <- exp(pred_log)
  }
  
  P1 <- pred_counts[, 1]
  P2 <- pred_counts[, 2]
  H1 <- pred_counts[, 3]
  H2 <- pred_counts[, 4]                        
  
  df <- data.frame(
    gene = rownames(fitObject@counts),
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
  colors[cis_x_trans_index] <- "forestgreen"
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
  cis_prop_reordered[cis_prop_reordered > 1] = cis_prop_reordered[cis_prop_reordered > 1] - 2.0
  cis_prop_reordered[cis_prop_reordered <= -1] <- cis_prop_reordered[cis_prop_reordered <= -1] + 2.0
  
  df$R <- R
  df$theta_scaled <- theta_scaled
  df$cis_prop_reordered <- cis_prop_reordered
  
  if (plot) {
      p1 <- ggplot2::ggplot(df, ggplot2::aes(x = Parlog2FC, y = Hyblog2FC)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.7, color = df$colors) +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "orangered") +
        ggplot2::geom_abline(intercept = 0, slope = .5, color = "skyblue") +
        ggplot2::geom_abline(intercept = 0, slope = -5, color = "forestgreen") +
        ggplot2::geom_vline(xintercept = 0) +
        ggplot2::geom_hline(yintercept = 0, color = "darkblue") +
        ggplot2::theme_bw() +
#         ggplot2::ggtitle("Parental vs Hybrid Fold Change") +
        ggplot2::ylab(expression(R[H]~"(log"[2]*" hybrid fold change)")) + 
        ggplot2::xlab(expression(R[P]~"(log"[2]*" parental fold change)")) +
        ggplot2::scale_color_discrete(name = "Assignments")

      p2 <- ggplot2::ggplot(df, ggplot2::aes(x = delta, y = H)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.7, color = df$colors) +
        ggplot2::geom_vline(xintercept = 0, color = "orangered") +
        ggplot2::geom_hline(yintercept = 0, color = "darkblue") +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "skyblue") +
        ggplot2::geom_abline(intercept = 0, slope = -1, color = "forestgreen") +
        ggplot2::theme_bw() + 
        ggplot2::xlab(expression(R[P] - R[H])) + 
        ggplot2::ylab(expression(R[H])) +
        ggplot2::scale_color_discrete(name = "Assignments")
#         ggplot2::theme(legend.position = "none")
#         ggplot2::ggtitle("R vs Hybrid Fold Change")
      
      x_ <- ggplot_build(p2)$layout$panel_params[[1]]$x.range
      p2 <- p2 + ggplot2::plot( color = "skyblue") +
    
      custom_breaks <- c(-1.0, -0.5, 0.0, 0.5, 1.0)
      custom_labels <- c(0.05, 0.0, 0.5, 1.0, 0.5)
      
      p3 <- ggplot2::ggplot(df, ggplot2::aes(x = cis_prop_reordered, y = P)) +
        ggplot2::geom_point(size = 0.8, alpha = 0.7, color = df$colors) +
        ggplot2::geom_vline(xintercept = c(-1, -0.5, 0, 0.5, 1), color = c("forestgreen", "darkblue", "skyblue", "orangered", "forestgreen")) +
        ggplot2::geom_hline(yintercept = 0, color = "black") +
        ggplot2::theme_bw() + 
        ggplot2::xlab("Proportion cis") + 
        ggplot2::ylab(expression(R[P])) +
        scale_x_continuous(breaks = custom_breaks, labels = custom_labels) +
        ggplot2::scale_color_discrete(name = "Assignments")
#         ggplot2::theme(legend.position = "none")
#         ggplot2::ggtitle("Polar Representation of Cis-Trans Balance")
      
      if (!is.null(combo)) {
          title <- grid::textGrob(combo,gp = grid::gpar(fontsize = 16))} else{
          title <- grid::textGrob("",gp = grid::gpar(fontsize = 16))
      }
      
      plot_panel <- do.call(gridExtra::arrangeGrob, c(list(p1, p2, p3), nrow = 3))

      p <- gridExtra::grid.arrange(
          title,
          plot_panel,
          nrow = 2,
          heights = c(0.08, 1)
        )
      
    } else { p<- NULL }
                       
  
  return(list(df = df, plot = p))
}
                            

#' Plot histogram of reegulatory assignments.
#'
#' Creates a histogram-style bar plot showing the number of genes assigned to each regulatory category:
#' \code{cis}, \code{trans}, \code{cisxtrans}, \code{cis+trans}, and \code{conserved}.
#'
#' @param df A data frame containing a column named \code{reg_assignment} with one of the five regulatory assignment values.
#' @param title A string to use as the plot title. Defaults to \code{"Regulatory Categories Histogram"}.
#'
#' @return A \code{ggplot} object representing the histogram.
#'
#' @import ggplot2
#' @export                            
plotRegulatoryHistogram <- function(df, title = "Regulatory Categories Histogram") {
    
  df <- as.data.frame(df)
  print(unique(df$reg_assignment))
    
  # define the order of categories and corresponding colors
  categories <- c("cis", "trans", "cisxtrans", "cis+trans", "conserved")
  colors <- c("orangered", "darkblue", "forestgreen", "skyblue", "gray")
  color_map <- setNames(colors, categories)
    
  # coerce the reg_assignment column to a factor with desired order
  df$reg_assignment <- factor(df$reg_assignment, levels = categories)
    
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = reg_assignment, fill =  reg_assignment)) +
    ggplot2::geom_bar() +
    ggplot2::geom_text(stat = "count", ggplot2::aes(label = ..count..), vjust = -0.3, size = 3) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Regulatory assignment") +
    ggplot2::ylab("Number of genes") +
    ggplot2::ggtitle(title) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())

  return(p)
}
