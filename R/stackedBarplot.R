#' Module-Correlation Stacked Barplot
#'
#' Creates a barplot illustrating eigengene modules, their gene types, and trait-module correlations.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param filename Filename to save the stacked barplot as PNG.
#' @param plot Logical. If TRUE, displays the plot in the R session.
#'
#' @return None (plots and saves the stacked barplot).
#' @export
#' @import WGCNA ggplot2
stackedBarplot <- function(lacenObject,
                           filename = "stackedplot_desk.png",
                           plot = TRUE) {
  UseMethod("stackedBarplot")
}

#' @export
stackedBarplot.lacen <- function(lacenObject,
                                 filename = "stackedplot_desk.png",
                                 plot = TRUE) {
  # Extract necessary data
  summary_df <- lacenObject$summdf
  expr_data <- lacenObject$datExpr
  traits <- lacenObject$datTraits$Trait

  # Prepare genes selected based on bootstrap
  genes_selected <- summary_df$gene_id[summary_df$cutBootstrap]

  # Filter expression data for selected genes
  expr_selected <- expr_data[, colnames(expr_data) %in% genes_selected]

  # Extract non-coding annotations
  nc_data <- summary_df$is_nc

  # Initialize data for stacked barplot
  barplot_data <- data.frame(Module = character(),
                             nGenes = integer(),
                             geneType = character(),
                             correlation = numeric(),
                             stringsAsFactors = FALSE)

  # Iterate over each module to compile data
  unique_modules <- sort(unique(summary_df$module))
  for (mod in unique_modules) {
    if (mod != 0) {
      # Subset genes in the current module
      genes_in_mod <- summary_df$module == mod
      module_genes <- summary_df$gene_id[genes_in_mod]

      # Count non-coding and protein-coding genes
      n_lnc <- sum(module_genes %in% summary_df$gene_id[summary_df$is_nc])
      n_pc <- length(module_genes) - n_lnc

      # Calculate module eigengenes and correlation with traits
      module_eigengenes <- WGCNA::moduleEigengenes(expr_selected, as.numeric(genes_in_mod))$eigengenes
      ordered_MEs <- WGCNA::orderMEs(module_eigengenes)
      module_trait_cor <- WGCNA::cor(ordered_MEs, traits, use = "p")
      module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nrow(expr_selected))
      correlation <- module_trait_cor[[2]]

      # Assign correlation to zero if p-value is not significant
      if (module_trait_pval[2, ] > 0.05) {
        correlation <- 0
      }

      # Append data to barplot_data
      barplot_data <- rbind(barplot_data, data.frame(
        Module = paste0("Module ", mod),
        nGenes = n_lnc,
        geneType = "lnc",
        correlation = correlation,
        stringsAsFactors = FALSE
      ))
      barplot_data <- rbind(barplot_data, data.frame(
        Module = paste0("Module ", mod),
        nGenes = n_pc,
        geneType = "pc",
        correlation = correlation,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Order modules based on correlation
  barplot_data <- barplot_data[order(barplot_data$correlation), ]
  barplot_data$Module <- factor(barplot_data$Module, levels = unique(barplot_data$Module))

  # Create separate columns for positive and negative correlations
  barplot_data$correlationPos <- ifelse(barplot_data$correlation > 0, barplot_data$correlation * max(barplot_data$nGenes), NA)
  barplot_data$correlationNeg <- ifelse(barplot_data$correlation < 0, abs(barplot_data$correlation) * max(barplot_data$nGenes), NA)

  # Define colors for modules
  module_colors <- c("grey", WGCNA::labels2colors(sort(unique_modules)))

  # Create the stacked barplot using ggplot2
  ggplot_obj <- ggplot2::ggplot(barplot_data, ggplot2::aes(x = "Module")) +
    ggplot2::geom_bar(aes(y = "nGenes", fill = "geneType"), stat = "identity") +
    ggplot2::scale_fill_manual(values = c("lnc" = "grey", "pc" = "blue")) +
    ggplot2::xlab("") +
    ggplot2::ylab("Number of Protein-coding (blue) and Non-coding Genes (grey)") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Gene Type")) +
    ggplot2::geom_point(aes(y = "correlationPos"), shape = 19, size = 3) +
    ggplot2::geom_point(aes(y = "correlationNeg"), shape = 1, size = 3) +
    ggplot2::scale_y_continuous(
      sec.axis = ggplot2::sec_axis(~ . / max(barplot_data$nGenes), name = "Trait-Module Correlation")
    ) +
    ggplot2::ggtitle("Transcript Type per Module and Trait-Module Correlation") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "right",
      text = ggplot2::element_text(size = 20),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title = ggplot2::element_text(size = 16)
    )

  # Save the plot as PNG
  grDevices::png(filename = filename, width = 1000, height = 600)
  print(ggplot_obj)
  grDevices::dev.off()

  # Optionally display the plot in the R session
  if (isTRUE(plot)) {
    print(ggplot_obj)
  }
}
