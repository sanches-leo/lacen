#' Heatmap of Top Connectivity
#'
#' Generates a heatmap of a specified module or submodule, showing interconnectedness based on the TOM.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param module Numeric value indicating the module to analyze.
#' @param submodule Optional. Numeric value indicating the submodule to analyze. Set to FALSE to ignore.
#' @param hmDimensions Numeric value indicating the number of genes to display in the heatmap.
#' @param filename Filename to save the heatmap as PNG.
#' @param removeNonDEG Logical. If TRUE, displays only differentially expressed genes.
#' @param outTSV Logical or character string. If FALSE, does not save. If a path is provided, saves the TSV with heatmap values.
#' @param ... Additional parameters for internal functions.
#'
#' @return None (plots and saves the heatmap).
#' @export
#' @import WGCNA igraph ggraph grid ggplot2
heatmapTopConnectivity <- function(lacenObject,
                                   module,
                                   submodule = FALSE,
                                   hmDimensions = FALSE,
                                   filename = FALSE,
                                   removeNonDEG = FALSE,
                                   outTSV = FALSE,
                                   ...) {
  UseMethod("heatmapTopConnectivity")
}

#' @export
heatmapTopConnectivity.lacen <- function(lacenObject,
                                          module,
                                          submodule = FALSE,
                                          hmDimensions = FALSE,
                                          filename = FALSE,
                                          removeNonDEG = FALSE,
                                          outTSV = FALSE,
                                          ...) {
  # Extract necessary data
  expr_data <- lacenObject$datExpr
  traits <- lacenObject$datTraits$Trait
  tom_matrix <- lacenObject$TOM
  summary_df <- lacenObject$summdf
  enrichment_list <- lacenObject$rrvgolist
  
  # Validate module parameter
  if (!is.numeric(module) || module <= 0) {
    stop("Invalid 'module' parameter. It must be a positive integer.")
  }
  
  if (!(module %in% unique(summary_df$module))) {
    available_modules <- sort(unique(summary_df$module))
    available_modules <- available_modules[available_modules > 0]
    stop("Invalid 'module'. Available modules are: ", paste(available_modules, collapse = ", "))
  }
  
  # Handle submodule parameter
  if (!isFALSE(submodule)) {
    if (!is.numeric(submodule) || submodule <= 0) {
      stop("Invalid 'submodule' parameter. It must be a positive integer or FALSE.")
    }
    possible_submodules <- unique(enrichment_list[[as.character(module)]]$cluster)
    if (!(submodule %in% possible_submodules)) {
      if (length(possible_submodules) == 0) {
        stop("This module has no enriched submodules.")
      } else {
        stop("Invalid 'submodule'. Available submodules are: ", paste(possible_submodules, collapse = ", "))
      }
    }
  }
  
  # Determine filename if not provided
  if (isFALSE(filename)) {
    if (isFALSE(submodule)) {
      filename <- paste0("7_heatmap_", module, ".png")
    } else {
      filename <- paste0("7_heatmap_", module, "_", submodule, ".png")
    }
  }
  
  # Identify genes in the specified module and submodule
  if (isFALSE(submodule)) {
    genes_in_module <- summary_df$module == module & !summary_df$is_nc
  } else {
    genes_in_module <- summary_df$module == module & summary_df[[as.character(submodule)]] & !summary_df$is_nc
  }
  module_genes <- summary_df$gene_id[genes_in_module]
  
  # Calculate median connectivity within the module
  module_connectivity <- summary_df$kWithin[summary_df$module == module]
  median_connectivity <- stats::median(module_connectivity, na.rm = TRUE)
  
  # Identify lncRNAs with connectivity above the median
  lnc_highly_connected <- summary_df$gene_id[summary_df$module == module & 
                                              summary_df$kWithin > median_connectivity & 
                                              summary_df$is_nc]
  
  # Subset TOM for these lncRNAs and protein-coding genes
  tom_subset <- tom_matrix[rownames(tom_matrix) %in% lnc_highly_connected,
                           colnames(tom_matrix) %in% module_genes]
  
  # Check if there are lncRNAs with high connectivity
  if (length(lnc_highly_connected) == 0) {
    stop("No lncRNAs in the module have connectivity above the median.")
  }
  
  # Ensure matrix dimensions are appropriate
  if (!(nrow(tom_subset) > 0 && ncol(tom_subset) > 0)) {
    stop("No lncRNAs with high connectivity found in the module.")
  }
  
  # Calculate module eigengenes and correlation with traits
  module_eigengenes <- WGCNA::moduleEigengenes(expr_data, as.numeric(genes_in_module))$eigengenes
  ordered_MEs <- WGCNA::orderMEs(module_eigengenes)
  module_trait_cor <- WGCNA::cor(ordered_MEs, traits, use = "p")
  module_trait_pval <- WGCNA::corPvalueStudent(module_trait_cor, nrow(expr_data))
  
  # Extract connectivity and DEG information for protein-coding genes
  pc_connectivity <- summary_df$kWithin[summary_df$gene_id %in% colnames(tom_subset)]
  names(pc_connectivity) <- summary_df$gene_id[summary_df$gene_id %in% colnames(tom_subset)]
  pc_connectivity <- sort(pc_connectivity, decreasing = TRUE)
  
  pc_DEG <- summary_df$gene_id[summary_df$gene_id %in% colnames(tom_subset)]
  pc_DEG <- pc_DEG[match(names(pc_connectivity), pc_DEG)]
  
  # Extract connectivity and DEG information for lncRNAs
  lnc_connectivity <- summary_df$kWithin[summary_df$gene_id %in% rownames(tom_subset)]
  names(lnc_connectivity) <- summary_df$gene_id[summary_df$gene_id %in% rownames(tom_subset)]
  lnc_connectivity <- sort(lnc_connectivity, decreasing = TRUE)
  
  lnc_DEG <- summary_df$gene_id[summary_df$gene_id %in% rownames(tom_subset)]
  lnc_DEG <- lnc_DEG[match(names(lnc_connectivity), lnc_DEG)]
  
  # Subset the TOM matrix
  tom_subset <- tom_subset[match(names(lnc_connectivity), rownames(tom_subset)),
                           match(names(pc_connectivity), colnames(tom_subset))]
  
  # Define colors based on connectivity and DEG status
  define_colors <- function(connectivity, DEG, mod) {
    colors <- rep("black", length(connectivity))
    # Highlight DEGs
    colors[DEG & abs(summary_df$log2FC) >= 1 & summary_df$pval <= 0.05] <- "#FF0000"  # Upregulated
    colors[DEG & summary_df$log2FC <= -1 & summary_df$pval <= 0.05] <- "#0000FF"  # Downregulated
    # Color based on connectivity
    max_conn <- max(connectivity, na.rm = TRUE)
    min_conn <- min(connectivity, na.rm = TRUE)
    pal <- grDevices::colorRampPalette(c("#c8fac8", "#003200"))(128)
    conn_scaled <- pmin(pmax(connectivity, min_conn), max_conn)
    conn_colors <- pal[findInterval(conn_scaled, seq(min_conn, max_conn, length.out = length(pal)), all.inside = TRUE)]
    colors[!DEG] <- conn_colors[!DEG]
    return(colors)
  }
  
  pc_colors <- define_colors(pc_connectivity, !is.na(summary_df$pval[match(names(pc_connectivity), summary_df$gene_id)]), module)
  lnc_colors <- define_colors(lnc_connectivity, !is.na(summary_df$pval[match(names(lnc_connectivity), summary_df$gene_id)]), module)
  
  # Prepare heatmap side colors
  colside <- cbind(DEG = pc_colors, Connectivity = pc_colors)
  rowside <- cbind(DEG = lnc_colors, Connectivity = lnc_colors)
  
  # Subset TOM based on hmDimensions if specified
  if (!isFALSE(hmDimensions)) {
    hm_dimensions <- as.integer(hmDimensions)
    if (is.na(hm_dimensions) || hm_dimensions <= 0) {
      stop("'hmDimensions' must be a positive integer.")
    }
    hm_dimensions_lnc <- min(hm_dimensions, nrow(tom_subset))
    hm_dimensions_pc <- min(hm_dimensions, ncol(tom_subset))
    tom_subset <- tom_subset[1:hm_dimensions_lnc, 1:hm_dimensions_pc]
    colside <- colside[1:hm_dimensions_lnc, ]
    rowside <- rowside[1:hm_dimensions_pc, ]
  }
  
  # Remove non-DEGs if requested
  if (isTRUE(removeNonDEG)) {
    tom_subset <- tom_subset[!(rownames(tom_subset) %in% summary_df$gene_id[!summary_df$is_nc & !is.na(summary_df$pval)]), 
                             !(colnames(tom_subset) %in% summary_df$gene_id[!summary_df$is_nc & !is.na(summary_df$pval)])]
    colside <- colside[!(rownames(colside) %in% summary_df$gene_id[!summary_df$is_nc & !is.na(summary_df$pval)]), ]
    rowside <- rowside[!(colnames(rowside) %in% summary_df$gene_id[!summary_df$is_nc & !is.na(summary_df$pval)]), ]
  }
  
  # Replace gene IDs with gene names for labels
  rownames(tom_subset) <- summary_df$gene_name[match(rownames(tom_subset), summary_df$gene_id)]
  colnames(tom_subset) <- summary_df$gene_name[match(colnames(tom_subset), summary_df$gene_id)]
  
  # Define color palettes for connectivity
  connectivity_palette <- rev(grDevices::grey.colors(100))
  
  # Open PNG device for heatmap
  grDevices::png(filename = filename, width = 1800, height = 1800, res = 150)
  
  # Generate the heatmap using custom heatmap.3 function or another suitable function
  # Here, using WGCNA's built-in heatmap functionality for simplicity
  WGCNA::plotTOM(tom_subset, main = paste("Top Connectivity Heatmap for Module", module),
                TOM.type = "unsigned", threshold = 0, palette = "RdBu",
                mar = c(10,10,3,3), dendro = "none", plotLegend = TRUE)
  
  # Close the PNG device
  grDevices::dev.off()
  
  # Optionally display the heatmap in the R session
  if (isTRUE(plot)) {
    WGCNA::plotTOM(tom_subset, main = paste("Top Connectivity Heatmap for Module", module),
                  TOM.type = "unsigned", threshold = 0, palette = "RdBu",
                  mar = c(10,10,3,3), dendro = "none", plotLegend = TRUE)
  }
  
  # Save the heatmap data as TSV if requested
  if (!isFALSE(outTSV)) {
    heatmap_values <- tom_subset
    # Optionally, transform the data or add gene annotations
    utils::write.table(heatmap_values, file = outTSV, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  }
}