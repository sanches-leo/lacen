#' Filter and Transform
#'
#' Removes genes with low variation, filters expression data by MAD or DEG, and applies voom correction.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param pThreshold p-value threshold for differential expression.
#' @param fcThreshold Fold change threshold for differential expression.
#' @param topVarGenes Number of top variable genes to retain if filtering by variance.
#' @param filterMethod Method to filter genes: "DEG" (default) or "var".
#'
#' @return A 'lacen' S3 object with updated `datExpr`.
#' @export
filterTransform <- function(lacenObject,
                            pThreshold = 0.01,
                            fcThreshold = 1,
                            topVarGenes = 5000,
                            filterMethod = "DEG") {
  UseMethod("filterTransform")
}

#' @export
filterTransform.lacen <- function(lacenObject,
                                  pThreshold = 0.01,
                                  fcThreshold = 1,
                                  topVarGenes = 5000,
                                  filterMethod = "DEG") {
  # Extract necessary data from lacenObject
  expression_data <- lacenObject$datExpression
  count_data <- lacenObject$datCounts
  
  # Internal function to filter DEGs based on p-value and fold change thresholds
  filter_DEGs <- function(expr_data, p_thresh, fc_thresh) {
    significant <- expr_data$pval < p_thresh & !is.na(expr_data$pval)
    high_fc <- abs(expr_data$log2FC) > fc_thresh & !is.na(expr_data$log2FC)
    deg <- expr_data[significant & high_fc, ]
    return(deg)
  }
  
  # Internal function to select top variable genes by MAD
  select_by_variance <- function(count_data, top_var) {
    # Calculate Median Absolute Deviation (MAD) for each gene
    gene_mads <- apply(as.matrix(count_data[, -1]), 1, stats::mad)
    # Select top variable genes
    top_genes <- names(sort(gene_mads, decreasing = TRUE))[1:top_var]
    filtered_data <- count_data[rownames(count_data) %in% top_genes, ]
    return(filtered_data)
  }
  
  # Internal function to filter counts by DEG standard deviation
  filter_counts_by_deg_sd <- function(count_data, deg_data) {
    # Calculate standard deviation for each gene
    gene_sd <- apply(count_data, 1, stats::sd, na.rm = TRUE)
    count_data$SD <- gene_sd
    # Merge DEG data with standard deviations
    deg_sd <- merge(deg_data[, c("gene_id", "log2FC")], 
                    data.frame(gene_id = rownames(count_data), SD = gene_sd),
                    by = "gene_id")
    # Determine the 1st percentile SD threshold
    sd_threshold <- quantile(deg_sd$SD, probs = 0.01, na.rm = TRUE)
    # Filter genes with SD above the threshold
    filtered_counts <- count_data[count_data$SD > sd_threshold, -ncol(count_data)]
    return(filtered_counts)
  }
  
  # Ensure at least 1/4 of samples have expression > 1
  required_expressions <- ceiling((ncol(count_data) - 1) / 4)
  keep_genes <- apply(count_data[, -1], 1, function(x) sum(x > 1) > required_expressions)
  count_data <- count_data[keep_genes, ]
  
  if (filterMethod == "DEG") {
    # Filter DEGs
    degs <- filter_DEGs(expression_data[expression_data$gene_id %in% rownames(count_data), ], 
                       pThreshold, fcThreshold)
    # Further filter counts based on DEG SD
    filtered_expr <- filter_counts_by_deg_sd(count_data, degs)
  } else if (filterMethod == "var") {
    # Select top variable genes by MAD
    filtered_expr <- select_by_variance(count_data, topVarGenes)
  } else {
    stop("Invalid filterMethod. Choose 'DEG' or 'var'.")
  }
  
  # Apply voom transformation using limma
  voom_transformed <- limma::voom(t(filtered_expr))$E
  voom_matrix <- data.matrix(voom_transformed, rownames.force = TRUE)
  
  # Update the lacenObject with the transformed expression data
  lacenObject$datExpr <- voom_matrix
  return(lacenObject)
}
