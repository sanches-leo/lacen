#' Select Outlier Sample
#'
#' Plots the sample cluster tree to help identify and exclude outlier samples.
#' Returns a vector of samples to keep based on the provided height threshold.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param plot Logical. If TRUE, plots the cluster tree in the R session.
#' @param filename Filename to save the cluster tree plot as PNG.
#' @param height Numeric value indicating the height threshold to separate outliers.
#'
#' @return None (plots and saves the cluster tree).
#' @export
#' @import WGCNA
selectOutlierSample <- function(lacenObject,
                                plot = TRUE,
                                filename = "clusterTree.png",
                                height = FALSE) {
  UseMethod("selectOutlierSample")
}

#' @export
selectOutlierSample.lacen <- function(lacenObject,
                                      plot = TRUE,
                                      filename = "clusterTree.png",
                                      height = FALSE) {
  # Extract expression data and traits
  expr_data <- lacenObject$datExpr
  traits <- lacenObject$datTraits$Trait
  
  # Open PNG device if a filename is provided
  grDevices::png(filename = filename)
  
  # Plot cluster tree without a height threshold
  if (isFALSE(height)) {
    WGCNA::plotClusterTreeSamples(datExpr = expr_data, y = traits)
  } else {
    # Plot cluster tree with a height threshold to identify outliers
    WGCNA::plotClusterTreeSamples(datExpr = expr_data, y = traits, abHeight = height)
  }
  
  # Close the PNG device
  grDevices::dev.off()
  
  # If height is provided, classify samples based on the height cutoff
  if (!isFALSE(height)) {
    # Perform hierarchical clustering
    sample_tree <- fastcluster::hclust(stats::dist(expr_data), method = "average")
    # Cut the tree at the specified height
    clusters <- WGCNA::cutreeStatic(sample_tree, cutHeight = height, minSize = 10)
    # Display cluster sizes
    print(table(clusters))
    # Determine samples to keep (cluster 1)
    keep_samples <- (clusters == 1)
    # Update the lacenObject with the samples to keep
    lacenObject$keepSamples <- keep_samples
    # Filter expression data and traits based on kept samples
    lacenObject$datExpr <- expr_data[keep_samples, ]
    lacenObject$datTraits <- lacenObject$datTraits[lacenObject$datTraits$Sample %in% rownames(lacenObject$datExpr), ]
    
    # Re-plot the cluster tree with highlighted outliers if requested
    if (isTRUE(plot)) {
      WGCNA::plotClusterTreeSamples(datExpr = lacenObject$datExpr, y = lacenObject$datTraits$Trait, abHeight = height)
    }
  }
}