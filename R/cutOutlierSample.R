#' Cut Outlier Sample
#'
#' Sets the cutoff value for sample clustering to exclude outliers.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param height Numeric value indicating the height cutoff for excluding outliers.
#'
#' @return A 'lacen' S3 object with updated `keepSamples` and filtered `datExpr` and `datTraits`.
#' @export
#' @import WGCNA
cutOutlierSample <- function(lacenObject,
                             height = FALSE) {
  UseMethod("cutOutlierSample")
}

#' @export
cutOutlierSample.lacen <- function(lacenObject,
                                   height = FALSE) {
  # Extract expression data and traits
  expr_data <- lacenObject$datExpr
  traits <- lacenObject$datTraits$Trait
  
  # Determine samples to keep
  if (isFALSE(height)) {
    keep_samples <- rep(TRUE, length(traits))  # Keep all samples
  } else {
    # Perform hierarchical clustering
    sample_tree <- fastcluster::hclust(stats::dist(expr_data), method = "average")
    # Cut the tree at the specified height
    clusters <- WGCNA::cutreeStatic(sample_tree, cutHeight = height, minSize = 10)
    # Display cluster sizes
    print(table(clusters))
    # Determine samples to keep (cluster 1)
    keep_samples <- (clusters == 1)
  }
  
  # Update the lacenObject with the samples to keep
  lacenObject$keepSamples <- keep_samples
  lacenObject$datExpr <- expr_data[keep_samples, ]
  lacenObject$datTraits <- traits[keep_samples]
  
  return(lacenObject)
}