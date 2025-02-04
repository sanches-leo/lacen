#' Plot Soft Threshold
#'
#' Generates a plot to help select the soft-thresholding power for network construction.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param filename Filename to save the soft threshold plot as PNG.
#' @param maxBlockSize Maximum number of genes processed at the same time.
#' @param plot Logical. If TRUE, displays the plot in the R session.
#'
#' @return None (plots and saves the soft threshold plot).
#' @export
#' @import WGCNA
plotSoftThreshold <- function(lacenObject,
                              filename = "indicePower.png",
                              maxBlockSize = 10000,
                              plot = TRUE) {
  UseMethod("plotSoftThreshold")
}

#' @export
plotSoftThreshold.lacen <- function(lacenObject,
                                    filename = "indicePower.png",
                                    maxBlockSize = 10000,
                                    plot = TRUE) {
  # Extract expression data
  expr_data <- lacenObject$datExpr
  
  # Define a range of soft-thresholding powers to evaluate
  powers <- 1:20
  
  # Perform soft-thresholding analysis
  soft_threshold <- WGCNA::pickSoftThreshold(expr_data, powerVector = powers, verbose = 3, blockSize = maxBlockSize)
  
  # Open PNG device to save the plot
  grDevices::png(filename = filename, units = "in", width = 9, height = 5, res = 96)
  graphics::par(mfrow = c(1, 2))  # Set plot layout
  
  # Plot Scale-Free Topology Fit Index
  plot(soft_threshold$fitIndices[, 1], -sign(soft_threshold$fitIndices[, 3]) * soft_threshold$fitIndices[, 2],
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R\u00B2",
       type = "n", main = "Scale Independence")
  graphics::text(soft_threshold$fitIndices[, 1], 
                -sign(soft_threshold$fitIndices[, 3]) * soft_threshold$fitIndices[, 2],
                labels = powers, cex = 0.9, col = "red")
  graphics::abline(h = 0.90, col = "red")  # Reference line for R2
  
  # Plot Mean Connectivity
  plot(soft_threshold$fitIndices[, 1], soft_threshold$fitIndices[, 5], log = 'y',
       xlab = "Soft Threshold (power)", 
       ylab = "Mean Connectivity",
       type = "n", main = "Mean Connectivity")
  graphics::text(soft_threshold$fitIndices[, 1], soft_threshold$fitIndices[, 5],
                labels = powers, cex = 0.9, col = "red")
  
  # Close the PNG device
  grDevices::dev.off()
  
  # If plotting is requested, display the plot in the R session
  if (isTRUE(plot)) {
    graphics::par(mfrow = c(1, 2))
    
    # Scale Independence Plot
    plot(soft_threshold$fitIndices[, 1], 
         -sign(soft_threshold$fitIndices[, 3]) * soft_threshold$fitIndices[, 2],
         xlab = "Soft Threshold (power)", 
         ylab = "Scale Free Topology Model Fit, signed R\u00B2",
         type = "n", main = "Scale Independence")
    graphics::text(soft_threshold$fitIndices[, 1], 
                  -sign(soft_threshold$fitIndices[, 3]) * soft_threshold$fitIndices[, 2],
                  labels = powers, cex = 0.9, col = "red")
    graphics::abline(h = 0.90, col = "red")
    
    # Mean Connectivity Plot
    plot(soft_threshold$fitIndices[, 1], 
         soft_threshold$fitIndices[, 5], log = 'y',
         xlab = "Soft Threshold (power)", 
         ylab = "Mean Connectivity",
         type = "n", main = "Mean Connectivity")
    graphics::text(soft_threshold$fitIndices[, 1], 
                  soft_threshold$fitIndices[, 5],
                  labels = powers, cex = 0.9, col = "red")
  }
}