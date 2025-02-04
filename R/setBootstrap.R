#' Set Bootstrap Cut-off
#'
#' Sets a cutoff value for bootstrap stability to filter out less stable genes.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param cutBootstrap Numeric value between 0 and 1 indicating the stability cutoff.
#'
#' @return A 'lacen' S3 object with updated `cutBootstrap`.
#' @export
setBootstrap <- function(lacenObject,
                         cutBootstrap) {
  UseMethod("setBootstrap")
}

#' @export
setBootstrap.lacen <- function(lacenObject,
                               cutBootstrap) {
  # Validate the cutoff value
  if (!is.numeric(cutBootstrap) || length(cutBootstrap) != 1 || 
      cutBootstrap < 0 || cutBootstrap > 1) {
    stop("'cutBootstrap' must be a single numeric value between 0 and 1.")
  }
  
  # Update the lacenObject with the bootstrap cutoff
  lacenObject$cutBootstrap <- cutBootstrap
  return(lacenObject)
}