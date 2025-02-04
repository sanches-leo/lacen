#' Select Soft Threshold
#'
#' Sets the chosen soft threshold power in the lacen object.
#'
#' @param lacenObject An object of class "lacen" created by [initLacen()].
#' @param indicePower Numeric value indicating the selected soft threshold power.
#'
#' @return A 'lacen' S3 object with updated `indicePower`.
#' @export
selectSoftThreshold <- function(lacenObject,
                                indicePower) {
  UseMethod("selectSoftThreshold")
}

#' @export
selectSoftThreshold.lacen <- function(lacenObject,
                                      indicePower) {
  # Validate that the provided power is numeric
  if (!is.numeric(indicePower) || length(indicePower) != 1) {
    stop("'indicePower' must be a single numeric value.")
  }
  
  # Update the lacenObject with the chosen soft threshold power
  lacenObject$indicePower <- indicePower
  return(lacenObject)
}