#' Lacen Constructor
#'
#' Create a new 'lacen' S3 object.
#'
#' @param datCounts Data counts dataframe with genes as rows and samples as columns.
#' The gene/transcript IDs should be in the row names and sample names/IDs in the column names.
#' @param datExpression Differential expression dataframe. Can be from Limma output or a dataframe with columns "ID", "log2FC", and "pvalue".
#' @param datTraits Conditions/Traits data. A two-column dataframe with sample IDs in "Sample" and condition codes in "Trait".
#' @param annotationData Annotation dataframe with "gene_id" and "gene_name" columns.
#' @param ncAnnotation Subset of "annotationData" containing only long non-coding RNAs.
#' @param datExpr Filtered and transformed expression data obtained via [filterTransform()].
#' @param keepSamples Boolean vector indicating samples to retain based on outlier removal.
#' @param height Numeric value to determine outlier samples from the cluster tree plot.
#' @param indicePower Soft threshold power used for network construction.
#' @param bootstrapStability Stability ratios from bootstrap iterations via [lacenBootstrap()].
#' @param cutBootstrap Threshold to discard low stability genes.
#' @param modGroups Module groups from bootstrap parameter set.
#' @param bootstrap Bootstrap parameters set by [lacenBootstrap()].
#' @param summdf Summarized results from [summarizeAndEnrichModules()].
#' @param rrvgolist Enriched and reduced terms from [summarizeAndEnrichModules()].
#' @param TOM Topological overlap matrix from [summarizeAndEnrichModules()].
#'
#' @return A 'lacen' S3 object containing all the provided parameters.
#' @export
initLacen <- function(datCounts = NULL,
                      datExpr = NULL,
                      datExpression = NULL,
                      datTraits = NULL,
                      annotationData = NULL,
                      ncAnnotation = NULL,
                      keepSamples = NULL,
                      height = NULL,
                      indicePower = NULL,
                      modGroups = NULL,
                      bootstrapStability = FALSE,
                      bootstrap = FALSE,
                      cutBootstrap = FALSE,
                      summdf = NULL,
                      rrvgolist = NULL,
                      TOM = NULL) {
  # Initialize the list with provided parameters
  lacen_list <- list(
    datCounts = datCounts,
    datExpr = datExpr,
    datExpression = datExpression,
    datTraits = datTraits,
    annotationData = annotationData,
    ncAnnotation = ncAnnotation,
    keepSamples = keepSamples,
    height = height,
    indicePower = indicePower,
    modGroups = modGroups,
    bootstrapStability = bootstrapStability,
    bootstrap = FALSE,
    cutBootstrap = cutBootstrap,
    summdf = summdf,
    rrvgolist = rrvgolist,
    TOM = TOM
  )
  
  # Assign 'lacen' class to the list
  structure(lacen_list, class = "lacen")
}