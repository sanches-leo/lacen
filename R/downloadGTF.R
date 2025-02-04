#' Download GTF File
#'
#' Downloads a GTF annotation file from a given URL and returns it as a dataframe.
#'
#' @param link URL of the GTF annotation file.
#'
#' @return A dataframe with "gene_id" and "gene_name" columns.
#' @export
#' @import rtracklayer
downloadGTF <- function(link) {
  # Set a longer timeout for large downloads
  options(timeout = 600)
  
  # Import GTF file
  gtf_gr <- rtracklayer::import(link, format = "gtf")
  gtf_df <- as.data.frame(gtf_gr)
  
  # Extract relevant columns
  annotation_df <- gtf_df[, c("gene_id", "gene_name")]
  # Remove version numbers from gene_id and gene_name if present
  annotation_df$gene_id <- sub("\\..*$", "", annotation_df$gene_id)
  annotation_df$gene_name <- sub("\\..*$", "", annotation_df$gene_name)
  
  # Remove duplicate gene_ids
  annotation_df <- annotation_df[!duplicated(annotation_df$gene_id), ]
  
  return(annotation_df)
}
