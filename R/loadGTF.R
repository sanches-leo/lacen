#' Load GTF File
#'
#' Loads a local GTF annotation file and returns it as a dataframe.
#'
#' @param file Path to the local GTF file.
#'
#' @return A dataframe with "gene_id" and "gene_name" columns.
#' @export
#' @importFrom rtracklayer import
loadGTF <- function(file) {
  # Import GTF file
  gtf_gr <- rtracklayer::import(file, format = "gtf")
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
