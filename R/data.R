#' BRCA Gene Expression Data (Tumor vs Normal)
#'
#' Paired tumor-adjacent tissue RNA-seq data from TCGA BRCA project.
#' Processed using TCGAbiolinks with STAR counts workflow.
#'
#' @format A data frame with 60,616 genes (rows) and 226 samples (columns):
#' \describe{
#'   \item{rownames}{ENSEMBL gene IDs (e.g., "ENSG00000136158")}
#'   \item{columns}{TCGA sample IDs (e.g., "TCGA-AC-A2FM-11B-32R-A19W-07")}
#'   \item{values}{Raw counts from RNA-seq}
#' }
#' @source Processed from TCGA-BRCA data
"raw_expression"

#' Sample Phenotype Data
#'
#' Metadata for samples in datExpr indicating tumor/normal status.
#'
#' @format A data frame with 226 rows and 2 columns:
#' \describe{
#'   \item{Sample}{TCGA sample ID (matches datExpr columns)}
#'   \item{Trait}{Integer: 1 = Primary Tumor, 2 = Solid Tissue Normal}
#' }
#' @source Derived from TCGA metadata during datExpr processing
"traits"

#' Differential Expression Results
#'
#' Tumor vs normal differential expression analysis results using limma-voom.
#'
#' @format A data frame with 45,384 rows and 3 columns:
#' \describe{
#'   \item{gene_id}{ENSEMBL gene ID (e.g., "ENSG00000136158")}
#'   \item{log2FC}{Log2 fold change (Tumor vs Normal)}
#'   \item{pval}{Adjusted p-value from limma-voom analysis}
#' }
#' @source Generated using `TCGAbiolinks::TCGAanalyze_DEA`
"expression_DGE"

#' GENCODE v36 Gene Annotations
#'
#' Gene ID to gene name mappings from GENCODE release 36 (GRCh38.p13).
#'
#' @format A data frame with 67,016 rows and 2 columns:
#' \describe{
#'   \item{gene_id}{ENSEMBL gene ID (e.g., "ENSG00000223972")}
#'   \item{gene_name}{Official gene symbol (e.g., "DDX11L1")}
#' }
#' @source Subset from gencode.v36.annotation.gtf.gz
#' @references \url{https://www.gencodegenes.org/human/release_36.html}
"annotation_data"

#' Long Non-coding RNA Annotations
#'
#' Subset of GENCODE v36 annotations containing only lncRNA genes.
#'
#' @format A data frame with 17,946 rows and 2 columns:
#' \describe{
#'   \item{gene_id}{ENSEMBL gene ID (e.g., "ENSG00000243485")}
#'   \item{gene_name}{lncRNA symbol (e.g., "MIR1302-2HG")}
#' }
#' @source Filtered from gencode.v36.long_noncoding_RNAs.gtf.gz
"nc_annotation"
