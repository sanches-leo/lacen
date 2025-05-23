% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{raw_expression}
\alias{raw_expression}
\title{BRCA Gene Expression Data (Tumor vs Normal)}
\format{
A data frame with 30,000 genes (rows) and 60 samples (columns):
\describe{
\item{rownames}{ENSEMBL gene IDs (e.g., "ENSG00000136158")}
\item{columns}{TCGA sample IDs (e.g., "TCGA-AC-A2FM-11B-32R-A19W-07")}
\item{values}{Raw counts from RNA-seq}
}
}
\source{
Processed from TCGA-BRCA data
}
\usage{
raw_expression
}
\description{
Random subset of Paired tumor-adjacent tissue RNA-seq data from TCGA BRCA project.
Processed using TCGAbiolinks with STAR counts workflow.
}
\keyword{datasets}
