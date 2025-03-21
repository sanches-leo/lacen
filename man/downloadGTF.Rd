% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/downloadGTF.R
\name{downloadGTF}
\alias{downloadGTF}
\title{Download GTF File}
\usage{
downloadGTF(link)
}
\arguments{
\item{link}{URL of the GTF annotation file.}
}
\value{
A dataframe with "gene_id" and "gene_name" columns.
}
\description{
Downloads a GTF annotation file from a given URL and returns it as a dataframe.
}
