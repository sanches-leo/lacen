% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selectOutlierSample.R
\name{selectOutlierSample}
\alias{selectOutlierSample}
\title{Select Outlier Sample}
\usage{
selectOutlierSample(
  lacenObject,
  plot = TRUE,
  filename = "clusterTree.png",
  height = FALSE
)
}
\arguments{
\item{lacenObject}{An object of class "lacen" created by \code{\link[=initLacen]{initLacen()}}.}

\item{plot}{Logical. If TRUE, plots the cluster tree in the R session.}

\item{filename}{Filename to save the cluster tree plot as PNG.}

\item{height}{Numeric value indicating the height threshold to separate outliers.}
}
\value{
None (plots and saves the cluster tree).
}
\description{
Plots the sample cluster tree to help identify and exclude outlier samples.
Returns a vector of samples to keep based on the provided height threshold.
}
