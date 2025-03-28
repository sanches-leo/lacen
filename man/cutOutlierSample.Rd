% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cutOutlierSample.R
\name{cutOutlierSample}
\alias{cutOutlierSample}
\title{Cut Outlier Sample}
\usage{
cutOutlierSample(lacenObject, height = FALSE)
}
\arguments{
\item{lacenObject}{An object of class "lacen" created by \code{\link[=initLacen]{initLacen()}}.}

\item{height}{Numeric value indicating the height cutoff for excluding outliers.}
}
\value{
A 'lacen' S3 object with updated \code{keepSamples} and filtered \code{datExpr} and \code{datTraits}.
}
\description{
Sets the cutoff value for sample clustering to exclude outliers.
}
