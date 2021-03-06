% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statepaintr.R
\name{decisionMatrix}
\alias{decisionMatrix}
\title{decisionMatrix}
\usage{
decisionMatrix(object)
}
\arguments{
\item{object}{of class decisionMatrix}
}
\value{
a matrix; The description of what combination of features are required to call a state.
}
\description{
decisionMatrix
}
\examples{
load(system.file("extdata", "poised.promoter.model.rda", package = "StatePaintR"))
poised.promoter.model
decisionMatrix(poised.promoter.model)
}
