% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rf_c50.R
\name{summary.stratified_rf}
\alias{summary.stratified_rf}
\title{Summary statistics from a model}
\usage{
\method{summary}{stratified_rf}(object, ...)
}
\arguments{
\item{object}{A stratified_rf model.}

\item{...}{other options (not currently used)}
}
\description{
Calculates error statistics for out-of-bag samples from a stratified_rf model.
}
\details{
Predictions for a class are made by averaging class probabilities across trees rather than by a majority vote. All trees are weighted equally.
}
\examples{
data(iris)
groups <- list(c("Sepal.Length","Sepal.Width"),c("Petal.Length","Petal.Width"))
mtry <- c(1,1)
m <- stratified_rf(iris,"Species",groups,mtry,ntrees=2,multicore=FALSE)
summary(m)
}
