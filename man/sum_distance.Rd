% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sum_distance.R
\name{sum_distance}
\alias{sum_distance}
\title{Sum the distance of all branches in a dendrogram}
\usage{
sum_distance(
  d,
  sum = 0,
  mapping = c("exponential", "linear"),
  scale = 0,
  max_height = 0,
  ...
)
}
\arguments{
\item{d}{dendrogram class object}

\item{sum}{the sum of distance}

\item{mapping}{in case of quantitative mode, either "linear" or "exponential" mapping}

\item{...}{ignored}
}
\value{
the sum of distances
}
\description{
This function takes a dendrogram class object as an input, and adds up all the distances of branches. 
This function is called recursively to adds up the sum. In case of exponential mapping for the quantitative
mode, the sum is in the exponential scale
}
\keyword{internal}
