% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flexwinn.R
\name{flexwinn}
\alias{flexwinn}
\title{flexwinn}
\usage{
flexwinn(data, graph = FALSE)
}
\arguments{
\item{data}{input dataset of metabolites as data frame.}

\item{graph}{TRUE or FALSE. If TRUE will display the graph of uncorrected
data with estimated change points and correction and the corrected signal
on a second plot.}
}
\value{
the corrected dataframe of metabolite(s).
}
\description{
A function for metabolite correction.
This function first estimate the change points
(jumps) and then correct each different segment
by detrending and normalizing when needed.
}
\examples{
met1<-rnorm(200)
flexwinn(as.data.frame(met1))
}
