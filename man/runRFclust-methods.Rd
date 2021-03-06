% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Interface.R
\docType{methods}
\name{runRFclust}
\alias{runRFclust}
\alias{runRFclust,RFclust.SGE-method}
\title{description of function runRFclust}
\usage{
runRFclust(x, ntree = 500, nforest = 500, name = "RFrun", force = FALSE)
}
\arguments{
\item{x}{the RFclust.SGE object}

\item{ntree}{the number of trees to grow}

\item{nforest}{the nuber of forests to create}

\item{name}{the name of the random forest clustering run (if you want to run multiple)}
}
\value{
a distRF object to be analyzed by pamNew
}
\description{
run the random forest calculations returning the density matrix

the clusters will be created for the columns of the data.frame
}
