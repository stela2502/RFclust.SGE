% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RFclust.SGE.R
\docType{methods}
\name{read.RF}
\alias{read.RF}
\alias{read.RF,RFclust.SGE-method}
\title{description of function read.RF}
\usage{
read.RF(x, name, max.wait = 20)
}
\arguments{
\item{max.wait}{maximum time to wait for the files to become accessable default= 20}

\item{files}{a list of files created by save.RF() default=c('')}
}
\description{
This function reads a set of outfiles and creates a summary Rf.data

object that can be processed using the RFdist() function.
}
