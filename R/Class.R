#' @name RFclust.SGE
#' @title RFclust.SGE
#' @docType package
#' @description  An S4 class to run unsupervides clustering based on random forrest predictors using the Sun Grid Engine.
#' @slot dat the data that should be clustered (data.frame with column == samples rows == observations)
#' @slot email the users email for the SGE report
#' @slot slices the number of slices the data should be analyzed in
#' @slot tmp.path the temp path for the clustering results
#' @slot SGE run using Sun Grid Engine (default=F)
#' @slot distRF the density data (internal use only)
#' @slot RFfiles an internal list of RF files (internal use only)
#' @slot name the name for this object that will be used to identify the data object in the spawned processes
#' @source \url{https://labs.genetics.ucla.edu/horvath/RFclustering/RFclustering.htm}
#' @exportClass RFclust.SGE
setClass( 
		Class='RFclust.SGE', 
		representation = representation ( 
				dat='data.frame',
				email='character',
				slices='numeric',
				tmp.path='character',
				SGE='logical',
				distRF='list',
				RFfiles='list',
				name='character'
		),
		prototype(tmp.path =NA_character_,  email = NA_character_, 
				slices =32, SGE=FALSE, distRF=list(), RFfiles=list()
		)
)
