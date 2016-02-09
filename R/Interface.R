#' @name RFclust.SGE
#' @aliases RFclust.SGE,RFclust.SGE-method
#' @rdname RFclust.SGE-methods
#' @docType methods
#' @description  create a new RFclust.SGE object This object is mainly used for plotting the
#' @description  data
#' @param dat data frame or matrix containing all expression data
#' @param tmp.path where to store the temporaray files
#' @param SGE whether to use the Sun Grid Engine to calcualate
#' @param slices how many threads to span
#' @return A new RFclust.SGE object
#' @title description of function RFclust.SGE
#' @export 
setGeneric('RFclust.SGE', ## Name
		function ( dat, tmp.path='', email='', slices=32, SGE=FALSE ) { ## Argumente der generischen Funktion
			standardGeneric('RFclust.SGE') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('RFclust.SGE', signature = c ('data.frame'),
		definition = function ( dat, tmp.path='', email='', slices=32, SGE=FALSE ) {
			if ( tmp.path == '' ){
				tmp.path = pwd()
			}
			ret <- new ( 'RFclust.SGE', dat= dat, email=email, tmp.path=tmp.path, slices=slices )
			} 
)

#' @name show
#' @aliases show,RFclust.SGE-method
#' @rdname show-methods
#' @docType methods
#' @description  print the object
#' @param object the RFclust.SGE object
#' @title description of function show
#' @export 
setMethod('show', signature = c ('RFclust.SGE'),
		definition = function ( object ) {
		print ( paste("data frame with",nrow(object@dat),"and",ncol(object@dat),"columns" ))
		if ( object@SGE ){ print ( paste("SGE will be used", "and the SGE will report to", object@email)) }
		else{print ( paste("SGE will NOT be used"))  }
		print (paste( "Number of cores to use:",object@slices ))
		print ( paste("files will be stored in", object@tmp.path))
	}
)

#' @name runRFclust
#' @aliases runRFclust,RFclust.SGE-method
#' @rdname runRFclust-methods
#' @docType methods
#' @description run the random forest calculations returning the density matrix
#' @description at the moment without SGE support and single core
#' @param x the RFclust.SGE object
#' @param ntree the number of trees to grow
#' @param nforest the nuber of forests to create
#' @title description of function runRFclust
#' @export 
setGeneric('runRFclust',
		function (x, ntree=500, nforest=500 ){
			standardGeneric('runRFclust')
		}
)
setMethod('runRFclust', signature = c ('RFclust.SGE'),
		definition = function ( x, ntree=500, nforest=500  ) {
			datRF = calculate.RF(x@dat,  ntree=ntree, no.rep=nforest )
			distRF = RFdist( datRF ,x@dat, imp=TRUE , no.tree=ntree )
			distRF
		}
)


