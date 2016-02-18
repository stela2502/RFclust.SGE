#' @name RFclust.SGE
#' @aliases RFclust.SGE,RFclust.SGE-method
#' @rdname RFclust.SGE-methods
#' @docType methods
#' @description  create a new RFclust.SGE object. The clustering will be performed on the columns of the data.frame.
#' @param dat data frame or matrix containing all expression data
#' @param tmp.path where to store the temporaray files
#' @param SGE whether to use the Sun Grid Engine to calcualate
#' @param slices how many threads to span
#' @return A new RFclust.SGE object
#' @title description of function RFclust.SGE
#' @export 
setGeneric('RFclust.SGE', ## Name
		function ( dat, tmp.path='', email='', slices=32, SGE=FALSE, name='RFclust' ) { ## Argumente der generischen Funktion
			standardGeneric('RFclust.SGE') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('RFclust.SGE', signature = c ('data.frame'),
		definition = function ( dat, tmp.path='', email='', slices=32, SGE=FALSE, name='RFclust' ) {
			if ( tmp.path == '' ){
				tmp.path = pwd()
			}
			if ( length(grep( '^/', tmp.path, perl=T)) == 0 ){
				stop( 'I need the absolute path for the temp path' )
			}
			if ( ! file.exists(tmp.path)){
				dir.create( tmp.path )
			}
			if ( SGE && email=='') {
				stop( "If you plan to use SGE I need an email from you!" )
			}
			ret <- new ( 'RFclust.SGE', dat= dat, email=email, tmp.path=tmp.path, slices=slices, SGE=SGE )
			
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
		if ( length(object@distRF) > 0) {
			print ( paste ("A total of",length(object@distRF),"different anaysis have been run:") )
			for ( i in 1:length(object@distRF ) ){
				print ( names(object@distRF)[i] )
			}
		}
		if ( length(object@RFfiles) > 0) {
			print ( paste ("Running analysis for ",length(object@RFfiles),"analysis runs:") )
			for ( i in 1:length(object@RFfiles ) ){
				print ( names(object@RFfiles)[i] )
			}
		}
	}
)

#' @name runRFclust
#' @aliases runRFclust,RFclust.SGE-method
#' @rdname runRFclust-methods
#' @docType methods
#' @description run the random forest calculations returning the density matrix
#' @description the clusters will be created for the columns of the data.frame
#' @param x the RFclust.SGE object
#' @param ntree the number of trees to grow
#' @param nforest the nuber of forests to create
#' @param name the name of the random forest clustering run (if you want to run multiple)
#' @title description of function runRFclust
#' @return a distRF object to be analyzed by pamNew
#' @export 
setGeneric('runRFclust',
		function (x, ntree=500, nforest=500, name="RFrun", force=FALSE ){
			standardGeneric('runRFclust')
		}
)
setMethod('runRFclust', signature = c ('RFclust.SGE'),
		definition = function ( x, ntree=500, nforest=500, name="RFrun", force=FALSE ) {
			## the most simple - one core no whistles
			run = TRUE
			if ( ! is.null(x@RFfiles[[name]])) {
				## OK - check if they are done and summarize the results
				notDone=FALSE
				for ( f in x@RFfiles[[name]] ){
					if ( locked(f) ) {
						notDone = TRUE
						break
					}
				}
				if ( notDone ){ stop( "Process has not finished!") }
				datRF = read.RF( x, name, 20 )
				x@distRF[[length(x@distRF) +1 ]] = RFdist( datRF ,t(x@dat), imp=TRUE , no.tree=ntree )
				names(x@distRF)[length(x@distRF) ] = name
				x@RFfiles[[name]] <- NULL
				run = FALSE
			}
			else if ( ! is.null( x@distRF[[ name ]] ) ) {
				stop( "This analysis has already been performed! - change name")
			}
			else {
				if ( x@slices == 1 && ! x@SGE ) {
					datRF = calculate.RF(data.frame(t(x@dat)),  no.tree=ntree, no.rep=nforest )
					x@distRF[[length(x@distRF) +1 ]] = RFdist( datRF ,t(x@dat), imp=TRUE , no.tree=ntree )
					names(x@distRF)[length(x@distRF) ] = name
				}
				else {
					## (1) create the RF object file
					srcObj = paste(sep='/', x@tmp.path,paste('RFclust.SGE.',x@name,'.RData', sep='')  )
					save( x, file= srcObj)
					## (2) create and run x@slices worker files
					this.forests = round(nforest/x@slices )
					scripts = vector('character', length= x@slices )
					for ( i in 1:x@slices ) {
						ret <- writeRscript( x, paste('runRFclust',name,i,sep='_'), ntree=ntree, nforest=this.forests,srcObj=srcObj, run = !x@SGE )
						if ( x@SGE ){
							writeSGEscript( x, paste('runRFclust',name,i,sep='_'), ret$cmd )
						}
						scripts[i] <- ret$data
					}
					x@RFfiles[[ length(x@RFfiles) +1 ]] <- scripts
					names(x@RFfiles)[ length(x@RFfiles) ] = name
					print (paste( name, ": The data is going to be analyszed now - re-run this function to check if the process has finished."))
					## now the data should become anayzed - re-running this function should then cluster the data
				}
			}
			x
		}
)



#' @name writeRscript
#' @aliases writeRscript,RFclust.SGE-method
#' @rdname writeRscript-methods
#' @docType methods
#' @description run the random forest calculations returning the density matrix
#' @description at the moment without SGE support and single core
#' @param filename the filename to save the R script to (has to be unique for the analysis!)
#' @param ntree the number of trees to grow
#' @param nforest the nuber of forests to create
#' @title description of function writeRscript
#' @return a filename for the expected data
#' @export 
setGeneric('writeRscript',
		function (x,filename, ntree=500, nforest=500, run=TRUE, srcObj ){
			standardGeneric('writeRscript')
		}
)
setMethod('writeRscript', signature = c ('RFclust.SGE'),
		definition = function ( x,filename, ntree=500, nforest=500, run=TRUE, srcObj  ) {
			#print ( paste( "Run =",run)) 
			wp <- paste(sep='/', x@tmp.path, filename )
			rscript <-  paste(wp, '.R', sep='')
			Rdata <-  paste(wp, '.RData', sep='')
			fileConn<-file( rscript )
			writeLines ( c( 'library(RFclust.SGE)', 
							paste('set.lock("',Rdata,'")',sep=''),
							paste('load("',srcObj,'")', sep='' ),
							'#reads object x',
							paste('datRF = calculate.RF(data.frame(t(x@dat)),  no.tree=',ntree,', no.rep=',nforest,' )'),
							paste('save( datRF, file="',Rdata,'")', sep='' ),
							paste('release.lock("',Rdata,'")',sep='')
					), con=fileConn )
			close(fileConn)
			cmd <- paste('R CMD BATCH --no-save --no-restore --no-readline --', rscript ) 
			if ( run ) {
				system( paste(cmd,"&" ) )
			}
			list( data=Rdata, script=rscript, cmd=cmd )
		}
)


#' @name writeSGEscript
#' @aliases writeSGEscript,RFclust.SGE-method
#' @rdname writeSGEscript-methods
#' @docType methods
#' @description run the random forest calculations returning the density matrix
#' @description at the moment without SGE support and single core
#' @param x the RFclust.SGE object
#' @param filename the base filename for the script (path and .sh will be added!)
#' @param cmd the command to include in the SGE script. Make sure, that all path entries are valid on the nodes!  
#' @title description of function writeSGEscript
#' @return a distRF object to be analyzed by pamNew
#' @export 
setGeneric('writeSGEscript',
		function ( x, filename, cmd ){
			standardGeneric('writeSGEscript')
		}
)
setMethod('writeSGEscript', signature = c ('RFclust.SGE'),
		definition = function ( x,filename, cmd  ) {
			wp <- paste(sep='/', x@tmp.path, filename )
			script <- paste(wp, '.sh', sep='')
			fileConn<-file( script )
			writeLines ( c("#!/bin/bash",
			"#$ -l mem_free=1G",
			"#$ -S /bin/bash",
			paste("#$ -M",x@email),
			"#$ -m eas" ,"#$ -pe orte 1",cmd
			), con=fileConn )
			close(fileConn)
		#	print ( script )
			system( paste("qsub",script) )
		}
)


#' @name createGroups
#' @aliases 'createGroups,RFclust.SGE-method
#' @rdname 'createGroups-methods
#' @docType methods
#' @description get a grouping table from the distRF object
#' @param x RFclust.SGE object after a call to runRFclust()
#' @param k the number of expected groupings or a vector of expected groupings
#' @param name the name of the rf analysis
#' @title description of function randomForest
#' @return a distRF object to be analyzed by pamNew
#' @export 
setGeneric('createGroups',
		function ( x, k, name='RFrun' ){
			standardGeneric('createGroups')
		}
)
setMethod('createGroups', signature = c ('RFclust.SGE'),
		definition = function (x, k, name='RFrun' ) {
			n = k[1]
			persistingCells <- colnames( x@dat )
			res = pamNew(x@distRF[[name]]$cl1, n )
			N <- names( res )
			N <- intersect( persistingCells, N )
			userGroups <- matrix(ncol=3, nrow=0)
			for ( a in 1:length(N) ){
				userGroups <- rbind (userGroups, c( N[a], 'no info', as.numeric(res[[N[a]]]) ) )
			}
			for ( i in 2:length(k) ) {
				if ( i > 1) {
					res = pamNew(x@distRF[[name]]$cl1, k[i] )
					n <- vector('numeric', length= length(N))
					for ( a in 1:length(N) ){
						n[a] <- as.numeric(res[[N[a]]])
					}
					userGroups <- cbind( userGroups, n )
				}
			}
			colnames(userGroups) <- c('cellName', 'userInput',  paste ('group n=', k) )
			userGroups
		}
)


