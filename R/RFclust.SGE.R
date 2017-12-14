#' @name Rand
#' @aliases Rand,RFclust.SGE-method
#' @rdname Rand-methods
#' @docType methods
#' @description The function computes the (adjusted) Rand index between two partitions 
#' @description Copyright Steve Horvath and Luohua Jiang, UCLA, 2003
#' @param tab the data table
#' @param adjust just leave it ... default=T
#' @title description of function Rand
setGeneric('Rand', ## Name
	function (tab,adjust=T) { ## Argumente der generischen Funktion
		standardGeneric('Rand') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('Rand', signature = c ('RFclust.SGE'),
	definition = function (tab,adjust=T) {
		
	# helper function
	choosenew <- function(n,k) {
		n <- c(n); out1 <- rep(0,length(n));
		for (i in c(1:length(n)) ){
			if ( n[i]<k ) {out1[i] <- 0}
			else {out1[i] <- choose(n[i],k) }
		}
		out1
	}
	
	a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
	n <- nrow(tab)
	for (i in 1:n) {
		for(j in 1:n) {
			a <- a+choosenew(tab[i,j],2)
			nj <- sum(tab[,j])
			c <- c+choosenew(nj,2)
		}
		ni <- sum(tab[i,])
		b <- b+choosenew(ni,2)
		nn <- nn+ni
	}
	if(adjust==T) {
		d <- choosenew(nn,2)
		adrand <- (a-(b*c/n)/d)/(0.5*(b+c/n)-(b*c/n)/d)
		adrand
	} else {
		b <- b-a
		c <- c/n-a
		d <- choosenew(nn,2)-a-b-c
		rand <- (a+d)/(a+b+c+d)
		rand
	}
} )
#' @name pamNew
#' @aliases pamNew,RFclust.SGE-method
#' @rdname pamNew-methods
#' @docType methods
#' @description A new pam clustering function which corrects the clustering membership based on the sillhouette strength.
#' @description The clustering membership of an observation with a negative sillhouette strength is reassigned to its
#' @description neighboring cluster.   
#' @description The inputs of the function are similar to the original 'pam' function.          
#' @description The function returns a vector of clustering labels.      
#' @description Copyright 2003 Tao Shi and Steve Horvath (last modified 10/31/03)  
#' @param x  the distRF object obtained from RFdist()
#' @param k  number of exprected clusters
#' @param diss1  TEXT MISSING default= inherits(x, "dist")
#' @param metric1  TEXT MISSING default= "euclidean"
#' @title description of function pamNew
#' @export 
setGeneric('pamNew', ## Name
	function (x, k, diss1 = inherits(x, "dist"), metric1 = "euclidean") { ## Argumente der generischen Funktion
		standardGeneric('pamNew') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('pamNew', signature = c ('matrix'),
	definition = function (x, k, diss1 = inherits(x, "dist"), metric1 = "euclidean") {
	
	if (diss1)
	{
		if (!is.null(attr(x, "Labels"))) { original.row.names <- attr(x, "Labels")}
		names(x) <- as.character(c(1:attr(x, "Size")))
	} 
	else
	{
		if(!is.null(dimnames(x)[[1]])) { original.row.names <- dimnames(x)[[1]]}
		row.names(x) <- as.character(c(1:dim(x)[[1]]))
	}
	pam1 <- cluster::pam(x,k,diss=diss1, metric=metric1)
	label2 <- pam1$clustering
	silinfo1 <- pam1$silinfo$widths
	index1 <- as.numeric(as.character(row.names(silinfo1)))
	silinfo2 <- silinfo1[order(index1),]
	labelnew <- ifelse(silinfo2[,3]<0, silinfo2[,2], silinfo2[,1])
	names(labelnew) <- original.row.names
	labelnew    
} )
#' @name collect.garbage
#' @aliases collect.garbage,RFclust.SGE-method
#' @rdname collect.garbage-methods
#' @docType methods
#' @description The following function collects garbage until the memory is clean.
#' @description Usage: 1. immediately call this function after you call a function or
#' @description        2. rm()
#' @title description of function collect.garbage
setGeneric('collect.garbage', ## Name
	function (x) { ## Argumente der generischen Funktion
		standardGeneric('collect.garbage') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('collect.garbage', signature = c ('missing'),
	definition = function (x='A') {

	while (gc()[2,4] != gc()[2,4]){}
} )


#' @name set.lock
#' @aliases set.lock,RFclust.SGE-method
#' @rdname set.lock-methods
#' @docType methods
#' @description set a lock for a file (threading)
#' @param filename  the locked file
#' @title description of function set_lock
#' @export 
setGeneric('set.lock', ## Name
	function ( filename ) { ## Argumente der generischen Funktion
		standardGeneric('set.lock') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('set.lock', signature = c ('character'),
	definition = function ( filename ) {
	system ( paste('touch ',filename,'.lock', sep='') )
} )
#' @name release.lock
#' @aliases release.lock,RFclust.SGE-method
#' @rdname release.lock-methods
#' @docType methods
#' @description releases the lock of a file (threading)
#' @param filename  the locked file
#' @title description of function release_lock
#' @export 
setGeneric('release.lock', ## Name
	function ( filename ) { ## Argumente der generischen Funktion
		standardGeneric('release.lock') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('release.lock', signature = c ('character'),
	definition = function ( filename ) {
	system ( paste('rm ',filename,'.lock', sep='') )
} )
#' @name locked
#' @aliases locked,RFclust.SGE-method
#' @rdname locked-methods
#' @docType methods
#' @description simple check for a file lock (threading)
#' @param filename  lock this file
#' @title description of function locked
#' @export 
setGeneric('locked', ## Name
	function ( filename ) { ## Argumente der generischen Funktion
		standardGeneric('locked') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('locked', signature = c ('character'),
	definition = function ( filename ) {
	ret = TRUE
	if (file.exists( filename )){
		ret <- file.exists( paste(filename,'.lock', sep='') )
	}
	ret
} )

#' @name save.RF
#' @aliases save.RF,RFclust.SGE-method
#' @rdname save.RF-methods
#' @docType methods
#' @description Store a calculate.RF result in a file for later usage (threading)
#' @param Rf.data  the return value from a calculate.RF call
#' @param fname  the file to store the data in
#' @title description of function save.RF
#' @export 
setGeneric('save.RF', ## Name
	function ( Rf.data , fname ) { ## Argumente der generischen Funktion
		standardGeneric('save.RF') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('save.RF', signature = c ('RFclust.SGE'),
	definition = function ( Rf.data , fname ) {
	save( Rf.data, file=fname )
} )


####################################################################
# Unsupervised randomForest function                               #
# Return a list "distRF" containing some of the following 6 fields #
#  depending on the options specified:                             #
#  (1) cl1:  addcl1 distance (sqrt(1-RF.proxAddcl1))               #
#  (2) err1: error rate                                            #
#  (3) imp1: variable importance for addcl1                        #
#  (4) prox1Conver: a matrix containing two convergence meausres   #
#                   for addcl1 proximity                           #
#                   a). max( abs( c(aveprox(N))-c(aveprox(N-1))))  #
#                   b). mean((c(aveprox(N))-c(aveprox(N-1)))^2)    #
#                   where N is number of forests (no.rep).         #
#  (5) cl2, (6) err2, (7)imp2 and (8) prox2Conver for addcl2       #
# Copyright Steve Horvath and Tao Shi (2004)                       #
####################################################################

#' @name calculate.RF
#' @aliases calculate.RF,RFclust.SGE-method
#' @rdname calculate.RF-methods
#' @docType methods
#' @description This function is the threadable random forest function.
#' @description The return of this fucntion should be stored to a file using save.RF
#' @description and later on during the clean up of the different threads be read using read.RF
#' @description The return list can be processed using the RFdist method
#' @param datRF the summary object is a list of random forest return values default= NULL
#' @param mtry1  Number of variables randomly sampled as candidates at each split. default=3
#' @param no.rep  how many times shoul the random forest function be called default= 20
#' @param no.tree  the numer of trees in the random forest function default= 500
#' @param addcl1  synthetic data created by the synthetic1 function default=TRUE
#' @param addcl2  synthetic data created by the synthetic2 function default=FALSE
#' @param imp  whether to calculate the importance default=T
#' @param oob.prox1  leave it like that default=T
#' @param max.syn  max synthetic data width default=50
#' @title description of function calculate.RF
#' @export 
setGeneric('calculate.RF', ## Name
	function ( datRF = NULL, mtry1=3, no.rep= 20, no.tree= 500, addcl1=TRUE, addcl2=FALSE,  imp=T, oob.prox1=T, max.syn=50) { ## Argumente der generischen Funktion
		standardGeneric('calculate.RF') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('calculate.RF', signature = c ('data.frame'),
	definition = function ( datRF = NULL, mtry1=3, no.rep= 20, no.tree= 500, addcl1=TRUE, addcl2=FALSE,  imp=T, oob.prox1=T, max.syn=50) {
	
	synthetic1 <- function(dat, syn.n=NULL) {
		sample1 <- function(X)   { sample(X, replace=T) } 
		g1      <- function(dat) { apply(dat,2,sample1) }
		nrow1 <- dim(dat)[[1]]
		yy <- rep(c(1,2),c(nrow1,nrow1) )
		data.frame(cbind(yy,rbind(dat,data.frame(g1(dat)))))
	}
	
	synthetic1.1 <- function(dat) {
		sample1 <- function(X)   { sample(X, replace=T) } 
		g1      <- function(dat) { apply(dat,2,sample1) }
		nrow1 <- dim(dat)[[1]]
		syn.n <- nrow1
		if ( syn.n > max.syn ) {
			syn.n = max.syn 
		}
		yy <- rep(c(1,2),c(nrow1,syn.n) )
		data.frame(cbind(yy,rbind(dat,data.frame(g1(dat)[1:syn.n,]))))
	}
	
	synthetic2 <- function(dat) {
		sample2 <- function(X)   { runif(length(X), min=min(X), max =max(X)) }
		g2      <- function(dat) { apply(dat,2,sample2) }
		nrow1 <- dim(dat)[[1]];
		yy <- rep(c(1,2),c(nrow1,nrow1) );
		data.frame(cbind(yy,rbind(dat,data.frame(g2(dat)))))
	}
	
	Rf.data <- vector('list', no.rep +1)
	syn.n <- nrow1 <- dim(datRF)[[1]]
	if ( syn.n > max.syn ) {
		syn.n = max.syn
	}
	ncol1 <- dim(datRF)[[2]]
	rep1 <- rep(-1,nrow1+syn.n) 
	
	if ( addcl1 && addcl2 ){
		stop( "Sorry you can not get an addc1 AND add2 distribution on the same time!")
	}
	if (addcl1) {
		for (i in c(0:no.rep))  {
			index1 <- sample(c(1:(nrow1+syn.n))) 
			rep1[index1] <-  c(1:(nrow1+syn.n)) 
			datRFsyn <- synthetic1(datRF,syn.n)[index1,] 
			yy <- datRFsyn[,1]
			RF1 <- randomForest::randomForest(factor(yy)~.,data=datRFsyn[,-1], ntree=no.tree, oob.prox=oob.prox1, proximity=TRUE,do.trace=F,mtry=mtry1,importance=imp)
			collect.garbage()
			RF1prox <- RF1$proximity[rep1,rep1]
			Rf.data[[i+1]] <- list(index1=index1, yy=yy, RF1prox=RF1prox, importance =RF1$importance, err.rate=RF1$err.rate )
		}
	}
	if (addcl2) { 
		for (i in c(0:no.rep))  {
			index1 <- sample(c(1:(2*nrow1))) 
			rep1[index1] <-  c(1:(2*nrow1)) 
			datRFsyn <- synthetic2(datRF)[index1,] 
			yy <- datRFsyn[,1] 
			RF1 <- randomForest::randomForest(factor(yy)~.,data=datRFsyn[,-1], ntree=no.tree, oob.prox=oob.prox1, proximity=TRUE,do.trace=F,mtry=mtry1,importance=imp) 
			collect.garbage()
			RF1prox <- RF1$proximity[rep1,rep1]
			Rf.data[[i+1]] <- list(index1=index1, yy=yy, RF1prox=RF1prox, importance =RF1$importance, err.rate=RF1$err.rate )
		}
	}
	Rf.data
} )

#' @name RFdist
#' @aliases RFdist,RFclust.SGE-method
#' @rdname RFdist-methods
#' @docType methods 
#' @description Initially this function did perform the unsupervised clustering,
#' @description but the most calculation has been exported to the threadable function
#' @description \link{calculate.RF}
#' @param Rf.data  TEXT MISSING
#' @param datRF  TEXT MISSING
#' @param imp  TEXT MISSING default=T
#' @param no.tree  TEXT MISSING
#' @param proxConver  TEXT MISSING default=F
#' @title description of function RFdist
#' @export 
setGeneric('RFdist', ## Name
		function (Rf.data, datRF, imp=T, no.tree, proxConver=F) { ## Argumente der generischen Funktion
			standardGeneric('RFdist') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)


setMethod('RFdist', signature = c ('list'),
		definition = function (Rf.data, datRF, imp=T, no.tree, proxConver=F) {
			
			####################################################################
# Unsupervised randomForest function                               #
# Return a list "distRF" containing some of the following 6 fields #
#  depending on the options specified:                             #
#  (1) cl1:  addcl1 distance (sqrt(1-RF.proxAddcl1))               #
#  (2) err1: error rate                                            #
#  (3) imp1: variable importance for addcl1                        #
#  (4) prox1Conver: a matrix containing two convergence meausres   #
#                   for addcl1 proximity                           #
#                   a). max( abs( c(aveprox(N))-c(aveprox(N-1))))  #
#                   b). mean((c(aveprox(N))-c(aveprox(N-1)))^2)    #
#                   where N is number of forests (no.rep).         #
#  (5) cl2, (6) err2, (7)imp2 and (8) prox2Conver for addcl2       #
# Copyright Steve Horvath and Tao Shi (2004)                       #
			####################################################################
			
			
			no.rep <- length(Rf.data)
			nrow1 <- dim(datRF)[[1]]
			ncol1 <- dim(datRF)[[2]]
			RFproxAddcl1 <- matrix(0,nrow=nrow1,ncol=nrow1)
			RFprox1Conver <- cbind(1:no.rep,matrix(0,(no.rep),3))
			RFimportance1 <- matrix(0, nrow=ncol1, ncol=4)
			RFerrrate1 <- 0
			rep1 <- rep(666,2*nrow1) 
			i = 0;
			while( length(Rf.data) > 0 ) {
				yy <- Rf.data[[1]]$yy
				importance <- Rf.data[[1]]$importance
				err.rate <- Rf.data[[1]]$err.rate
				RF1prox <- Rf.data[[1]]$RF1prox
				if (i > 0) { 
					if (i > 1){
						xx <- ((RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]))/i) - (RFproxAddcl1/(i-1))
						yy <- mean( c(as.dist((RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]))/i))) 
						RFprox1Conver[i,2] <- max(abs(c(as.dist(xx))))
						RFprox1Conver[i,3] <- mean((c(as.dist(xx)))^2)
						RFprox1Conver[i,4] <- yy
					}
					RFproxAddcl1 <- RFproxAddcl1 + (RF1prox[c(1:nrow1),c(1:nrow1)]) 
					if(imp) { RFimportance1 <- RFimportance1+ 1/no.rep*(importance) }
					RFerrrate1 <- RFerrrate1+ 1/no.rep*(err.rate[no.tree])
				}
				Rf.data[[1]] <- NULL
				i = i +1
			}
			
#			cleandist <- function(x) { 
#				x1 <- as.dist(x)
#				x1[x1<=0] <- 0.0000000001
#				as.matrix(x1)
#			}
#			distRFAddcl1 <- cleandist(sqrt(1-RFproxAddcl1/no.rep))
			#distRF$cl1 <- cleandist(sqrt(1-distRF$cl1/no.rep))
			distRF <- list(cl1=NULL, err1=NULL, imp1=NULL, prox1Conver=NULL, RFproxAddcl1 = RFproxAddcl1,
					cl2=NULL, err2=NULL, imp2=NULL, prox2Conver=NULL)
			
			#distRF$cl1 <- distRFAddcl1
			distRF$err1 <- RFerrrate1
			if(imp) distRF$imp1 <- RFimportance1 
			if(proxConver) distRF$prox1Conver <- RFprox1Conver
			
			distRF
		} )
		
#' @name read.RF
#' @aliases read.RF,RFclust.SGE-method
#' @rdname read.RF-methods
#' @docType methods
#' @description This function reads a set of outfiles and creates a summary Rf.data
#' @description object that can be processed using the RFdist() function.
#' @param files  a list of files created by save.RF() default=c('')
#' @param max.wait maximum time to wait for the files to become accessable default= 20
#' @title description of function read.RF
#' @export 
setGeneric('read.RF', ## Name
		function (x, name, max.wait = 20 ) { ## Argumente der generischen Funktion
			standardGeneric('read.RF') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('read.RF', signature = c ('RFclust.SGE'),
		definition = function (x, name, max.wait = 20 ) {
			returnRF <- NULL
			waited = 0
			read = 0
			files <- x@RFfiles[[name]]
			x@RFfiles <- lapply(  x@RFfiles, function( oldF ) { file.path( x@tmp.path, basename(oldF) ) } )
			print (paste ("Reading",length(files),"result files"))
			for ( i in 1:length(files) ){
				if ( ! locked( files[i]) ) {
					if ( i == 1 ){
						load(files[i])
						returnRF <- datRF
					}
					else {
						load(files[i])
						for( n in names( datRF )){
							returnRF[[n]] = returnRF[[n]] + datRF[[n]]
						}
						cat (".")
					}
							
				}else {
					stop( paste("Error: File", files[i], "not finished" ) )
				}
			}
					
			print ( "files read")
			
			
			returnRF$err1 =  returnRF$err1 / length(files)
			returnRF
		} )

#' @name pwd
#' @aliases pwd
#' @rdname pwd-methods
#' @docType methods
#' @description  uses the linux pwd command to determin the working directory 
#' @return A string containing the working directory 
#' @title description of function pwd
setGeneric('pwd', ## Name
		function ( a ) { ## Argumente der generischen Funktion
			standardGeneric('pwd') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('pwd', signature = c () ,
		definition = function ( a ) {
			rm(a)
			system( 'pwd > __pwd' )
			t <- read.delim( file = '__pwd', header=F)
			t <- as.vector(t[1,1])
			t <- paste(t,"/",sep='')
			unlink( '__pwd')
			t
		})
