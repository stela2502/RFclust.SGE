
system( 'rm -f *.RData' )
rf2 <- RFclust.SGE ( dat=dat, SGE=F, slices=2 )
rf2 <- runRFclust ( rf2, nforest=10, name="RFrun2")

files = rf2@RFfiles[["RFrun2"]]

locked = length(files)

while ( locked != 0 ) {
	locked = 0
	for ( f in files ){
		if (locked( f )) {
			locked = locked +1;
		}
	}
	Sys.sleep( locked )
}

data <- read.RF( rf2, 'RFrun2', 20 )
expect_equal ( round(data$err1, 3), 0.282 )
expect_equal ( dim(data$RFproxAddcl1), c(100,100))

rf2 <- runRFclust ( rf2, nforest=10, name="RFrun2")

expect_equal(dim(rf2@distRF[[1]]$cl1),c(100,100))


groups <- createGroups( rf2, c( 2,4,6,8), name="RFrun2" )

expect_equal( dim(groups), c(100,6))


#rf2 <- RFclust.SGE ( dat=dat, SGE=F, slices=2 )
#rf2 <- runRFclust ( rf, nforest=5 )
#groups <- createGroups( rf, c( 2,4,6,8) )
#
#
#rf2 <- RFclust.SGE ( dat=dat, SGE=T, slices=2, email='someMail@somewhere.se' )
#rf2 <- runRFclust ( rf, nforest=5 )
#groups <- createGroups( rf, c( 2,4,6,8) )
