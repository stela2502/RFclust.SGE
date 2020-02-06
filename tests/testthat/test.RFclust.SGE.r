
rf <- RFclust.SGE ( dat=dat[,1:99], SGE=F, slices=1 )

# rf <- RFclust.SGE ( dat=cellexalObj@data, SGE=F, slices=1 )
#rf <- runRFclust ( rf, nforest=5)
#groups <- createGroups( rf, c( 2,4,6,8), name="RFrun" )

#expect_equal( dim(groups), c(100,6))

#rf <-  RFclust.SGE ( dat=dat, SGE=F, slices=2, name="SLURM", settings = list( A = "lsens2017-3-2", t="00:20:00" , slurm=T) )
rf <- runRFclust ( rf, nforest=500, ntree=500)

#rf <- RFclust.SGE ( dat=dat, SGE=F, slices=2 )
#rf <- runRFclust ( rf, nforest=5 )
g = c( 2,4,6,8,12)
groups <- createGroups( rf, g )
expect_equal( dim(groups), c(99,length(g)+2),label='group dimensions')

expect_equal( as.vector(apply( groups[,3:length(g)+2], 2, function(x) length(unique(x)) )), g, label="group complexity")

rf <- runRFclust ( rf, nforest=500, ntree=500, name="TEST")
groups2 <- createGroups( rf, g, name="TEST")

complexity = 7

A = groups[,complexity]
B = groups2[,complexity]

RandG <- function( a, b , n=1000) {
	lapply( 1:n, function(N){
		max(table( sample(a) ,sample(b)))
	})
}

t=table( A, B )
max = max(t)
randMax = unlist(RandG( A, B ))
p = length(which( randMax  >= max )) / 1000

expect_true( p < 3/1000 ,label="group overlap is fine for random data")



if ( FALSE ){
	rf <- RFclust.SGE ( dat=cellexalObj@data, SGE=F, slices=1 )

	start.time <- Sys.time()
	rf <- runRFclust ( rf, nforest=50, ntree=500, name="TEST")
	end.time <- Sys.time()

	time.taken <- end.time - start.time

	time.taken
	#Time difference of 6.941505 mins

	groups <- createGroups( rf, 40 , name="TEST")

	start.time <- Sys.time()
	rf <- runRFclust ( rf, nforest=50, ntree=500)
	end.time <- Sys.time()

	time.taken <- end.time - start.time

	time.taken
	#3.669606 mins

	groups2 <- createGroups( rf, 40 )

	A= groups[,3]
	B= groups2[,3]

	t=table( A, B )
	max = max(t)
	randMax = unlist(RandG( A, B ))
	p = length(which( randMax  >= max )) / 1000
	expect_true( p < 3/1000 ,label="group overlap is fine for random data")
}