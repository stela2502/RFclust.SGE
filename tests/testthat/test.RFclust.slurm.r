
rf <- RFclust.SGE ( dat=dat, SGE=F, slices=3, slurm=T, settings=list( 'A' = "lsens2017-3-2", t="00:20:00" ,p="dell" ) )

rf <- runRFclust ( rf, nforest=5)

for ( i in 1:3 ) {
	fn <- paste( 'runRFclust_RFrun_',i,'.sh', sep="" )
	expect_true ( file.exists( fn ), info=fn)
}

## this will break everywhere but on aurora-ls2.lunarc.lu.se

#expect_equal(dim(rf@distRF[[1]]$cl1),c(100,100))

#groups <- createGroups( rf, c( 2,4,6,8), name="RFrun"  )


#expect_equal( dim(groups), c(100,6))


#rf <- RFclust.SGE ( dat=dat, SGE=F, slices=2 )
#rf <- runRFclust ( rf, nforest=5 )
#groups <- createGroups( rf, c( 2,4,6,8) )
#
#
#rf <- RFclust.SGE ( dat=dat, SGE=T, slices=2, email='someMail@somewhere.se' )
#rf <- runRFclust ( rf, nforest=5 )
#groups <- createGroups( rf, c( 2,4,6,8) )
