
rf <- RFclust.SGE ( dat=dat, SGE=F, slices=1 )

#rf <- runRFclust ( rf, nforest=5)
#groups <- createGroups( rf, c( 2,4,6,8), name="RFrun" )

#expect_equal( dim(groups), c(100,6))

rf <-  RFclust.SGE ( dat=dat, SGE=F, slices=2, name="SLURM",  A = "lsens2017-3-2", t="00:20:00" , slurm=T )
rf <- runRFclust ( rf, nforest=5)

#rf <- RFclust.SGE ( dat=dat, SGE=F, slices=2 )
#rf <- runRFclust ( rf, nforest=5 )
#groups <- createGroups( rf, c( 2,4,6,8) )
#
#
#rf <- RFclust.SGE ( dat=dat, SGE=T, slices=2, email='someMail@somewhere.se' )
#rf <- runRFclust ( rf, nforest=5 )
#groups <- createGroups( rf, c( 2,4,6,8) )
