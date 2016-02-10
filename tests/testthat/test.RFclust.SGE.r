
rf <- RFclust.SGE ( dat=dat, SGE=F, slice=1 )
rf <- runRFclust ( rf, nforest=5 )
groups <- createGroups( distRF, rf@dat, c( 2,4,6,8) )


rf <- RFclust.SGE ( dat=dat, SGE=F, slice=2 )
rf <- runRFclust ( rf, nforest=5 )
groups <- createGroups( rf, c( 2,4,6,8) )
