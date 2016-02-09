
rf <- RFclust.SGE ( dat=dat, SGE=F, slice=1 )
distRF <- runRFclust ( rf, nforest=5 )
