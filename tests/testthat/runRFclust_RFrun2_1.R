library(RFclust.SGE)
set.lock("/home/med-sal/git_Projects/RFclust.SGE/tests/testthat//runRFclust_RFrun2_1.RData")
load("/home/med-sal/git_Projects/RFclust.SGE/tests/testthat//.RData")
#reads object x
datRF = calculate.RF(data.frame(t(x@dat)),  no.tree= 500 , no.rep= 5  )
datRF = RFdist( datRF ,t(x@dat), imp=TRUE , no.tree= 500  )
save( datRF, file="/home/med-sal/git_Projects/RFclust.SGE/tests/testthat//runRFclust_RFrun2_1.RData")
release.lock("/home/med-sal/git_Projects/RFclust.SGE/tests/testthat//runRFclust_RFrun2_1.RData")
