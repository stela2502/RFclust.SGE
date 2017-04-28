
source("https://bioconductor.org/biocLite.R")
biocLite()

install.packages( 
		lib  = lib <- .libPaths()[1],
		pkgs = as.data.frame(installed.packages(lib), stringsAsFactors=FALSE)$Package,
		type = 'source',  
		repos='https://ftp.acc.umu.se/mirror/CRAN/'
)
biocLite("BiocUpgrade") 

if (!library("devtools", quietly = TRUE,logical.return=TRUE )) {
	install.packages(c('devtools'),  repos='https://ftp.acc.umu.se/mirror/CRAN/')
	library(devtools)
}