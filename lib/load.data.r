# Loads a map and otu file (taxa or pathway), normalizes the OTU table and returns the OTU table, 
# map and kegg descriptions (if any) in order of the remaining samples
#
# mapfile: mapping metadata file
# otufile: OTU table in classic format
# normalize: normalizes, transforms, then collapses OTU table; set this to false when working with merged L1,..,L6 taxa files 
# minOTUInSamples: drop OTUs in less than this ratio of samples (usually .001)
# returns otu, map, and kegg vector containing kegg descriptions (named by whatever KEGG level was passed in)
load.data<-function(mapfile, otufile, minOTUInSamples=NA, minPrevalence = .10, normalize=TRUE)
{
	source("/Users/pvangay/Dropbox/UMN/Rscripts/collapse-features.r")
	require(biom)
	
	map <- read.table(mapfile,sep='\t',head=T,row=1,comment='')
#	rownames(map) <- make.names(rownames(map)) # this is automatically done for otu table, so do this here to make sure sample ids are exactly the same 

	fnlength <- nchar(otufile)
	if(substr(otufile, fnlength-4, fnlength) == ".biom") {
		otu0 <- as.matrix(biom_data(read_biom(otufile)))
	} else {
		line1 <-readLines(otufile,n=1)
		if(line1=="# Constructed from biom file") {
			otu0 <-read.table(otufile,sep='\t',head=T,row=1,comment='',quote="",skip=1)
		} else {
			otu0 <-read.table(otufile,sep='\t',head=T,row=1,comment='',quote="")
		}
	}

	# pathway files from picrust usually have an additional last column with descriptions, just drop this for now
	KEGG <- NULL
	if(colnames(otu0)[ncol(otu0)]=="KEGG_Pathways"){
		KEGG <- setNames(otu0[, ncol(otu0)], rownames(otu0))
		otu0 <- otu0[,-ncol(otu0)]
	}
	
	otu <- t(otu0)

	if(normalize==TRUE){
		 otu <- sweep(otu, 1, rowSums(otu), '/')

        if(!is.na(minOTUInSamples)){
     		otu <- otu[, colMeans(otu) > minOTUInSamples, drop=FALSE]
     	}

 		prevalences <- apply(otu, 2, function(bug.col) mean(bug.col > 0))
 		otu <- otu[, prevalences >= minPrevalence]

		otu <- asin(sqrt(otu))

		ret <- collapse.by.correlation(otu, .95)
		otu <- otu[, ret$reps]
	
	}
	full <- merge(otu, map, by=0)	
	otu <- otu[full$Row.names,]
	map <- map[full$Row.names,]
	
	if(length(KEGG) > 0){
		KEGG <- KEGG[colnames(otu)] #use the same order as our final otu table
	}
	
	list(otu=otu, map=map, kegg=KEGG)
}
