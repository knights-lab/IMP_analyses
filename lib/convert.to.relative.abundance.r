require(biom)
	
convert.to.relative.abundance <- function(otufile, outputfile)
{
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
	taxonomy <- otu0$taxonomy
    names(taxonomy) <- rownames(otu0)
	
	otu <- t(otu0[,1:(ncol(otu0)-1)])

    otu <- sweep(otu, 1, rowSums(otu), '/')

    # turn all rel abundances to abs counts based on the smallest rel abundance value
    x <- min(otu[otu!=0])    
    otu <- otu*10^(zeros_after_period(x) + 2)
    
    otu <- as.data.frame(t(otu))
    
    otu[,"taxonomy"] <- taxonomy[rownames(otu)]
    
    cat("#OTU ID\t", file=outputfile)
    write.table(otu, outputfile,sep="\t",quote=F, append=T)
}

zeros_after_period <- function(x) {
    if (isTRUE(all.equal(round(x),x))) return (0) # y would be -Inf for integer values
    y <- log10(abs(x)-floor(abs(x)))   
    ifelse(isTRUE(all.equal(round(y),y)), -y-1, -ceiling(y)) # corrects case ending with ..01
}
