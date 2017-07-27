#Simple way to process the embalmer (embalmulate) output for this study:

#### Loading ####
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/sequences_062717")


summarize.taxa.embalmer("taxatable.txt", 7, "taxatable_L7.txt")
summarize.taxa.embalmer("taxatable.txt", 6, "taxatable_L6.txt")
summarize.taxa.embalmer("taxatable.txt", 2, "taxatable_L2.txt")

# level = 7 for species, 6 for genus, 2 for phyla 
summarize.taxa.embalmer <- function(taxafn, level, outputfn)
{
    taxa <- read.delim(taxafn,row=1)
    
    split <- strsplit(rownames(taxa),";")           
    taxaStrings <- sapply(split,function(x) paste(x[1:level],collapse=";"))
    taxaStrings <- gsub("NA", "Other", taxaStrings) 
    taxa <- rowsum(taxa,taxaStrings) # rowsum and group by taxastrings

    cat("#Taxonomy\t", file=outputfn)
    write.table(taxa,file=outputfn,quote=F,sep="\t",append = T)
}

format.otu.table <- function(otufn, outputfn)
{
    otu <- read.delim(otufn,row=1)
    old.names <- rownames(otu)
    
    old.names <- sub(" ", "*", old.names) #replace only the first space with special char
    split.names <- strsplit(old.names, split="*", fixed=T)
    split.df <- data.frame(matrix(unlist(split.names), ncol=2, byrow=T))
    rownames(otu) <- sub("_", " ", split.df[,1])
    otu <- cbind(otu, taxonomy=split.df[,2])
    
    cat("#OTU ID\t", file=outputfn)
    write.table(otu, file=outputfn, quote=F, sep="\t", append=T)
    
}