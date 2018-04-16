# Run this on teraminx (takes forever) - Imputes zeros for CLR transformation
# CLR cannot deal with tree metrics yet. therefore for this script we'll only use the taxa table for everything.

library('optparse')
require('robCompositions')


# impute zeros for taxa (samples in rows, taxa in cols)
impute.clr <- function(taxa, datafn)
{
    taxa <- taxa*1.0 
    # dl can be 1 or 2 (this means that every cell in taxa needs to be greater than this number. filter first)
    # maxit maximum likelihood iterations that we'll need it 
    
    # save the taxa table for debug purposes only
    #write.table(taxa, paste0(datafn,".inputtaxa.txt"), sep="\t", quote=F)
    
    imputed = impRZilr(taxa, maxit = 3, method = "lm", dl = rep(2, ncol(taxa)), verbose = T) 
    save(imputed, file = datafn)
}

# remove singleton OTUs (OTUs that appear in only 1 sample) and low depth samples
# taxa.count = integer specifying the min number of samples an OTU must appear in (1 = drop singleton OTUs)
# min.sample.depth = integer specifying min depth per sample (inclusive)
filter.taxa <- function(taxa, min.prevalence=0, sample.depth=0)
{
    taxa <- taxa[rowSums(taxa) >= sample.depth, ]
    prevalences <- apply(taxa, 2, function(bug.col) mean(bug.col > 0))
    taxa <- taxa[, prevalences > min.prevalence]
    
    # for the purposes of speeding up the clr imputation, remove all bugs that have prevalence < 10%
    taxa
}

summarize.taxa <- function(taxa, L)
{
    split = strsplit(rownames(taxa),";");             # Split and rejoin at desired level
    taxaStrings = sapply(split,function(x) paste(x[1:L],collapse=";"));
    for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
    taxa = rowsum(taxa,taxaStrings);                  # Collapse by taxonomy name
    # Uncomment the next line to add NA's up to 'L' levels; replace L with desired level
    # rownames(taxa) = sapply(strsplit(rownames(taxa),";"),function(x) paste(x[1:L],collapse=";"));
    # taxa = sweep(taxa,2,colSums(taxa),'/');
    taxa
}



 	option_list <- list(
 		make_option(c("-t", "--taxafile"), type="character", 
 			help="taxa file [REQUIRED]"),
        make_option(c("-f", "--outputfolder"), type="character", 
            help="full path to output folder [REQUIRED]"),
        make_option(c("-d", "--sampledepth"), type="integer", 
            help="minimum sample depth [REQUIRED]")
 			)

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
 
    taxafile <- opts$taxafile
    outfolder <- opts$outputfolder
    sample.depth <- opts$sampledepth   # sample.depth = 4584
    
    # read in files
    taxa0 <- read.delim(taxafile, row=1)
        
    # filter taxa table - transpose
    taxa.f <- filter.taxa(t(taxa0), min.prevalence=.10, sample.depth=sample.depth)
    
    for(L in 2:7)
    {
        # summarize taxa table
        this.taxa <- t(summarize.taxa(t(taxa.f), L))
        
        # after summarizing, some taxa could still be filled with only 0s or 1s - drop these
        this.taxa.f <- this.taxa[, !(apply(this.taxa, 2, function(v) all(v < 2)))]

        impute.clr(this.taxa.f, datafn = paste0("impute.L", L ,".rdata"))
    }
    # save sample names for later because we lose them during the imputation step
    write.table(rownames(this.taxa), file="sample.names.txt", sep="\t", quote=F, col.names=F, row.names=F)



    

    
    
    
    
    
    
    
    
    
    
    