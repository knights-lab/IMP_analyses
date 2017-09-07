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
    imputed = impRZilr(taxa, maxit = 3, method = "lm", dl = rep(2,ncol(taxa)), verbose = T) 
    save(imputed, file = datafn)
}

# remove singleton OTUs (OTUs that appear in only 1 sample) and low depth samples
# taxa.count = integer specifying the min number of samples an OTU must appear in (1 = drop singleton OTUs)
# min.sample.depth = integer specifying min depth per sample (inclusive)
filter.taxa <- function(taxa, min.prevalence, sample.depth)
{
    taxa <- taxa[rowSums(taxa) >= sample.depth, ]
    prevalences <- apply(taxa, 2, function(bug.col) mean(bug.col > 0))
    taxa <- taxa[, prevalences > min.prevalence]
    
    # for the purposes of speeding up the clr imputation, let's remove all bugs that 
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
            help="minimum sample depth [REQUIRED]"),
        make_option(c("-m", "--method"), type="character", 
            help="[impute, transform] [REQUIRED]")
 			)

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
 
    taxafile <- opts$taxafile
    outfolder <- opts$outputfolder
    sample.depth <- opts$sampledepth   # sample.depth = 4584
    method <- opts$method
    
    # read in files
    taxa0 <- read.delim(taxafile, row=1)
        
    # filter taxa table - transpose
    taxa.f <- filter.taxa(t(taxa0), min.prevalence=.10, sample.depth=sample.depth)
    
    # summarize taxa table
    taxa.list <- lapply(2:7, function(L) t(summarize.taxa(t(taxa.f), L)))
    
    # impute zeros - this takes a long time
    for(i in 1:6)   impute.clr(taxa.list[[i]], datafn = paste0("impute.L", i+1 ,".rdata"))

    # save sample names for later because we lose them during the imputation step
    write.table(rownames(taxa.list[[1]]), file="sample.names.txt", sep="\t", quote=F, col.names=F, row.names=F)



    

    
    
    
    
    
    
    
    
    
    
    