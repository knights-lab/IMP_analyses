library('optparse')

# level = 7 for species, 6 for genus, 2 for phyla 
summarize.taxa.embalmer <- function(taxafn, level, outputfn)
{
    taxa <- read.delim(taxafn,row=1)
    
    split <- strsplit(rownames(taxa),";")           
    taxaStrings <- sapply(split,function(x) paste(x[1:level],collapse=";"))
    taxaStrings <- gsub("NA", "", taxaStrings) # remove NAs
    taxa <- rowsum(taxa,taxaStrings) # rowsum and group by taxastrings

    cat("#Taxonomy\t", file=outputfn)
    write.table(taxa,file=outputfn,quote=F,sep="\t",append = T)
}

#make sure to turn underscores to spaces in otutable and get rid of all text after the ID itself. also there seems to be a blank for anything that was interpolated as something between archaea and bacteria — so fill that in with something like k_transkingdom instead
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

# make option list and parse command line
 	option_list <- list(
 		make_option(c("-i", "--taxafile"), type="character", 
 			help="taxa file to be summarized [REQUIRED]")
 			)

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
 
    taxafile <- opts$taxafile

summarize.taxa.embalmer(taxafile, 7, gsub(".txt","_L7.txt",taxafile))
summarize.taxa.embalmer(taxafile, 6, gsub(".txt","_L6.txt",taxafile))
summarize.taxa.embalmer(taxafile, 2, gsub(".txt","_L2.txt",taxafile))

