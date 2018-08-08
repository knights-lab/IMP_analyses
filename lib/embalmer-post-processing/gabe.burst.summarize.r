# taken from gabe - 8/10/17
summarize.taxa.embalmer <- function(taxafn, L, outputfn, do.counts=FALSE)
{
    line1 <-readLines(taxafn,n=1)
    if(line1=="# Constructed from biom file") skip <- 1
    else skip <- 0

    taxa = read.delim(taxafn,row=1, skip=skip);

    split = strsplit(rownames(taxa),";");             # Split and rejoin at desired level
    taxaStrings = sapply(split,function(x) paste(x[1:L],collapse=";"));
    for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
    taxa = rowsum(taxa,taxaStrings);                  # Collapse by taxonomy name
    # Uncomment the next line to add NA's up to 'L' levels; replace L with desired level
#    rownames(taxa) = sapply(strsplit(rownames(taxa),";"),function(x) paste(x[1:L],collapse=";"));
    if(!do.counts) taxa = sweep(taxa,2,colSums(taxa),'/');
    sink(outputfn); cat("#Taxonomy\t");
    write.table(taxa,file=outputfn,quote=F,sep="\t",append = T);
    sink(NULL)
}

library('optparse')
# make option list and parse command line
 	option_list <- list(
 		make_option(c("-i", "--taxafile"), type="character", 
 			help="taxa file to be summarized [REQUIRED]"),
        make_option(c("-a", "--absolute"), default=FALSE, 
 			help="keep taxa table in absolute counts [OPTIONAL]"),
        make_option(c("-o", "--outputfile"), type="character", default=NULL,
 			help="outputfile path with prepended file name [OPTIONAL]")
 			)

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
 
    taxafile <- opts$taxafile
    do.counts <- opts$absolute
    outfile <- opts$outputfile

for(level in c(2,6,7)){
    ext <- paste0("_L", level, ".txt")
    summarize.taxa.embalmer(taxafile, level, ifelse(is.null(outfile), gsub(".txt", ext, taxafile), paste0(outfile, ext)), do.counts=do.counts)
}
# summarize.taxa.embalmer(taxafile, 7, ifelse(is.null(outfile), gsub(".txt","_L7.txt",taxafile), paste0(outfile, "_7.txt")), do.counts=do.counts)
# summarize.taxa.embalmer(taxafile, 6, gsub(".txt","_L6.txt",taxafile), do.counts=do.counts)
# summarize.taxa.embalmer(taxafile, 2, gsub(".txt","_L2.txt",taxafile), do.counts=do.counts)
