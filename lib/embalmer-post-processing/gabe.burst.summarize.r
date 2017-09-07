# taken from gabe - 8/10/17
summarize.taxa.embalmer <- function(taxafn, L, outputfn)
{
    taxa = read.delim(taxafn,row=1);

    split = strsplit(rownames(taxa),";");             # Split and rejoin at desired level
    taxaStrings = sapply(split,function(x) paste(x[1:L],collapse=";"));
    for (i in 1:7) taxaStrings = gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T) # clean tips
    taxa = rowsum(taxa,taxaStrings);                  # Collapse by taxonomy name
    # Uncomment the next line to add NA's up to 'L' levels; replace L with desired level
    rownames(taxa) = sapply(strsplit(rownames(taxa),";"),function(x) paste(x[1:L],collapse=";"));
    taxa = sweep(taxa,2,colSums(taxa),'/');
    sink(outputfn); cat("#Taxonomy\t");
    write.table(taxa,file=outputfn,quote=F,sep="\t",append = T);
    sink(NULL)
}

library('optparse')
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
