# after generating imputation Rdata files, run this to actually transform the tables 
# takes in input and output folders - assumes file names are hard coded

require('optparse')
require('robCompositions')
require('vegan')

# do clr transform, samples in rows, taxa in cols
transform.clr <- function(datafn, samplenamesfn)
{
    load(file = datafn)
    taxa.clr = cenLR(imputed$x)$x.clr # clr tranformation of imputed table
    rownames(taxa.clr) <- read.table(samplenamesfn, sep="\t", colClasses="character")[,1]
    taxa.clr
}

#inputfolder should contain all imputed data files and taxa names
transform.clr.folder <- function(inputfolder)
{  
    taxa.clr.list <- NULL
    for(i in 1:6) 
    {
        taxa.clr.list[[i]] <- transform.clr(datafn = paste0(inputfolder, "/impute.L", i+1 ,".rdata"), 
                                samplenamesfn = paste0(inputfolder, "/sample.names.txt"))
        outputfn <- paste0(inputfolder, "/taxa.clr.L", i+1 ,".txt")
        sink(outputfn); 
        cat("#Taxonomy\t");
        write.table(t(taxa.clr.list[[i]]),file=outputfn,quote=F,sep="\t",append = T);
        sink(NULL)
    }
    invisible(taxa.clr.list)
}    


option_list <- list(
 		make_option(c("-i", "--inputfolder"), type="character", 
 			help="taxa file [REQUIRED]")
        )

 	opts <- parse_args(OptionParser(option_list=option_list), 
 		args=commandArgs(trailing=TRUE))
 
    inputfolder <- opts$inputfolder

    transform.clr.folder(inputfolder)