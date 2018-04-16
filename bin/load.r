source(paste0(LIBDIR,"alphadiv.r"))
source(paste0(LIBDIR,"boxplot.by.group.x.bmi.r"))
source(paste0(LIBDIR,"relative.distance.r"))
source(paste0(LIBDIR,"fraction.hits.r"))
source(paste0(LIBDIR,"load.data.r"))
source(paste0(LIBDIR,"b.p.ratio.r"))
source(paste0(LIBDIR,"group.distances.r"))
source(paste0(LIBDIR,"nutrients.r"))
source(paste0(LIBDIR,"pcoa.color.by.nutrients.r"))
source(paste0(LIBDIR,"alphadiv.r"))
source(paste0(LIBDIR,"pcoa.r"))
source(paste0(LIBDIR,"body.trends.r"))
source(paste0(LIBDIR,"predict.years.r"))
source(paste0(LIBDIR,"utils.r"))
source(paste0(LIBDIR,"nutrient.mb.heatmap.r"))
source(paste0(LIBDIR,"taxa.summary.r"))
source(paste0(LIBDIR,"differential.taxa.r"))
source(paste0(LIBDIR,"heatmap.r"))

## set file locations
    datadirs <- list()
    #     datadirs[["datadir_gg97"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/gg97"
    #     datadirs[["datadir_dada2"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/dada2"
    #     datadirs[["datadir_refseq"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/refseq"
    #     datadirs[["datadir_openref"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/openref"
    #     datadirs[["datadir_akronymer"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/akronymer" 
    datadirs[["ref"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/refseq_allruns"
    datadirs[["denovo"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/denovo"
    datadirs[["shotgun"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/shotgun"
    datadirs[["denovo.clr"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/denovo/clr"
    datadirs[["shotgun.clr"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/shotgun/clr"
    datadirs[["food"]] <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/food"

    config <- data.frame(datadir=matrix(unlist(datadirs), nrow=length(datadirs), byrow=T),stringsAsFactors=FALSE, row.names=names(datadirs))
    config$is.phylo <- FALSE
    config[c("denovo","ref"), "is.phylo"] <- TRUE # these are the only two with phylogenetic metrics
    config$is.clr <- FALSE
    config$is.clr[grep("clr$", rownames(config))] <- TRUE

    is.phylo <- config[datadirname, "is.phylo"]
    is.clr <- config[datadirname, "is.clr"]
    datadir <- config[datadirname, "datadir"]

    food_datadir <- datadirs[["food"]]

    taxa_L7_fn <- paste(datadir,"taxatable_L7.txt",sep="/")
    taxa_L6_fn <- paste(datadir,"taxatable_L6.txt",sep="/")
    taxa_L2_fn <- paste(datadir,"taxatable_L2.txt",sep="/")

    otu_fn <- paste(datadir,"final_otu.txt",sep="/")

    wuf_dm_fn <- paste(datadir,"weighted_unifrac_dm.txt",sep="/")
    uwuf_dm_fn <- paste(datadir,"unweighted_unifrac_dm.txt",sep="/")
    bc_dm_fn <- paste(datadir,"bray_curtis_dm.txt",sep="/")

    alpha_fn <- paste(datadir,"alpha.txt",sep="/")

    nutrientsfn <- paste(food_datadir,"nutrients.txt",sep="/")
    foodgroupfn <- paste(food_datadir,"foodgroups.txt",sep="/")
    dietmap_fn <- paste(food_datadir,"sampleid-to-dietid.txt",sep="/")
    
    mapfile <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mapping.txt"

## set constants
    bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
    prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"
    alpha.phylo <- c("PD_whole_tree", "observed_otus", "shannon") 
    alpha.nonphylo <- c("chao1", "observed_otus",  "shannon")
    if(is.phylo) alpha.metrics <- alpha.phylo else alpha.metrics <- alpha.nonphylo

# load data
    # only load these files if this is not a CLR dataset
    if(!is.clr){
        alphadiv <- read.table(alpha_fn, sep="\t", quote="", row=1, head=T, comment.char="")
        bc_dm <- read.table(bc_dm_fn, sep="\t", quote="", row=1, head=T)
        otu <- load.data(mapfile, otufile=otu_fn, normalize=F)$otu
        if(is.phylo){
            wuf_dm <- read.table(wuf_dm_fn, sep="\t", quote="", row=1, head=T)
            uwuf_dm <- read.table(uwuf_dm_fn, sep="\t", quote="", row=1, head=T)
        }
    }
    # for all taxa and otu, always skip normalizing, filtering, or collapsing - do these only when needed
    # taxa and otus have already had singletons removed (taxa/otus present in only 1 sample, i.e. chimeras)
    taxa_L7 <- load.data(mapfile, otufile=taxa_L7_fn, normalize=F)$otu
    taxa_L2 <- load.data(mapfile, otufile=taxa_L2_fn, normalize=F)$otu
    ret <- load.data(mapfile, otufile=taxa_L6_fn, normalize=F)
    map_all <- ret$map
    taxa_L6 <- ret$otu
    taxa <- taxa_L6 # save main taxa and map files based on L6
        
# reformat mapping variables
    # only work with samples that are in both the mapping AND the dm/alpha files
    # turns out that 1 sample that has been sequenced should be excluded (abx)
    map <- map_all[map_all$Exclude != "Y",]

    source("bin/format.map.r")
    
# add food metadata to mapping 
    # let's add some diet metadata into the mapping file 
    diet_map <- read.table(dietmap_fn, sep="\t", header=T, check.names=F, as.is=T,row=1)
    merged <- merge(map, diet_map[,c("SuperTracker.DATE", "Diet.Month", "Diet.ID", "Diet.Date")], by=0)
    # check on rows
    #cat("map: ", nrow(map), " diet_map: ", nrow(diet_map), " merged: ", nrow(merged), "\n")
    
    map <- data.frame(merged[,-1], row.names=merged[,1])

    # calculate Days since arrival for sample and diet dates
    map$Sample.Day.Since.Arrival <- as.numeric(as.Date(map$Sample.Date, format="%m/%d/%y") - as.Date(map$Arrival.in.US, format="%m/%d/%y"))
    map$Diet.Day.Since.Arrival <- as.numeric(as.Date(map$Diet.Date, format="%m/%d/%y") - as.Date(map$Arrival.in.US, format="%m/%d/%y") )
    map$Sample.Day.Since.Arrival <- as.numeric(as.Date(map$Sample.Date, format="%m/%d/%y") - as.Date(map$Arrival.in.US, format="%m/%d/%y"))
 
    # note that we expect n=4 samples filtered out of regular taxa/otu files (too low of depth)
    # check map_all for the original unfiltered mapping file 
    # these contain sample ids found in diet but not in mapping -- note that most of these were samples that were not sequenced for whatever reason
    # check <- data.frame(Sample.ID=setdiff(rownames(diet_map), rownames(map)), Status=NA)
    # check[grep("^[IL].*", check$Sample.ID), "Status"] <- "No Diet or No Sample"
    # check[grep("^TFS.*", check$Sample.ID), "Status"] <- "Low Depth"
    # check[check$Sample.ID=="T.CS.060", "Status"] <- "Low Depth"
    # check[check$Sample.ID=="T.CS.046", "Status"] <- "Excluded"