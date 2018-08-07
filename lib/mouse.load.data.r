library(lme4)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(cowplot)
library(scales)
library(ggsignif)
library(reshape)
library(car)

source(paste0(LIBDIR,"load.data.r"))
source(paste0(LIBDIR,"utils.r"))
source(paste0(LIBDIR,"mouse.immune.cells.r"))
source(paste0(LIBDIR,"mouse.villus.crypt.r"))
source(paste0(LIBDIR,"mouse.pcoa.r"))
source(paste0(LIBDIR,"mouse.utils.r"))
source(paste0(LIBDIR,"mouse.weights.r"))
source(paste0(LIBDIR,"mouse.glucose.r"))
source(paste0(LIBDIR,"mouse.food.r"))

### load & prep all mouse data

    datadir <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse"
    datadirdonors <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse/with_donors"

    # file names
        mapfile <- paste(datadir,"mapping.txt",sep="/")
        mapfile_donors <- paste(datadirdonors,"mapping_donors.txt",sep="/")
        vc_fn <- paste(datadir,"all-cpsr-measurements.txt",sep="/")
        immune_fn <- paste(datadir,"immune-cell-counts.txt",sep="/")
        immune_apc_fn <- paste(datadir,"immune-cell-counts-apc.txt",sep="/")
        immune_live_fn <- paste(datadir,"immune-live-cells.txt",sep="/")
        body_comp_fn <- paste(datadir,"Body Composition.txt",sep="/")
        mucus_fn <- paste(datadir,"mucus.txt",sep="/")

        datadir <- datadirdonors
        taxa_L7_fn <- paste(datadir,"taxatable_L7.txt",sep="/")
        taxa_L6_fn <- paste(datadir,"taxatable_L6.txt",sep="/")
        taxa_L2_fn <- paste(datadir,"taxatable_L2.txt",sep="/")
        otu_fn <- paste(datadir,"final_otu.txt",sep="/")
        wuf_dm_fn <- paste(datadir,"weighted_unifrac_dm.txt",sep="/")
        uwuf_dm_fn <- paste(datadir,"unweighted_unifrac_dm.txt",sep="/")
        bc_dm_fn <- paste(datadir,"bray_curtis_dm.txt",sep="/")
        alpha_fn <- paste(datadir,"alpha.txt",sep="/")

    ## load MB data
        alphadiv <- read.table(alpha_fn, sep="\t", quote="", row=1, head=T, comment.char="")
        bc_dm <- read.table(bc_dm_fn, sep="\t", quote="", row=1, head=T)
        otu <- read.table(otu_fn, sep="\t", quote="", row=1, head=T, comment="")
        wuf_dm <- read.table(wuf_dm_fn, sep="\t", quote="", row=1, head=T)
        uwuf_dm <- read.table(uwuf_dm_fn, sep="\t", quote="", row=1, head=T)

        # for all taxa and otu, always skip normalizing, filtering, or collapsing - do these only when needed
        # taxa and otus have already had singletons removed (taxa/otus present in only 1 sample, i.e. chimeras)
        taxa_L7 <- load.data(mapfile, otufile=taxa_L7_fn, normalize=F)$otu
        taxa_L2 <- load.data(mapfile, otufile=taxa_L2_fn, normalize=F)$otu
        ret <- load.data(mapfile, otufile=taxa_L6_fn, normalize=F)
        taxa_L6 <- ret$otu
        taxa <- taxa_L6 # save main taxa and map files based on L6
        map_mb <- ret$map
        map_mb$Mouse.ID <- as.character(map_mb$Mouse.ID)      
        colnames(map_mb) <- gsub("\\.\\.g\\.", "", colnames(map_mb))
        
    ## load and format the full map file
        # DO NOT use map from loading MB data because only sequenced samples are present!
        map_donor <- read.table(file=mapfile_donors, sep="\t", header=T, as.is=T, quote="",row=1)
        map <- read.table(file=mapfile, sep="\t", header=T, as.is=T, quote="",row=1)
        map$Mouse.ID <- as.character(map$Mouse.ID)      
        colnames(map) <- gsub("\\.\\.g\\.", "", colnames(map))
        
        map$Group.Start <- factor(map$Group.Start, 
                    levels=c("Thai.HighFiber", "Thai.LowFiber", "US.HighFiber", "US.LowFiber"))

        map$Group.End <- factor(map$Group.End, 
                    levels=c(levels(map$Group.Start), c("Thai.HighFiber.Cohoused", "US.HighFiber.Cohoused", "Thai.LowFiber.Cohoused", "US.LowFiber.Cohoused")))

# global vars
GROUP.FILLS <<- c("#865596","#865596", "#F26230","#F26230","#865596","#F26230","#865596","#F26230") # purple orange
GROUP.COLORS <<- c("#865596","#865596", "#F26230","#F26230",rep("black",4)) # purple orange
GROUP.ALPHAS <<- c(rep(.8,4), rep(.5,4))
GROUP.SHAPES <<- c(16,17,16,17,21,21,24,24) #circle and triangle

GROUP.NAMES.SHORT <- c("TH", "TL", "UH", "UL", "THC", "UHC", "TLC", "ULC")
names(GROUP.NAMES.SHORT) <- levels(map$Group.End)

names(GROUP.COLORS) <- levels(map$Group.End)
names(GROUP.FILLS) <- levels(map$Group.End)
names(GROUP.ALPHAS) <- levels(map$Group.End)
names(GROUP.SHAPES) <- levels(map$Group.End)


    ## load other data
        villus_crypt <- read.table(file=vc_fn, sep="\t", header=T, as.is=T)
       
        mucus <- read.table(file=mucus_fn, sep="\t", header=T, as.is=T)