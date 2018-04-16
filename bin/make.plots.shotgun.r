datadirname <- "shotgun"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("bin/load.r")

### heatmap - look only at BP strains
    ret <- load.data(mapfile, otufile=taxa_L7_fn, normalize=F)
    taxa_L7 <- ret$otu
    bpnames <- c(grep("k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides.*",colnames(taxa_L7), value=T),
                grep("k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella.*",colnames(taxa_L7), value=T))
          
    make.heatmap(taxa_L7[,bpnames], map, .25, presence.absence=T, baseline.groups="HmongThai", outputfn="heatmap.h.BP.pdf") 
    make.heatmap(taxa_L7[,bpnames], map, .25, presence.absence=T, baseline.groups="Control", outputfn="heatmap.c.BP.pdf") 
