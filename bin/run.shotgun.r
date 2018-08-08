# change these next two lines to run on your local machine
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

datadirname <- "shotgun"
source("bin/load.r")

library(tidyr)
library(car)

### BCOV by sample

    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.top5.pdf", max.features=5, per="all")
    #with sample-labels
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bpgenomes-coverage-per-sample.txt", 
                                rescale="none", outputfn="bubble.BP.top5-xlabs.pdf", max.features=5, per="all", show.samplenames=TRUE)
    # cut off at 10% mean coverage overall
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.mincov10.pdf", min.coverage=.10, per="all")
    # cut off at 25% mean coverage per group
    # this gives roughly the same as the overall 10% coverage, but smaller subset
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.mincov25group.pdf", min.coverage=.25, per="group")

    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.top3.per.group.pdf", max.features=3, min.coverage=.1, per="group")

    # useful for allowing thai to contribute bacteroides even if the coverage is low
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.top2.per.person.pdf", max.features=2, min.coverage=.4, per="person", baseheight=8)

    # useful for cutting off only high coverage genomes
    # note that this results in nearly all thai groups with 0 contributions to bacteroides
    ### USE this for manuscript
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.cov50.per.person.pdf", max.features=NULL, min.coverage=.5, per="person", baseheight=4)

    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                                rescale="none", outputfn="bubble.BP.cov50.per.person.xlab.pdf", max.features=NULL, min.coverage=.5, per="person", baseheight=4, show.samplenames=T)


### plot DBCAN heatmap

    hmm.sample.df <- read.table("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/project_043_shotgun/CAzyme/sample.x.cazyme.txt", sep="\t", head=T, row=1)

    make.heatmap.traditional(otu0=hmm.sample.df, map0=map, outputfn="heatmap.cazymes.q0.05.pdf", 
                             min.prev=0.1, sig.level=.05, show.rownames=T, rescale="zero-mean")
    make.heatmap.traditional(otu0=hmm.sample.df, map0=map, outputfn="heatmap.cazymes.q0.10.pdf", 
                             min.prev=0.1, sig.level=.10, show.rownames=T, rescale="zero-mean")
    # remove longitdinal
    make.heatmap.traditional(otu0=hmm.sample.df, map0=map[map$Sample.Group != "Karen1st",], outputfn="heatmap.cazymes.HT-H1-C.q0.05.pdf", 
                             min.prev=0.1, sig.level=.05, show.rownames=T, rescale="zero-mean", filter.mode="nonparametric")

    make.heatmap.traditional(otu0=hmm.sample.df, map0=map[map$Sample.Group != "Karen1st",], outputfn="heatmap.cazymes.HT-H1-C.q0.01.pdf", 
                             min.prev=0.1, sig.level=.01, show.rownames=T, rescale="zero-mean", filter.mode="nonparametric")

    # make sure to comment back in get.color.list in heatmap for longitudinal
    # nothing significantly different here
    #     mapk1 <- map[map$Sample.Group == "Karen1st",]
    #     mapk1$Sample.Month <- as.factor(mapk1$Sample.Month)
    #     make.heatmap.traditional(otu0=hmm.sample.df, map0=mapk1, group.var="Sample.Month", outputfn="heatmap.cazymes.K1.q0.05.pdf", 
    #                              min.prev=0.1, sig.level=.05, show.rownames=T, rescale="zero-mean", filter.mode="nonparametric", paired=TRUE)


        
### take a look at worms database
    ret <- load.data(mapfile, otufile="data/shotgun/worms/imp_worms_taxa.txt", normalize=F)
    worms.otu <- ret$otu
    
    # rel abundance of OTUs that hit
    worms.otu <- sweep(worms.otu, 1, rowSums(worms.otu), '/')

    # filter out rare OTUs
	worms.otu <- worms.otu[, colMeans(worms.otu) > .001, drop=FALSE]

    # filter by prevalence
    prevalences <- apply(worms.otu, 2, function(bug.col) mean(bug.col > 0))
    worms.otu <- worms.otu[, prevalences >= 0.10]

    ret <- collapse.by.correlation(worms.otu, .95)
    worms.otu <- worms.otu[, ret$reps]

    make.heatmap.traditional(otu0=worms.otu, map0=map, outputfn="heatmap.worms.pdf", min.prev=0.01, sig.level=1, show.rownames=T, rescale="standard", filter.mode="none")
    make.heatmap.traditional(otu0=worms.otu, map0=map, outputfn="heatmap.worms2.pdf", min.prev=0.01, sig.level=1, show.rownames=T, rescale="zero-mean", filter.mode="none")


### BCOV by entire dataset
#     bcovbin <- read.table("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/project_043_shotgun/burst-humarine-plus/bcov98unique_binary.txt", comment="", sep="\t", head=T)
#     colnames(bcovbin) <- c("id","coverage")
# 
#     bcovbin <- bcovbin[order(bcovbin$coverage),]
# 
#     bcovbin_80 <- bcovbin[bcovbin$coverage > .8,]
# 
#     bp <- grep("Bacteroides|Prevotella", bcovbin_80$id, value=T)
# 
#     ret <- load.data(mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", normalize=F)
#     otu <- ret$otu
# 
#     otu_bp80 <- otu[,bp]
# 
#     colnames(otu_bp80) <- shorten.taxonomy(colnames(otu_bp80))
# 
#     make.heatmap.traditional(otu0=otu_bp80, map0=map, outputfn="heatmap.trad.BP80.pdf",min.prev=0.01, sig.level=1, show.rownames=T)


### filter by high coverage based on Ben's coverage python code
# 
#         b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bpgenomes-coverage-per-sample.txt", 
#                                     rescale="none", outputfn="heatmap.bubble.BPnorm.top-overall.pdf", max.features=5, min.coverage.per.group=NULL)
#         # with sample-labels
#         b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bpgenomes-coverage-per-sample.txt", 
#                                     rescale="none", outputfn="heatmap.bubble.BPnorm.top-overall-labels.pdf", max.features=5, min.coverage.per.group=NULL)
# 
#         b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bpgenomes-coverage-per-sample.txt", 
#                                     rescale="standard", outputfn="heatmap.bubble.BP.top-overall.pdf", max.features=5, min.coverage.per.group=NULL)

### heatmap - look only at BP strains -- this doesn't take into account coverage - NIX
#     ret <- load.data(mapfile, otufile="data/shotgun/humarine-plus/taxatable.txt", normalize=F)
#     otu <- ret$otu
#     
#     # prevalence only      
#     #make.heatmap(taxa_L7[,bpnames], map, .25, presence.absence=T, baseline.groups="HmongThai", outputfn="heatmap.h.BP.pdf") 
#     #make.heatmap(taxa_L7[,bpnames], map, .25, presence.absence=T, baseline.groups="Control", outputfn="heatmap.c.BP.pdf") 
# 
#     # rel abundance of OTUs that hit
#      otu <- sweep(otu, 1, rowSums(otu), '/')
# 
#     prevalences <- apply(otu, 2, function(bug.col) mean(bug.col > 0))
#     otu <- otu[, prevalences >= 0.10]
# 
#     ret <- collapse.by.correlation(otu, .95)
#     
#     otu <- otu[, ret$reps]
# 
#     bpnames <- c(grep("k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides.*",colnames(otu), value=T),
#                 grep("k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella.*",colnames(otu), value=T))
# 
#     bptaxa <- otu[,bpnames]
#     colnames(bptaxa) <- shorten.taxonomy(colnames(bptaxa))
#     make.heatmap.traditional(otu0=bptaxa, map0=map, outputfn="heatmap.trad.BP.pdf",min.prev=0.01, sig.level=1, show.rownames=T)