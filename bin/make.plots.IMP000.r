datadirname <- "datadir_denovo"
is.phylo <- TRUE # data supports phylogenetic measures

setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("bin/load.r")

# set up data
    map.000 <- map[map$Subject.ID == "IMP.000",]
    map.000$Sample.Day.Since.First.Sample <- as.numeric(as.Date(map.000$Sample.Date, format="%m/%d/%y") - as.Date("08/01/16", format="%m/%d/%y")) # hard code first sample date
    map.000[,"travel.phase"] <- "Traveling"
    map.000[map.000$Sample.Day.Since.First.Sample <= 2,"travel.phase"] <- "Pre"
    map.000[map.000$Sample.Day.Since.First.Sample >= 28,"travel.phase"] <- "Post"
    map.000$travel.phase <- factor(map.000$travel.phase, levels=c("Pre", "Traveling", "Post")) 

# alpha
    multiplot.alphadiv.L(map.000, alphadiv, alpha.metrics, "alphadiv.IMP000.pdf", x.var="Sample.Day.Since.First.Sample")
    
# TAXA summaries
    ssamples <- map[map$Subject.ID == "IMP.000", "Sample.Order", drop=F]
    ssamples.ordered <- rownames(ssamples[ order(ssamples[,1]), , drop=F]) 
    day.labels <- as.character(sort(ssamples[,1])-3)
    day.labels[(length(day.labels)-3):length(day.labels)] <- paste0("P",1:4)
    plot.taxa.summary(map0=map[ssamples.ordered,], otu=taxa, fn="taxa.summary.IMP000.pdf", sample.order=ssamples.ordered,
                            x.labels = day.labels)
    plot.taxa.summary.L(taxa0=taxa, map0=map.000, x.var="Sample.Day.Since.First.Sample", outputfn="taxa.summary.IMP000.pdf")
    
# PCOA
    plot.pcoa.long(map, rownames(map[map$Subject.ID == "IMP.000",]), bc_dm, "BC - IMP.000")

# RELATIVE
    # IMP000 to Day 0    
    # relative to first day
    multiplot.relative.L(mains=mains, dms=dms, map0=map.000, ylab="Distance to Self at Day 0", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref.sample.order = 1, outputfn="IMP000_to_Day0.pdf", override.cols.by = "travel.phase")
    # overall variability -- relative to previous sample
    multiplot.relative.L(mains=mains, dms=dms, map0=map.000, ylab="Distance to Previous Sample", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref.sample.order = 1, outputfn="IMP000_to_Previous.pdf", to.previous=T, override.cols.by="travel.phase")                
    # relative to Hmong Thai
    multiplot.relative.L(mains=mains, dms=dms, map0=map.000, ylab="Distance to Hmong Thai", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref_samples=hmongthai, outputfn="IMP000_to_HmongThai.pdf", override.cols.by="travel.phase")
    # relative to 2nd Gen Hmong
    multiplot.relative.L(mains=mains, dms=dms, map0=map.000, ylab="Distance to 2nd Gen Hmong", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref_samples=hmong_secondgen_cs, outputfn="IMP000_to_Hmong2nd.pdf", override.cols.by="travel.phase")


# ----------- CLR

    # IMP000 to Day 0    
    # relative to first day
    p <- plot.relative.L(map0 = map.000, otu0=taxa_clr_L7, main="Euclidean", ref.sample.order = 1, xlab="Day", 
                    ylab="Distance to Self at Day 0", x.var="Sample.Day.Since.First.Sample", override.cols.by = "travel.phase")
    save_plot("IMP000_to_Day0.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    # overall variability -- relative to previous sample
    p <- plot.relative.L(map0 = map.000, otu0=taxa_clr_L7, main="Euclidean", ref.sample.order = 1, xlab="Day", 
                    ylab="Distance to Self at Day 0", x.var="Sample.Day.Since.First.Sample", override.cols.by = "travel.phase", to.previous=T)
    save_plot("IMP000_to_Previous.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    # relative to Hmong Thai
    p <- plot.relative.L(map0 = map.000, otu0=taxa_clr_L7, main="Euclidean", ref_samples=hmongthai, xlab="Day", 
                    ylab="Distance to Hmong Thai", x.var="Sample.Day.Since.First.Sample", override.cols.by = "travel.phase")
    save_plot("IMP000_to_HmongThai.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
    
    # relative to 2nd Gen Hmong
    p <- plot.relative.L(map0 = map.000, otu0=taxa_clr_L7, main="Euclidean", ref_samples=hmong_secondgen_cs, xlab="Day", 
                    ylab="Distance to 2nd Gen Hmong", x.var="Sample.Day.Since.First.Sample", override.cols.by = "travel.phase")
    save_plot("IMP000_to_Hmong2nd.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
