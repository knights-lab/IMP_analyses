# change these next two lines to run on your local machine
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

# load and format all data
source(paste0(LIBDIR,"mouse.load.data.r")) 

### PCOA
    map_mb_last <- map_mb[map_mb$Week %in% c(8,10),]

    map_pcoa <- map_mb[-which(map_mb$Mouse.ID %in% c("M35","M36","M37","M38")),]
    map_pcoa_last <- map_pcoa[map_pcoa$Week %in% c(8,10),]

    map_pcoa_nobase1_end <- map_pcoa[map_pcoa$Donor %in% c("TFSCS023","IMP.263"),]    
    map_pcoa_nobase1_last <- map_pcoa_last[map_pcoa_last$Donor %in% c("TFSCS023","IMP.263"),]
    map_pcoa2_end <- map_pcoa[map_pcoa$Donor %in% c("TFSCS026","IMP.264"),]
    map_pcoa2_last <- map_pcoa_last[map_pcoa_last$Donor %in% c("TFSCS026","IMP.264"),]

    # ALL Points
    plot.mouse.pcoa(map_mb, dm=uwuf_dm, outputfn="pcoa.uwuf.pdf", type="allpoints")
    plot.mouse.pcoa(map_mb, dm=wuf_dm, outputfn="pcoa.wuf.pdf", type="allpoints")

    # AlL endpoints    
    plot.mouse.pcoa(map0=map_mb, dm=uwuf_dm, outputfn="pcoa.endpoints.uwuf.pdf", type="endpoints")
    plot.mouse.pcoa(map0=map_mb, dm=wuf_dm, outputfn="pcoa.endpoints.wuf.pdf",  type="endpoints")

    # LAST POINTS
    plot.mouse.pcoa(map_mb_last, dm=uwuf_dm, outputfn="pcoa.uwuf.last.pdf", type="allpoints")
    plot.mouse.pcoa(map_mb_last, dm=wuf_dm, outputfn="pcoa.wuf.last.pdf", type="allpoints")


    # ALL POINTS, without infected mouse
    plot.mouse.pcoa(map_pcoa,  dm=uwuf_dm, outputfn="pcoa.uwuf.ex.pdf", type="allpoints")
    plot.mouse.pcoa(map_pcoa, dm=wuf_dm, outputfn="pcoa.wuf.ex.pdf", type="allpoints")

    # LAST POINTS, without infected mouse
    plot.mouse.pcoa(map_pcoa_last, dm=uwuf_dm, outputfn="pcoa.uwuf.last.ex.pdf", type="allpoints")
    plot.mouse.pcoa(map_pcoa_last, dm=wuf_dm, outputfn="pcoa.wuf.last.ex.pdf", type="allpoints")

    # endpoints only, without infected mouse
    plot.mouse.pcoa(map_pcoa, dm=uwuf_dm, outputfn="pcoa.endpoints.uwuf.ex.pdf", type="endpoints")
    plot.mouse.pcoa(map_pcoa, dm=wuf_dm, outputfn="pcoa.endpoints.wuf.ex.pdf",  type="endpoints")


    # Donor Pair 1 - ENDPOINTS ONLY
    plot.mouse.pcoa(map_pcoa_nobase1_end, dm=uwuf_dm, outputfn="pcoa.pair1.endpoints.uwuf.pdf", type="endpoints")
    plot.mouse.pcoa(map_pcoa_nobase1_end, dm=wuf_dm, outputfn="pcoa.pair1.endpoints.wuf.pdf",  type="endpoints")

    # Donor Pair 1 - last points only
    plot.mouse.pcoa(map_pcoa_nobase1_last, dm=uwuf_dm, outputfn="pcoa.pair1.lastpoints.uwuf.pdf", type="allpoints")
    plot.mouse.pcoa(map_pcoa_nobase1_last, dm=wuf_dm, outputfn="pcoa.pair1.lastpoints.wuf.pdf",  type="allpoints")

    # Donor Pair 1 - last points WITH DONORS
    plot.mouse.pcoa.donors(map_pcoa_nobase1_last, dm=uwuf_dm, outputfn="pcoa.pair1.lastpoints.donor.uwuf.pdf")
    plot.mouse.pcoa.donors(map_pcoa_nobase1_last, dm=wuf_dm, outputfn="pcoa.pair1.lastpoints.donor.wuf.pdf")
 
    plot.mouse.pcoa.donors(map_pcoa_nobase1_last[map_pcoa_nobase1_last$Donor=="TFSCS023",], dm=uwuf_dm, outputfn="pcoa.pair1.lastpoints.thaidonor.uwuf.pdf")
    plot.mouse.pcoa.donors(map_pcoa_nobase1_last[map_pcoa_nobase1_last$Donor=="TFSCS023",], dm=wuf_dm, outputfn="pcoa.pair1.lastpoints.thaidonor.wuf.pdf")

    plot.mouse.pcoa.donors(map_pcoa_nobase1_last[map_pcoa_nobase1_last$Donor=="IMP.263",], dm=uwuf_dm, outputfn="pcoa.pair1.lastpoints.usdonor.uwuf.pdf")
    plot.mouse.pcoa.donors(map_pcoa_nobase1_last[map_pcoa_nobase1_last$Donor=="IMP.263",], dm=wuf_dm, outputfn="pcoa.pair1.lastpoints.usdonor.wuf.pdf")

    # Donor Pair 2 - ENDPOINTS ONLY
    plot.mouse.pcoa(map_pcoa2_end, dm=uwuf_dm, outputfn="pcoa.pair2.endpoints.uwuf.pdf", type="endpoints")
    plot.mouse.pcoa(map_pcoa2_end, dm=wuf_dm, outputfn="pcoa.pair2.endpoints.wuf.pdf",  type="endpoints")

    # Donor Pair 2 - last points only
    plot.mouse.pcoa(map_pcoa2_last, dm=uwuf_dm, outputfn="pcoa.pair2.lastpoints.uwuf.pdf", type="allpoints")
    plot.mouse.pcoa(map_pcoa2_last, dm=wuf_dm, outputfn="pcoa.pair2.lastpoints.wuf.pdf",  type="allpoints")

    # Donor Pair 2 - last points WITH DONORS
    plot.mouse.pcoa.donors(map_pcoa2_last, dm=uwuf_dm, outputfn="pcoa.pair2.lastpoints.donor.uwuf.pdf")
    plot.mouse.pcoa.donors(map_pcoa2_last, dm=wuf_dm, outputfn="pcoa.pair2.lastpoints.donor.wuf.pdf")

    # Donor Pair 2 - last points WITH DONORS, with infected
    plot.mouse.pcoa.donors(map_mb_last[map_mb_last$Donor %in% c("TFSCS026","IMP.264"),], dm=uwuf_dm, outputfn="pcoa.pair2.last.infected.donor.uwuf.pdf")
    plot.mouse.pcoa.donors(map_mb_last[map_mb_last$Donor %in% c("TFSCS026","IMP.264"),], dm=wuf_dm, outputfn="pcoa.pair2.last.infected.donor.wuf.pdf")

    # Donor Pair 1 - ENDPOINTS ONLY
    plot.mouse.pcoa(map_mb[map_mb$Cohoused==F,], dm=uwuf_dm, outputfn="test.pdf", type="endpoints")

        # SKIP, this looks too messy
        #Donor Pair 1 - PATH
        #plot.mouse.pcoa(map_pcoa, samples=rownames(map_pcoa), dm=uwuf_dm, outputfn="pcoa.pair1.uwuf.pdf", type="path")
        #plot.mouse.pcoa(map_pcoa, samples=rownames(map_pcoa), dm=wuf_dm, outputfn="pcoa.pair1.wuf.pdf", type="path")

### check taxa summary for cage swap errors for M1-M10 at Week 8
    test.map <- map_mb[map_mb$Donor %in% c("TFSCS023","IMP.263") & map_mb$Diet.Type=="LowFiber" & map_mb$Week %in% 1:8, ]
    test.map$Mouse.ID <- factor(test.map$Mouse.ID, levels=paste0("M",1:10))
    plot.taxa.summary.L(taxa0=taxa, map0=test.map, outputfn="M1-10.taxa.summary.pdf", max.taxa=30, x.var="Week", subject.var="Mouse.ID", grid.ncol=5)

### check taxa summary for infection blooms
    test2.map <- map_mb[which(map_mb$Mouse.ID %in% c("M35","M36","M37","M38")),]
    test2.map <- map_mb[map_mb$Donor == "TFSCS026",]
    plot.taxa.summary.L(taxa0=taxa_L7, map0=test2.map, outputfn="infected.mice.taxa.summary.pdf", max.taxa=15, x.var="Week", subject.var="Mouse.ID", grid.ncol=4)


### VILLI & CRYPT
    plot.vc(merge(villus_crypt, map[map$timepoint=="endpoint",], by="Mouse.ID"), outputfn="villus-crypt-ratio.pdf")
    plot.v.or.c(merge(villus_crypt, map[map$timepoint=="endpoint",], by="Mouse.ID"), type="v", outputfn="villus.pdf")        

    plot.vc(merge(villus_crypt, map[!(map$Mouse.ID %in% c("M35","M36","M37","M38")) & map$Cohoused==F & map$timepoint=="endpoint",], by="Mouse.ID"), outputfn="villus-crypt-ratio-noncohoused-noninfected.pdf")

    plot.vc(merge(villus_crypt, map[map$Cohoused==F & map$timepoint=="endpoint",], by="Mouse.ID"), outputfn="villus-crypt-ratio-non-cohoused.pdf")
    plot.v.or.c(villus.crypt.df=merge(villus_crypt, map[map$Cohoused==F & map$timepoint=="endpoint",], by="Mouse.ID"), type="v", outputfn="villus-non-cohoused.pdf")        
    plot.v.or.c(villus.crypt.df=merge(villus_crypt, map[map$Cohoused==F & map$timepoint=="endpoint",], by="Mouse.ID"), type="c", outputfn="crypt-non-cohoused.pdf")        

### MUCUS THICKNESS


    mucus.n <- normalize.mucus(mucus)

    mucus.df <- merge(mucus.n, map[map$timepoint=="endpoint",], by="Mouse.ID")
 
    mucus.df <- merge(mucus.df, unique(mucus[,c("Mouse.ID","Reembedded")]), by="Mouse.ID")
    
   
    p.reemb.all <- plot.mucus(ggdata=mucus.df[mucus.df$Week==8,], highlight.reembedded=TRUE)

    p.reemb.two <- plot.mucus(ggdata=mucus.df[mucus.df$Group.End %in% c("Thai.LowFiber","US.HighFiber") & mucus.df$Week==8,], highlight.reembedded=TRUE)
    
    p.mucus.all <- plot.mucus(ggdata=mucus.df)
 
    # remove suspected infected mouse
    mucus.noninf.df <- mucus.df[-which(mucus.df$Mouse.ID %in% c("M35","M36","M37","M38")),]
    
    p.mucus.noninf <- plot.mucus(ggdata=mucus.noninf.df[mucus.noninf.df$Week==8,])
    save_plot("mucus-noninfected.pdf", p.mucus.noninf.all)
    
    p.mucus.noninf.all <- plot.mucus(ggdata=mucus.noninf.df)
    save_plot("mucus-noninfected-allgroups.pdf", p.mucus.noninf.all)
   
    # try donor 1 only 
    plot.mucus(mucus.df[mucus.df$Donor %in% c("TFSCS023","IMP.263"),], outputfn="mucus-donor1.pdf")
    
### WEIGHT
    plot.weights(ggdata=map, group="Group.End")
   
### GLUCOSE
    plot.glucose(map[!is.na(map$Fasting.Glucose),], group="Group.End", add.pval=TRUE, outputfn="mouse-glucose-foldchange.pdf")

### FOOD CONSUMPTION
    plot.food(map, group="Group.Start")  
    
### FEED EFFICIENCY - percent weight gain by food consumed
    # A high FE means that they are more efficient users of food --> less food required to gain weight
    plot.feed.efficiency.L(map, group="Group.Start")
    plot.feed.efficiency.wk8(map, group="Group.End", outputfn="feed-efficiency-endpoint.pdf")
    
### body composition
    colnames(body_comp) <- gsub("\\.\\.g\\.", "", colnames(body_comp))
    body_comp <- read.table(file=body_comp_fn, sep="\t", header=T, as.is=T)    
    gg_bodycomp <- merge(body_comp[,c("Mouse.ID","Percent.Fat","Percent.Lean")], map[map$timepoint=="endpoint",], by="Mouse.ID")
    
    p_lean <- mouse.boxplot(y=gg_bodycomp$Percent.Lean, Group=gg_bodycomp$Group.End, 
            main="% Lean", facet.var=NULL, add.pval=TRUE,
            ylab="", strip.text.size=10, group.vars.df=gg_bodycomp[,c("Donor.Type","Diet.Type")])
    save_plot("mouse-percent-lean.pdf", p_lean, base_aspect_ratio = 1)
    p_fat <- mouse.boxplot(y=gg_bodycomp$Percent.Fat, Group=gg_bodycomp$Group.End, 
            main="% Fat", facet.var=NULL, add.pval=TRUE,
            ylab="", strip.text.size=10,group.vars.df=gg_bodycomp[,c("Donor.Type","Diet.Type")])
    save_plot("mouse-percent-fat.pdf", p_fat, base_aspect_ratio = 1)

### plot all variations of immune cells

    ## last TH group removed (infected)
    immune <- read.table(file=immune_fn, sep="\t", header=T, as.is=T, check.names=F, quote="", comment="")
    immune_apc <- read.table(file=immune_apc_fn, sep="\t", header=T, as.is=T, check.names=F, quote="", comment="")
    immune_live <- read.table(file=immune_live_fn, sep="\t", header=T, as.is=T, check.names=F, quote="", comment="")

    # remove suspected infected mouse
    immune <- immune[-which(immune$Mouse.ID %in% c("M35","M36","M37","M38")),]
    
    # merge percent live with both full immune and apc files
    immune.live <- merge(immune, immune_live, by="Sample")

    # get rid of the old one since it's not really the Live, rename
    immune.live[,"Percent.Live.CD45+.x"] <- immune.live[,"Percent.Live.CD45+.y"]
    immune.live <- immune.live[,-which(colnames(immune.live)=="Percent.Live.CD45+.y")]
    immune.live <- immune.live[,-which(colnames(immune.live)=="Percent.Live.Cells")]
    colnames(immune.live) <- gsub("\\.[xy]", "", colnames(immune.live))
    
    # using % dead-and-alive CD45 counts
        # all groups w/ pvals
        make.immune.plots(ggdata=map[map$timepoint=="endpoint" & map$Date != max(map$Date),], 
                immune=immune, immune_apc=immune_apc, 
                group.var="Group.End", grouping.vars=c("Donor.Type","Diet.Type","Cohoused"),
                add.pval=TRUE, outdir="immune-cells-all/")
        
        # non-cohoused groups only w/ pvals
        make.immune.plots(ggdata=map[map$timepoint=="endpoint" & map$Date != max(map$Date) & map$Cohoused==FALSE,], 
                immune=immune, immune_apc=immune_apc, 
                group.var="Group.End", grouping.vars=c("Donor.Type","Diet.Type"),
                add.pval=TRUE, outdir="immune-cells-noncohoused/")
    
    # using live CD45 counts only
        # all groups w/ pvals
        make.immune.plots(ggdata=map[map$timepoint=="endpoint" & map$Date != max(map$Date),], 
                immune=immune.live, immune_apc=immune_apc, 
                group.var="Group.End", grouping.vars=c("Donor.Type","Diet.Type","Cohoused"),
                add.pval=TRUE, outdir="immune-cells-all-live/")
        
        # non-cohoused groups only w/ pvals
        make.immune.plots(ggdata=map[map$timepoint=="endpoint" & map$Date != max(map$Date) & map$Cohoused==FALSE,], 
                immune=immune.live, immune_apc=immune_apc, 
                group.var="Group.End", grouping.vars=c("Donor.Type","Diet.Type"),
                add.pval=TRUE, outdir="immune-cells-noncohoused-live/")






















