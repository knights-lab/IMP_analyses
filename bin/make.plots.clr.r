## plot stuff using CLR generated taxa tables

# pcoa with euclidean and bray curtis
# differential taxa

taxa_clr_L7_fn <- paste(datadir,"taxa.clr.L7.txt",sep="/")
taxa_clr_L6_fn <- paste(datadir,"taxa.clr.L6.txt",sep="/")
taxa_clr_L2_fn <- paste(datadir,"taxa.clr.L2.txt",sep="/")
taxa_clr_L3_fn <- paste(datadir,"taxa.clr.L3.txt",sep="/")
taxa_clr_L5_fn <- paste(datadir,"taxa.clr.L5.txt",sep="/")

taxa_clr_L7 <- load.data(mapfile, otufile=taxa_clr_L7_fn, normalize=F)$otu
taxa_clr_L6 <- load.data(mapfile, otufile=taxa_clr_L6_fn, normalize=F)$otu
taxa_clr_L2 <- load.data(mapfile, otufile=taxa_clr_L2_fn, normalize=F)$otu
taxa_clr_L3 <- load.data(mapfile, otufile=taxa_clr_L3_fn, normalize=F)$otu
taxa_clr_L5 <- load.data(mapfile, otufile=taxa_clr_L5_fn, normalize=F)$otu



### PCOA

    # remove IMP.000 from any dm calcs
    plot.pcoa.long(map[map$Subject.ID != "IMP.000",], samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), otu0 = taxa_clr_L7, method="euclidean", plot.title="Euclidean - CLR - L.outlined", convex.hull=TRUE)

    # let's flip the axis so that Years in US goes left to right, BMI/Age goes bottom to up
    plot.pcoa(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR", flip.axis=2)    

    # this isn't super useful but helps get the point across -- PCOA with Environment Vectors (simply maps the vectors onto an unconstrained ordination)
    plot.pcoa(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - with Vectors", env.vars=c("BMI","Years.in.US","Age"), flip.axis=2)    

    # RDA: constrain ordination with selected env variables
    # here we want to look at the % relative variation explained by the constrained model VS the unconstrained model
    plot.constrained.ordination(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR", env.vars=c("Years.in.US","BMI","Age"))    

    #plot.constrained.ordination(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR, Fraction.Life", env.vars=c("Fraction.Life.in.US"))    

    plot.constrained.ordination(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Unconstrained")    
    # RDA: with 1st gen only
    plot.constrained.ordination(map[firstgen_cs,], otu0=taxa_clr_L7[firstgen_cs,], method="euclidean", plot.title="Euclidean - CLR - 1stgen", env.vars=c("Years.in.US","BMI","Age"))    
    plot.constrained.ordination(map[firstgen_cs,], otu0=taxa_clr_L7[firstgen_cs,], method="euclidean", plot.title="Euclidean - CLR - 1stgen - Unconstrained")    


### Intra-inter group variabilities
    plot.within.group.distances(map0=map[cs,], otu0=taxa_clr_L7[cs,], fn="within.group.euc.pdf", ylab="Euclidean distance")
    plot.between.group.distances(map0=map[cs,], otu0=taxa_clr_L7[cs,], fn="between.group.euc.pdf", ylab="Euclidean distance")

### Relative Distance - Longitudinal Plots

    # by days since arrival
    # by days to first sample        
    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=1, xlab="Days Since US Arrival", 
                    ylab="Dissimilarity to First Sample", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_Day0.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=6, xlab="Days Since US Arrival", 
                    ylab="Dissimilarity to Last Sample", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_DayM6.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=1, to.previous=T, xlab="Days Since US Arrival", 
                    ylab="Dissimilarity to Previous Sample", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_Previous.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    # L to KarenThai
    l.samples <- rownames(map)[map$Sub.Study == "L" & map$Subject.ID != "IMP.000"]
    p <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=karenthai, xlab="Days Since US Arrival", 
                    ylab="Distance to Karen in Thailand", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_KarenThai.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
    # L to 2ndGenHmong
    p <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=hmong_secondgen_cs, xlab="Days Since US Arrival", 
                    ylab="Distance to 2nd-Gen", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_USborn.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)


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


### Relative Distance - Cross-Sectional Plots

    xlab="Years in US"
    x.var="Years.in.US"

    # all 1st gen to all thai
    plot.relative.CS(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_clr_L7, main="Euclidean", ref_samples=c(karenthai,hmongthai), xlab=xlab, x.var=x.var, ylab="Distance to all Thai", outputfn="allfirstgen_to_Thai.pdf")

    # 1stKaren KarenThai
    plot.relative.CS(map0=map[karen_firstgen_cs,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=karenthai, xlab=xlab, x.var=x.var, ylab="Distance to Karen Thai", outputfn="Karen_CS_to_Thai.pdf")

    # 1stKaren to 2ndGenHmong
    plot.relative.CS(map0=map[karen_firstgen_cs,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=hmong_secondgen_cs, xlab=xlab, x.var=x.var, ylab="Distance to 2nd-gen", outputfn="Karen_CS_to_USborn.pdf")

    # 1stHmong to HmongThai
    plot.relative.CS(map0=map[hmong_firstgen_cs,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=hmongthai, xlab=xlab, x.var=x.var, ylab="Distance to Hmong Thai", outputfn="Hmong_CS_to_Thai.pdf")

    # 1stHmong to 2ndHmong
    plot.relative.CS(map0=map[hmong_firstgen_cs,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=hmong_secondgen_cs, xlab=xlab, x.var=x.var, ylab="Distance to US born Hmong", outputfn="Hmong_CS_to_USbornHmong.pdf")


### Differential Taxa in Lean vs Obese
    # try all the data, but add ethnicity as a control
    plot.diff.taxa(lean.obese.cs.map, taxa_clr_L7, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="BMI.all.L7_clr", sig.level=.10, do.sqrt=FALSE)

    plot.diff.taxa(lean.obese.cs.map, taxa_clr_L6, x.var="BMI.Class", 
        control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L6_clr", sig.level=.10, do.sqrt=FALSE)

    plot.diff.taxa(lean.obese.cs.map, taxa_clr_L2, x.var="BMI.Class", 
        control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L2_clr", sig.level=.1, do.sqrt=FALSE)

    plot.diff.taxa(lean.obese.cs.map, taxa_clr_L6, x.var="Waist.Height.Ratio", 
            control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="WHR.all.L6_clr", sig.level=.10, do.sqrt=FALSE)


    # look for inflammatory bugs over time in US (within ethnicity) (for R01)    
    inflam_taxa <- cbind(taxa_clr_L3[,"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria"],
    taxa_clr_L5[,"k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae"],
    taxa_clr_L6[,c("k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Mailhella",
    "k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfovibrionales;f__Desulfovibrionaceae;g__Desulfovibrio",
    "k__Bacteria;p__Proteobacteria;c__Deltaproteobacteria;o__Desulfovibrionales")])
    colnames(inflam_taxa)[1:2] <- c("k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria", "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae")

    map_inf <-  map
Years.in.US.Factor <- map$Years.in.US
Years.in.US.Factor[Years.in.US.Factor == 0] <- 50 # 2nd gen are 0, let's set them to something higher like 50??
Years.in.US.Factor[is.na(Years.in.US.Factor)] <- 0 # Thai are NA, let's set them to 0
Years.in.US.Factor <- factor(cut(Years.in.US.Factor,c(-1, .0001, seq(5,35,5), 41, 50)), ordered=T)
map_inf[,"Years.in.US.Factor"] <- Years.in.US.Factor

    plot.diff.taxa(map_inf[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], inflam_taxa, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US.Factor"), outputfn.prepend="inflamed_hmong", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)

    plot.diff.taxa(map_inf[c(karenthai,karen_firstgen_cs),], inflam_taxa, x.var="BMI.Class", 
               control.vars=c("Age","Years.in.US.Factor"), outputfn.prepend="inflamed_karen", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)

    plot.diff.taxa(map_inf[c(hmong_firstgen_cs),], inflam_taxa, x.var="Years.in.US", 
               control.vars="Age", outputfn.prepend="inflamed_hmong1st", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)

    plot.diff.taxa(map_inf[c(karen_firstgen_cs),], inflam_taxa, x.var="Years.in.US", 
               control.vars=c("Age"), outputfn.prepend="inflamed_karen1st", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)


    # lean obese only
    plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Ethnicity=="Hmong",], inflam_taxa, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US.Factor"), outputfn.prepend="inflamed_hmong", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)

    plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Ethnicity=="Karen",], inflam_taxa, x.var="BMI.Class", 
               control.vars=c("Age","Years.in.US.Factor"), outputfn.prepend="inflamed_karen", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)

    plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Sample.Group=="Hmong1st",], inflam_taxa, x.var="Years.in.US", 
               control.vars="Age", outputfn.prepend="inflamed_hmong1st", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)

    plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Sample.Group=="Karen1st",], inflam_taxa, x.var="Years.in.US", 
               control.vars=c("Age"), outputfn.prepend="inflamed_karen1st", sig.level=1, do.sqrt=FALSE, do.filter=FALSE)


### MB americanization (relative distance to reference group) - Lean vs. Obese only
    query_samples_list <- list(c(hmong_secondgen_cs, hmong_firstgen_cs), c(hmong_secondgen_cs,hmong_firstgen_cs, karen_firstgen_cs), karen_firstgen_cs)
    ref_samples_list <- list(hmongthai, c(hmongthai, karenthai), karenthai)
    fn_ext <- c("Hmong", "All", "Karen")
    for(i in 1:length(query_samples_list))
    {
        map_temp <- map[query_samples_list[[i]],]
        map_temp <- map_temp[map_temp$BMI.Class %in% c("Lean","Obese"),]
        new_query_samples <- rownames(map_temp)
        # calculate the relative distances to reference groups as response variable
        all.samples <- c(new_query_samples, ref_samples_list[[i]])
        ddm <- vegdist(taxa_clr_L7[all.samples,], method="euc")
        dm <- as.matrix(ddm)
        rel.dist <- get.relative.distance(new_query_samples, ref_samples_list[[i]], dm)
        p <- plot.boxplot.by.group.x.bmi(map00=map_temp, y=rel.dist, ylab="Distance to Thai Reference", main="Euclidean")
        save_plot(paste0("boxplot-BMI-x-RelativeDistance-", fn_ext[i], "-CLR.pdf"), p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
    }

### M1 vs M6

    s1 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 1 & map$Subject.ID != "IMP.000",]))
    s2 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 6 & map$Subject.ID != "IMP.000",]))
    taxa0 <- taxa_clr_L6[c(s1,s2),]
    ret <- collapse.by.correlation(taxa0, .95)
    taxa0 <- taxa0[, ret$reps]

    ret <- test.features.paired(taxa0, samples1=s1, samples2=s2, sig.level=.10, parametric=TRUE)    
    #                             k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Micrococcaceae;g__Rothia 
    #                                                                                                               0.04085939 

    boxplot(c(rep(1,length(s1)), rep(2,length(s2))), taxa0[,"k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Micrococcales;f__Micrococcaceae;g__Rothia"])
    
### heatmap

    make.heatmap(taxa_clr_L6, map[cs,], .25, presence.absence=F, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.all.thaibaseline.pdf")
    
    
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.pdf")   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.karen.pdf") 

    make.heatmap(otu, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.all.thaibaseline.ordernumotus.pdf", order.by.num.otus=TRUE)
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.ordernumotus.pdf", order.by.num.otus=TRUE)   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.karen.ordernumotus.pdf", order.by.num.otus=TRUE) 


    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("HmongThai"), outputfn="heatmap.hmong.thaibaseline.pdf")   
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("Hmong2nd"), outputfn="heatmap.hmong.2ndgenbaseline.pdf")   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=c("KarenThai"), outputfn="heatmap.karen.thaibaseline.pdf") 

    make.heatmap(taxa, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.taxa.all.thaibaseline.pdf")
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("HmongThai"), outputfn="heatmap.taxa.hmong.thaibaseline.pdf")   
    make.heatmap(taxa, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("Hmong2nd"), outputfn="heatmap.taxa.hmong.2ndgenbaseline.pdf")   
    make.heatmap(taxa, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=c("KarenThai"), outputfn="heatmap.taxa.karen.thaibaseline.pdf") 

    
    
    
    
    
    
    