setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
datadirname <- "datadir_denovo"
source("bin/load.r")

taxa_clr_L7_fn <- paste(datadir,"taxa.clr.L7.txt",sep="/")
taxa_clr_L6_fn <- paste(datadir,"taxa.clr.L6.txt",sep="/")
taxa_clr_L2_fn <- paste(datadir,"taxa.clr.L2.txt",sep="/")
taxa_clr_L7 <- load.data(mapfile, otufile=taxa_clr_L7_fn, normalize=F)$otu
taxa_clr_L6 <- load.data(mapfile, otufile=taxa_clr_L6_fn, normalize=F)$otu
taxa_clr_L2 <- load.data(mapfile, otufile=taxa_clr_L2_fn, normalize=F)$otu


### PCOA

    # remove IMP.000 from any dm calcs
    plot.pcoa.long(map[map$Subject.ID != "IMP.000",], samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), otu0 = taxa_clr_L7, method="euclidean", plot.title="Euclidean - CLR - L.outlined", convex.hull=TRUE)

    # let's flip the axis so that Years in US goes left to right, BMI/Age goes bottom to up
    plot.pcoa(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR", show.stats=F)    

    plot.pcoa(map, otu0=taxa_clr_L7[rownames(map),], method="euclidean", plot.title="Euclidean - CLR", show.stats=F)    


    # this isn't super useful but helps get the point across -- PCOA with Environment Vectors (simply maps the vectors onto an unconstrained ordination)
    plot.pcoa(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - PCOA with Vectors", env.vars=c("BMI","Years.in.US","Age"), flip.axis=2)    

    # RDA: constrain ordination with selected env variables
    # here we want to look at the % relative variation explained by the constrained model VS the unconstrained model
    plot.constrained.ordination(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Constrained", env.vars=c("Years.in.US","BMI","Age"))    

    #plot.constrained.ordination(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR, Fraction.Life", env.vars=c("Fraction.Life.in.US"))    

    plot.constrained.ordination(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Unconstrained")    
    # RDA: with 1st gen only
    plot.constrained.ordination(map[firstgen_cs,], otu0=taxa_clr_L7[firstgen_cs,], method="euclidean", plot.title="Euclidean - CLR - 1stgen", env.vars=c("Years.in.US","BMI","Age"))    
    plot.constrained.ordination(map[firstgen_cs,], otu0=taxa_clr_L7[firstgen_cs,], method="euclidean", plot.title="Euclidean - CLR - 1stgen - Unconstrained")    


    # try plotting constrained ordination with 1st 5 diet PCs as included vars
    
    # generate food PCs here (weighted unifrac)
    plot.pcoa(map_all[cs_all,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", save.pc=TRUE, axis1 = 1, axis2 = 5) # put axis 5 so we get all 5 PCs back
    # load in food PCs
    food.pc <- read.table("Food Weighted Unifrac-PC.txt", sep="\t", skip=1, row=1)
    # note that IMP.000 was sequenced, but we don't want to include it for food analyses
    valid_cs <- cs[cs != "IMP000.L.T30"]
    colnames(food.pc) <- paste0("Food.PC", 1:5)
    # redo constrained ordination with just metadata vars of interest
    plot.constrained.ordination(data.frame(map[valid_cs,], food.pc[valid_cs,]), otu0=taxa_clr_L7[valid_cs,], method="euclidean", plot.title="Euclidean - CLR - Constrained Vars", env.vars=c("Years.in.US","BMI","Age", "Ethnicity"))    
    # redo constrained ordination with metadata vars AND food PCs
    plot.constrained.ordination(data.frame(map[valid_cs,], food.pc[valid_cs,]), otu0=taxa_clr_L7[valid_cs,], method="euclidean", plot.title="Euclidean - CLR - Constrained Vars and Food PCs", env.vars=c("Years.in.US","BMI", "Ethnicity","Age",colnames(food.pc)))    
    # redo constrained ordination with food PCs only
    plot.constrained.ordination(data.frame(map[valid_cs,], food.pc[valid_cs,]), otu0=taxa_clr_L7[valid_cs,], method="euclidean", plot.title="Euclidean - CLR - Food PCs", env.vars=colnames(food.pc))    
    
    
### Intra-inter group variabilities
    plot.within.group.distances(map0=map[cs,], otu0=taxa_clr_L7[cs,], fn="within.group.euc.pdf", ylab="Euclidean distance")
    plot.between.group.distances(map0=map[cs,], otu0=taxa_clr_L7[cs,], fn="between.group.euc.pdf", ylab="Euclidean distance")

### Relative Distance - Longitudinal Plots

    # by days since arrival
    # by days to first sample        
    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=1, xlab="Days Since US Arrival", 
                    ylab="Dissimilarity to First Sample", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_Day0.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=1, xlab="Sample Number", 
                    ylab="Dissimilarity to First Sample", x.var="Sample.Order")
    save_plot("Karen_L_to_Day0_month.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=6, xlab="Days Since US Arrival", 
                    ylab="Dissimilarity to Last Sample", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_DayM6.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=6, xlab="Sample Number", 
                    ylab="Dissimilarity to Last Sample", x.var="Sample.Order")
    save_plot("Karen_L_to_DayM6_month.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)


    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=1, to.previous=T, xlab="Days Since US Arrival", 
                    ylab="Dissimilarity to Previous Sample", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_Previous.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
    p <- plot.relative.L(map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], otu0=taxa_clr_L7, main="Euclidean", ref.sample.order=1, to.previous=T, xlab="Sample Number", 
                    ylab="Dissimilarity to Previous Sample", x.var="Sample.Order")
    save_plot("Karen_L_to_Previous_month.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    # L to KarenThai
    l.samples <- rownames(map)[map$Sub.Study == "L" & map$Subject.ID != "IMP.000"]
    p <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=karenthai, xlab="Days Since US Arrival", 
                    ylab="Distance to Karen in Thailand", x.var="Sample.Day.Since.Arrival")
    save_plot("Karen_L_to_KarenThai.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    p <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=karenthai, xlab="Sample Number", 
                    ylab="Distance to Karen in Thailand", x.var="Sample.Order")
    save_plot("Karen_L_to_KarenThai_month.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)

    # L to Controls
    p <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=controls, xlab="Sample Number", 
                    ylab="Distance to Controls", x.var="Sample.Order")
    save_plot("Karen_L_to_Controls_month.pdf", p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)


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

    # 1stHmong to Controls
    plot.relative.CS(map0=map[hmong_firstgen_cs,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=controls, xlab=xlab, x.var=x.var, ylab="Distance to Controls", outputfn="Hmong_CS_to_Controls.pdf")

    plot.relative.CS(map0=map[karen_firstgen_cs,], otu0=taxa_clr_L7, main="Euclidean", ref_samples=controls, xlab=xlab, x.var=x.var, ylab="Distance to Controls", outputfn="Karen_CS_to_Controls.pdf")

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
        save_plot(paste0("boxplot-BMI-x-RelativeDistance-", fn_ext[i], "-CLR.pdf"), p$p, ncol = 1, nrow = 1, base_aspect_ratio = 1.3)
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

    make.heatmap(taxa_clr_L7, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.all.thaibaseline.pdf")
    
    
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.pdf")   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.karen.pdf") 

    make.heatmap(otu, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.all.thaibaseline.ordernumotus.pdf", order.by.num.otus=TRUE)
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.ordernumotus.pdf", order.by.num.otus=TRUE)   
#    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs, controls),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.ordernumotus.pdf", order.by.num.otus=TRUE)   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.karen.ordernumotus.pdf", order.by.num.otus=TRUE) 


    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("HmongThai"), outputfn="heatmap.hmong.thaibaseline.pdf")   
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("Hmong2nd"), outputfn="heatmap.hmong.2ndgenbaseline.pdf")   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=c("KarenThai"), outputfn="heatmap.karen.thaibaseline.pdf") 

    make.heatmap(taxa, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.taxa.all.thaibaseline.pdf")
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("HmongThai"), outputfn="heatmap.taxa.hmong.thaibaseline.pdf")   
    make.heatmap(taxa, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("Hmong2nd"), outputfn="heatmap.taxa.hmong.2ndgenbaseline.pdf")   
    make.heatmap(taxa, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=c("KarenThai"), outputfn="heatmap.taxa.karen.thaibaseline.pdf") 


    ret <- load.data(mapfile, otufile=taxa_L7_fn, normalize=F)
    taxa_L7 <- ret$otu
    bpnames <- c(grep("k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides.*",colnames(taxa_L7), value=T),
                grep("k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella.*",colnames(taxa_L7), value=T))
          
    make.heatmap(taxa_L7[,bpnames], map, .25, presence.absence=T, baseline.groups="HmongThai", outputfn="heatmap.h.BP.pdf") 
    make.heatmap(taxa_L7[,bpnames], map, .25, presence.absence=T, baseline.groups="Control", outputfn="heatmap.c.BP.pdf") 


### Predict Years in US using different covariates 
#   1) Covariates only (Age, BMI, Waist, Breastfed, Ethnicity, Medical.Assistance, Public.Housing, Children.Free.Lunch, Highest.Education, Religion)
#   2) Covariates + MB    
#   3) Diet? Diet + MB?    
    
    mapcs <- map[cs,]
    
    # remove "IMP000" since food doesn't exist for it
    mapcs <- mapcs[mapcs$Subject.ID != "IMP.000",]
    
    vars <- c("Age", "BMI", "Waist.Height.Ratio", "Breastfed", "Ethnicity", "Medical.Assistance", "Public.Housing", "Children.Free.Lunch", "Highest.Education", "Religion")

    mapcs[,vars[-(1:3)]] <- apply(mapcs[,vars[-(1:3)]], 1:2, as.character)
    varset1 <- mapcs[,vars]
    varset1[is.na(varset1)] <- "NA"
    templist <- lapply(varset1[,-(1:3)], as.factor)
    
    # for variable sets here
    vars_only <- cbind(varset1[,1:3], do.call(cbind.data.frame, templist))

    samples <- rownames(mapcs)
    
    mb_only <- as.data.frame(taxa_clr_L7[samples,])
    
    vars_mb <- cbind(vars_only, mb_only)
    
    food_only <- food_otu_L3[samples,]
    
    # this takes extremely long -- save data from here and run this line on teraminx instead
    # preds <- compare.rf.vars(mapcs$Years.in.US, list(vars_only, mb_only, vars_mb, food_only), n=10)
    
    preds <- list()
    varlist <- list(vars_only, mb_only, vars_mb, food_only)
    for(i in 1:length(varlist))
    {
        preds[[i]] <- compare.rf.vars(mapcs$Years.in.US, varsdf=varlist[[i]], n=1)
    }
    
    preds.df <- do.call(cbind.data.frame, preds)
    colnames(preds.df) <- c("Vars.Only", "MB.Only", "Vars.MB", "Food.Only")
    preds.df <- cbind(Sample.ID = samples, Actual = mapcs$Years.in.US, preds.df)

    ggdata <- melt(preds.df, id.vars=c("Actual","Sample.ID"))
    
    pre <- rownames(mapcs[mapcs$Years.in.US == 0,])
    secondgen <- rownames(mapcs[mapcs$Years.in.US == 50,])
    firstgen <- rownames(mapcs[!(mapcs$Years.in.US %in% c(0,50)),])
    
    cols <- c("#5fa55a", "#01b4bc", "#fa8925", "#fa5457")
    cols2 <- alpha(c("#5fa55a", "#01b4bc", "#fa8925", "#fa5457"), .5)
    names(cols) <- c("Vars.Only", "MB.Only", "Vars.MB", "Food.Only")
    names(cols2) <- c("Vars.Only", "MB.Only", "Vars.MB", "Food.Only")
        
    p.first <- ggplot(ggdata[ggdata$Sample.ID %in% firstgen,], aes(x = Actual, y = value, group=variable, color = variable)) +  
           geom_point() + geom_line() + scale_color_manual(values=cols2) + theme(legend.position='none') + xlab("Actual") + ylab("Predicted") + ggtitle("1st-Gen")
               
    p.pre <- ggplot(ggdata[ggdata$Sample.ID %in% pre,], aes(x = variable, y = value, group=variable, color = variable)) +  
            geom_boxplot(aes(fill=variable), colour="black", width=.9) + 
            scale_fill_manual(values = cols) + 
            ggtitle("Pre") + ylab("Predicted") + xlab("") + 
            theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), strip.background = element_blank(),
            strip.text = element_blank()) + coord_flip() + geom_hline(yintercept = 0, color="red")  + theme(legend.position='none')


    p.second <- ggplot(ggdata[ggdata$Sample.ID %in% secondgen,], aes(x = variable, y = value, group=variable, color = variable)) +  
            geom_boxplot(aes(fill=variable), colour="black", width=.9) + 
            scale_fill_manual(values = cols) +
            ggtitle("2nd-Gen") + ylab("Predicted") + xlab("") + 
            theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), strip.background = element_blank(),
            strip.text = element_blank()) + coord_flip() + geom_hline(yintercept = 50, color="red")

    leg <- get_legend(p.second)
    p.second <- p.second + theme(legend.position='none')
    
    p.combined <- plot_grid(plot_grid(p.pre, p.first, p.second, nrow=1, rel_widths=c(1,3,1)), leg, nrow=1, rel_widths=c(4,1))
    
    save_plot("rf.vars.pdf", p.combined, base_aspect_ratio = 1.7)
    
### DIFF TAXA --> SKIP, only do diff taxa tests using CLR transformed data and parametric tests
    # between M1 vs M6
#     s1 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 1,]))
#     s2 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 6,]))
#     taxa0 <- taxa[c(s1,s2),]
#     prevalences <- apply(taxa0, 2, function(bug.col) mean(bug.col > 0))
#     taxa0 <- taxa0[, prevalences >= .10]
#     ret <- collapse.by.correlation(taxa0, .95)
#     taxa0 <- taxa0[, ret$reps]
# 
#     ret <- test.features.paired(taxa0, samples1=s1, samples2=s2, sig.level=.10, parametric=FALSE)
    ### nothing significant
                   
### CS Differential Taxa --> SKIP, only do diff taxa tests using CLR transformed data and parametric tests
        # by BMI classes
        #         ret <- load.data(mapfile, otufile=taxa_L2_fn, normalize=F)
        #         taxa_L2 <- ret$otu
        #         ret <- load.data(mapfile, otufile=taxa_L6_fn, normalize=F)
        #         taxa_L6 <- ret$otu
        # 
        #         #iterate through all combos of looking at diff taxa in lean vs. obese by subgroup
        #         groups <- c("HmongThai","Hmong1st","Hmong2nd","KarenThai","Karen1st")
        #         taxatables <- list(taxa, taxa_L2, taxa_L6)
        #         names(taxatables) <- c("L6", "L2", "L6_rare")
        #         combos <- expand.grid(names(taxatables),groups, stringsAsFactors=F)
        #         colnames(combos) <- c("taxa","group")
        #         combos[combos$group %in% c("KarenThai","HmongThai", "Hmong2nd"), "controls"] <- "Age"
        #         combos[combos$group %in% c("Karen1st","Hmong1st"), "controls"] <- c("Age,Years.in.US")
        # 
        #         sig.level=.25    
        #         lean.obese.cs.map <- map[map$BMI.Class %in% c("Lean", "Obese") & (is.na(map$Sample.Order) | map$Sample.Order==1),]
        # 
        #         for(i in 1:nrow(combos))
        #         {
        #             combo.name <- paste("BMI", paste(unlist(combos[i,c("group","taxa")]), collapse="."), sep=".")
        #             print(combo.name) # so we know which iteration we are in
        #             plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Sample.Group == combos$group[i],], taxatables[[combos$taxa[i]]], x.var="BMI.Class", 
        #                 control.vars=unlist(strsplit(combos$controls[i],",")), outputfn.prepend=combo.name, sig.level=sig.level)
        #         }
        # 
        #         # create factor version of years in us so that data isn't omitted for NAs or 0s - to be used only when dealing with the entire dataset
        #         Years.in.US.Factor <- lean.obese.cs.map$Years.in.US
        #         Years.in.US.Factor[Years.in.US.Factor == 0] <- 50 # 2nd gen are 0, let's set them to something higher like 50??
        #         Years.in.US.Factor[is.na(Years.in.US.Factor)] <- 0 # Thai are NA, let's set them to 0
        # 
        #         Years.in.US.Factor <- factor(cut(Years.in.US.Factor,c(-1, .0001, seq(5,35,5), 41, 50)), ordered=T)
        # 
        #         lean.obese.cs.map[,"Years.in.US.Factor"] <- Years.in.US.Factor
        # 
        #         # try all the data, but add ethnicity as a control
        #         plot.diff.taxa(lean.obese.cs.map, taxa, x.var="BMI.Class", 
        #                 control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="BMI.all", sig.level=.10)
        # 
        #         plot.diff.taxa(lean.obese.cs.map, taxa_L6, x.var="BMI.Class", 
        #             control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L6", sig.level=.10)
        # 
        #         plot.diff.taxa(lean.obese.cs.map, taxa_L2, x.var="BMI.Class", 
        #             control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L2", sig.level=1)
        # 
        # by WHR
        #         plot.diff.taxa(lean.obese.cs.map, taxa, x.var="Waist.Height.Ratio", control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="WHR.all", sig.level=.10)
