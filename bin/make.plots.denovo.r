datadirname <- "denovo"
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("bin/load.r")

# ----- CROSS-SECTIONAL PLOTS ----- 

### alpha diversity
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(alphadiv[cs, alpha.metrics]), ylabs=rep("",length(alpha.metrics)), mains=alpha.metrics, outputfn="boxplot.alphadiv.bmi.pdf")
    
### body trends
    plot.body.trends(map[firstgen_cs,])
    plot.WHR.boxplot(map_all[c(karenthai_all,karen_firstgen_cs_all),], bins=0:10, fn="WHR_boxplot_Karen.pdf")   
    plot.WHR.boxplot(map_all[c(hmongthai_all,hmong_firstgen_cs_all,hmong_secondgen_cs_all),], bins=seq(0,45,5), fn="WHR_boxplot_Hmong.pdf")
    
### taxa summaries
    plot.taxa.summary(map0=map[cs,], otu=taxa, fn="taxa.summary.pdf")    

### PCOA 
    plot.pcoa(map[cs,], dm=wuf_dm, plot.title="Weighted Unifrac", flip.axis=1:2)
    plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="Unweighted Unifrac")
    plot.pcoa(map[cs,], dm=bc_dm, plot.title="Bray Curtis")

    # plot PCOA of firstgen only and color and fill by different metadata vars
    maptemp <- map[c(hmong_firstgen_cs, karen_firstgen_cs),]
    maptemp$Arrived.As <- ifelse(maptemp$Age.at.Arrival >= 18, "Adult", "Child")
    p <- plot.pcoa.by(map0=maptemp, dm=bc_dm, fn=, fill.by="Years.in.US", shape.by="Arrived.As")
    save_plot(p, "pcoa-BC-Years-x-AgeArrived.Binary.pdf", base_aspect_ratio=1.3)
    
    # plot PCOA by alpha diversity
    p <- plot.pcoa.by(map0=cbind(map[cs,], alphadiv[cs, "PD_whole_tree",drop=F]), dm=bc_dm, fill.by="PD_whole_tree")
    save_plot(p, "pcoa-BC-All-PD_whole_tree.pdf", base_aspect_ratio=1.3)

    pdiv <- plot.pcoa.by(map0=cbind(maptemp, alphadiv[rownames(maptemp), "PD_whole_tree",drop=F]), dm=bc_dm, fill.by="PD_whole_tree")
    pdiv.leg <- get_legend(pdiv)
    pdiv <- pdiv + theme(legend.position='none')
    pyrs <- plot.pcoa.by(map0=maptemp, dm=bc_dm, fill.by="Years.in.US", fill.color="#00441B")
    pyrs.leg <- get_legend(pyrs)
    pyrs <- pyrs + theme(legend.position='none')
    
    multiplot <- plot_grid(pdiv, pyrs, plot_grid(pdiv.leg, pyrs.leg, ncol=1, align="v"), ncol=3, rel_widths=c(1,1,.5))
    save_plot("pcoa-BC-PD_whole_tree,Yrs.in.US.pdf", multiplot, ncol = 2, nrow = 1, base_aspect_ratio = 1.3)

    
    # RDA: constrain ordination with selected env variables using Unweighted Unifrac DM
    #plot.constrained.ordination(map[cs,], dm0=uwuf_dm, plot.title="Unweighted Unifrac", env.vars=c("Years.in.US","BMI","Age"))    
    #plot.constrained.ordination(map[cs,], dm0=uwuf_dm, plot.title="Unweighted Unifrac - Full Model")    
    
### Bacteroides-Prevotella
    plot.b.p.ratio(map[firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.pdf")
    plot.b.p.ratio(map[hmong_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.hmong.pdf")
    plot.b.p.ratio(map[karen_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.karen.pdf")

    # BP ratio over time, but colored by BMI
    plot.b.p.ratio.x.bmi(map[hmong_firstgen_cs,,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.hmong.x.BMI.pdf")
    plot.b.p.ratio.x.bmi(map[karen_firstgen_cs,,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.karen.x.BMI.pdf")
    
    # ALL samples 
    plot.b.p.ratio.all(map[cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.pdf")
    
    # BP Ratio - Lean vs. Obese
    bp <- get.log.taxa.ratio(taxa, bug1=bacteroides, bug2=prevotella) # calculate bp for everyone
    bp <- bp[cs]
    bp <- bp[is.finite(bp)]
    
    p <- plot.boxplot.by.group.x.bmi(map00=map[names(bp),], y=bp, ylab="B-P Ratio", main="") # plot cross sectional only
    save_plot("boxplot-BMI-x-BPRatio.pdf", p$p, base_aspect_ratio = 1.3)

    # prob will have NAs here, need to filter out
        p <- plot.boxplot.by.group.x.bmi(map[c(karenthai,karen_firstgen_cs),], bp[c(karenthai,karen_firstgen_cs)], "B-P Ratio", main="") # plot cross sectional only
        save_plot("boxplot-BMI-x-BPRatio-Karen.pdf", p$p, base_aspect_ratio = 1.3)
        # prob will have NAs here, need to filter out
        p <- plot.boxplot.by.group.x.bmi(map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], bp[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs)], "B-P Ratio", main="") # plot cross sectional only
        save_plot("boxplot-BMI-x-BPRatio-Hmong.pdf", p$p, base_aspect_ratio = 1.3)

    # look at BP Ratio correlations with Age
        p <- plot.XY(merge(map[hmong_secondgen_cs,], bp, by=0)[,c("Age","y")])
        save_plot("BPRatio-x-Age-Hmong2nd.pdf", p, base_aspect_ratio = 1.3)

        p <- plot.XY(merge(map[hmongthai,], bp, by=0)[,c("Age","y")])
        save_plot("BPRatio-x-Age-HmongThai.pdf", p, base_aspect_ratio = 1.3)

        p <- plot.XY(merge(map[karenthai,], bp, by=0)[,c("Age","y")])
        save_plot("BPRatio-x-Age-KarenThai.pdf", p, base_aspect_ratio = 1.3)

    # look at BP Ratio correlations with Age AND Years.in.US
        p <- plot.XY(merge(map[hmong_firstgen_cs,], bp, by=0)[,c("Age", "Years.in.US", "y")])
        save_plot("BPRatio-x-Years.Age-Hmong1st.pdf", p, base_aspect_ratio = 1.3)
            # what we see here is that the interaction is NOT significant between Age and Years
            # which means that being older + living here longer doesn't necessarily mean you'll have more BP
        p <- plot.XY(merge(map[karen_firstgen_cs,], bp, by=0)[,c("Age", "Years.in.US", "y")])
        save_plot("BPRatio-x-Years.Age-Karen1st.pdf", p, base_aspect_ratio = 1.3)

    # BP with Age.at.Arrival
        p <- plot.XY(merge(map[hmong_firstgen_cs,], bp, by=0)[,c("Age.at.Arrival", "Years.in.US", "y")])
        save_plot("BPRatio-x-Years.AgeArrival-Hmong1st.pdf", p, base_aspect_ratio = 1.3)

        p <- plot.XY(merge(map[rownames(map) %in% hmong_firstgen_cs & map$Years.in.US > 20,], bp, by=0)[,c("Age.at.Arrival", "y")])
        save_plot("BPRatio-x-AgeArrival-Hmong1st-YearsGT20.pdf", p, base_aspect_ratio = 1.3)

        p <- plot.XY(merge(map[rownames(map) %in% hmong_firstgen_cs & map$Years.in.US > 30,], bp, by=0)[,c("Age.at.Arrival", "y")])
        save_plot("BPRatio-x-AgeArrival-Hmong1st-YearsGT30.pdf", p, base_aspect_ratio = 1.3)

        p <- plot.XY(merge(map[karen_firstgen_cs,], bp, by=0)[,c("Age.at.Arrival", "Years.in.US", "y")])
        save_plot("BPRatio-x-Years.AgeArrival-Karen1st.pdf", p, base_aspect_ratio = 1.3)
        p <- plot.XY(merge(map[karen_firstgen_cs,], bp, by=0)[,c("Age.at.Arrival", "y")])
        save_plot("BPRatio-x-AgeArrival-Karen1st.pdf", p, base_aspect_ratio = 1.3)

        p <- plot.XY(merge(map[c(karen_firstgen_cs,hmong_firstgen_cs),], bp, by=0)[,c("Age.at.Arrival", "Years.in.US", "y")])
        save_plot("BPRatio-x-Years.AgeArrival-HmongKaren1st.pdf", p, base_aspect_ratio = 1.3)

### Intra-inter group variabilities - consider doing in CLR only
    bc_ret <- plot.within.group.distances(map0=map[cs,], dm=bc_dm, fn="within.group.bc.pdf", ylab="Bray-Curtis distance")
    wuf_ret <- plot.within.group.distances(map0=map[cs,], dm=wuf_dm, fn="within.group.wuf.pdf", ylab="Weighted Unifrac distance")
    uwuf_ret <- plot.within.group.distances(map0=map[cs,], dm=uwuf_dm, fn="within.group.uwuf.pdf", ylab="Unweighted Unifrac distance")

    plot.between.group.distances(map0=map[cs,], dm=bc_dm, fn="between.group.bc.pdf", ylab="Bray-Curtis distance")
    plot.between.group.distances(map0=map[cs,], dm=wuf_dm, fn="between.group.wuf.pdf", ylab="Weighted Unifrac distance")

### Relative Distance - consider doing in CLR only
    mains <- c("Unweighted Unifrac","Weighted Unifrac", "Bray-Curtis")
    dms <- list(uwuf_dm,wuf_dm,bc_dm)
    xlab="Years in US"
    x.var="Years.in.US"

    # all 1st gen to all thai
    multiplot.relative.CS(mains, dms, map[c(hmong_firstgen_cs, karen_firstgen_cs),], ylab="Distance to all Thai", xlab="BMI", x.var="BMI", 
                        ref_samples=c(karenthai,hmongthai), outputfn="BMI - allfirstgen_to_Thai.pdf")

    # all 1st gen to all thai
    multiplot.relative.CS(mains, dms, map[c(hmong_firstgen_cs, karen_firstgen_cs),], ylab="Distance to all Thai", xlab=xlab, x.var=x.var, 
                        ref_samples=c(karenthai,hmongthai), outputfn="allfirstgen_to_Thai.pdf")
    # 1stKaren KarenThai
    multiplot.relative.CS(mains, dms, map[karen_firstgen_cs,], ylab="Distance to Karen Thai", xlab=xlab, x.var=x.var,
                        ref_samples=karenthai, outputfn="Karen_CS_to_Thai.pdf")
    # 1stKaren to 2ndGenHmong
    multiplot.relative.CS(mains, dms, map[karen_firstgen_cs,], ylab="Distance to 2nd-Generation Hmong", xlab=xlab, x.var=x.var, 
                        ref_samples=hmong_secondgen_cs, outputfn="Karen_CS_to_USbornHmong.pdf")
    # 1stHmong to HmongThai
    multiplot.relative.CS(mains, dms, map[hmong_firstgen_cs,], ylab="Distance to Hmong Thai", xlab=xlab, x.var=x.var,
                        ref_samples=hmongthai, outputfn="Hmong_CS_to_Thai.pdf")
    # 1stHmong to 2ndHmong
    multiplot.relative.CS(mains, dms, map[hmong_firstgen_cs,], ylab="Distance to US-born Hmong", xlab=xlab, x.var=x.var,
                        ref_samples=hmong_secondgen_cs, outputfn="Hmong_CS_to_USbornHmong.pdf")

    # all Hmong to Caucasians
    multiplot.relative.CS(mains, dms, map[c(hmong_firstgen_cs, hmong_secondgen_cs),], ylab="Distance to Caucasians", xlab=xlab, x.var=x.var,
                        ref_samples=controls, outputfn="Hmong_CS_to_Controls.pdf")

    # Lean vs. Obese only
    # Hypothesis: are obese individuals more different from their groups?
    # 1. Hmong1st and Hmong2nd to HmongThai    
    # 2. all US to all Thai
    # 3. Karen1st to KarenThai
    query_samples_list <- list(c(hmong_secondgen_cs, hmong_firstgen_cs), c(hmong_secondgen_cs,hmong_firstgen_cs, karen_firstgen_cs), karen_firstgen_cs, c(hmong_firstgen_cs, hmong_secondgen_cs, hmongthai))
    ref_samples_list <- list(hmongthai, c(hmongthai, karenthai), karenthai, controls)
    fn_ext <- c("Hmong", "All", "Karen", "Controls")
    for(i in 1:length(query_samples_list))
    {
        map_temp <- map[query_samples_list[[i]],]
        map_temp <- map_temp[map_temp$BMI.Class %in% c("Lean","Obese"),]
        new_query_samples <- rownames(map_temp)
        # calculate the relative distances to reference groups as response variable
        rel.dists <- lapply(dms, function(xx) get.relative.distance(new_query_samples, ref_samples_list[[i]], xx))    
        multiplot.boxplot.by.group.x.bmi(map_temp, rel.dists, rep("Distance to Reference",3), mains, outputfn = paste0("boxplot-BMI-x-RelativeDistance-",fn_ext[i],".pdf"))
    }

### Heatmap
    # see CLR for traditional heatmap
    make.heatmap.binary(otu0=taxa, map0=map[c(hmongthai,hmong_firstgen_cs,hmong_secondgen_cs),], min.prevalence=0.5, baseline.groups="Hmong2nd", show.colnames=F, outputfn="heatmap.taxa.hmongthai.binary.pdf")
    make.heatmap.binary(otu0=taxa, map0=map[c(hmongthai,hmong_firstgen_cs),], min.prevalence=0.5, baseline.groups="HmongThai", show.colnames=F, outputfn="heatmapbinary.taxa.hmong1st.thai.pdf")
    make.heatmap.binary(otu0=taxa, map0=map[cs,], min.prevalence=0.5, baseline.groups=c("KarenThai","HmongThai"), show.colnames=F, outputfn="heatmapbinary.taxa.all.thai.pdf")

    make.heatmap.binary(otu0=taxa, map0=map[c(hmongthai,hmong_firstgen_cs,hmong_secondgen_cs),], min.prevalence=0.5, baseline.groups="Hmong2nd", show.colnames=F, outputfn="heatmap.taxa.hmongthai.binary.pdf")
    
# ----- LONGITUDINAL PLOTS ----- 
    
    map_L <- map[map$Sub.Study == "L",]
    # IMP.049 and IMP.050 were recruited at 2 months and 3 months respectively - these might be throwing things off
    map_L <- map_L[!(map_L$Subject.ID %in%  c("IMP.049", "IMP.050")),]

### alpha diversity, one per subject
    multiplot.alphadiv.L(map_L, alphadiv, alpha.metrics, outputfn="alphadiv.L.Day.pdf") # detaults to by day since arrival
    multiplot.alphadiv.L(map_L, alphadiv, alpha.metrics, outputfn="alphadiv.L.Num.pdf", x.var="Sample.Order", xlab="Sample Number") # by sample number

### taxa summaries
    # stream plots are useful for longitudinal changes over time
    plot.taxa.summary.L(taxa0=taxa, map0=map[map$Sub.Study=="L", ], outputfn="taxa.summary.L.pdf", grid.ncol=2)

### PCOA - consider using CLR
    plot.pcoa.long(map, samples=rownames(map_L), dm=uwuf_dm, plot.title="UnWeighted Unifrac - L")
    plot.pcoa.long(map, samples=rownames(map_L), dm=wuf_dm, plot.title="Weighted Unifrac - L")

    plot.pcoa.long(map, samples=rownames(map_L), dm=uwuf_dm, plot.title="UnWeighted Unifrac - L.outlined", convex.hull=TRUE)
    plot.pcoa.long(map, samples=rownames(map_L), dm=wuf_dm, plot.title="Weighted Unifrac - L.outlined", convex.hull=TRUE)
    plot.pcoa.long(map, samples=rownames(map_L), dm=bc_dm, plot.title="Bray-Curtis - L.outlined", convex.hull=TRUE)

### BP Ratio
    # all months
    plot.b.p.ratio.L(map_L, taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.L.pdf")
    # first available month to last available month
    plot.b.p.ratio.L(map_L, taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.L.first.last.pdf", num.clip.months=1)
    # first 2 available months to last 2 available months
    plot.b.p.ratio.L(map_L, taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.L.first2.last2.pdf", num.clip.months=2)
    
    
### Relative Distance - consider using CLR

    # by days since arrival
    xlab<-"Days Since US Arrival"
    
    # by days to first sample        
    multiplot.relative.L(mains=mains, dms=dms, map0 = map_L, ylab="Dissimilarity to First Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Karen_L_to_Day0.pdf")
    
    # by days to last sample    
    multiplot.relative.L(mains=mains, dms=dms, map0 = map_L, ylab="Dissimilarity to Last Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=6, outputfn="Karen_L_to_Day_M6.pdf")

    # days to first sample - day-to-day
    multiplot.relative.L(mains=mains, dms=dms, map0 = map_L, ylab="Dissimilarity to Previous Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Karen_L_to_Previous.pdf", to.previous=T)
    
    # L to KarenThai
    multiplot.relative.L(mains=mains, dms=dms, map0 = map_L, ylab="Distance to Karen in Thailand", 
                        xlab=xlab, x.var="Sample.Day.Since.Arrival", ref_samples=karenthai, outputfn="Karen_L_to_KarenThai.pdf")
    # L to 2ndGenHmong
    multiplot.relative.L(mains=mains, dms=dms, map0 = map_L, ylab="Distance to 2nd-Generation", 
                        xlab=xlab, x.var="Sample.Day.Since.Arrival", ref_samples=hmong_secondgen_cs, outputfn="Karen_L_to_USborn.pdf")

### DIFF TAXA --> do this using only CLR transformed data and parametric tests
