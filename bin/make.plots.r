setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
use.rare <- FALSE

# load all source files, data, and prep any variables for reuse
source("bin/load.r")

# load some extra vars here
# TODO: alpha div should probably just be done with rarefied tables to get rid of biases towards higher depth samples
    
    alpha.metrics <- c("observed_otus", "shannon", "simpson") #c("PD_whole_tree", "shannon", "simpson")
    
    # let's make an extra var for IMP000
    map.000 <- map[map$Subject.ID == "IMP.000",]
    map.000$Sample.Day.Since.First.Sample <- as.numeric(as.Date(map.000$Sample.Date, format="%m/%d/%y") - as.Date("08/01/16", format="%m/%d/%y")) # hard code first sample date
    map.000[,"travel.phase"] <- "Traveling"
    map.000[map.000$Sample.Day.Since.First.Sample <= 2,"travel.phase"] <- "Pre"
    map.000[map.000$Sample.Day.Since.First.Sample >= 28,"travel.phase"] <- "Post"
    map.000$travel.phase <- factor(map.000$travel.phase, levels=c("Pre", "Traveling", "Post")) 


### alpha diversity
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(alphadiv0[cs, alpha.metrics]), ylabs=rep("",3), mains=alpha.metrics, outputfn="boxplot.alphadiv.bmi.pdf")
    
    # plot alpha for longitudinal subjects only
    multiplot.alphadiv.L(map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], alphadiv0, alpha.metrics, "alphadiv.L.pdf")
    multiplot.alphadiv.L(map.000, alphadiv0, alpha.metrics, "alphadiv.IMP000.pdf", x.var="Sample.Day.Since.First.Sample")
    
### body trends
    plot.body.trends(map[firstgen_cs,])

    plot.WHR.boxplot(map_all[c(karenthai_all,karen_firstgen_cs_all),], bins=0:10, fn="WHR_boxplot_Karen.pdf")   
    plot.WHR.boxplot(map_all[c(hmongthai_all,hmong_firstgen_cs_all,hmong_secondgen_cs_all),], bins=seq(0,45,5), fn="WHR_boxplot_Hmong.pdf")
    
    # plot everyone, not just those that have been sequenced
    plot.BMI.barplot(map_all[c(karenthai_all,karen_firstgen_cs_all),], bins=0:10, fn="BMI_barplot_Karen.pdf", freq=F)   
    plot.BMI.barplot(map_all[c(hmongthai_all,hmong_firstgen_cs_all,hmong_secondgen_cs_all),], bins=seq(0,45,5), fn="BMI_barplot_Hmong.pdf", freq=F)
    plot.BMI.barplot(map_all[c(karenthai_all,karen_firstgen_cs_all),], bins=0:10, fn="BMI_barplot_Karen_abs.pdf", freq=T)   
    plot.BMI.barplot(map_all[c(hmongthai_all,hmong_firstgen_cs_all,hmong_secondgen_cs_all),], bins=seq(0,45,5), fn="BMI_barplot_Hmong_abs.pdf", freq=T)
    
    plot.BMI.barplot(map_all[cs_all,], bins=seq(0,45,5), fn="BMI_barplot_ALL.pdf", freq=F)   

### taxa summaries
    # stream plots are useful for longitudinal changes over time
    plot.taxa.summary.L(taxa0=taxa, map0=map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000", ], outputfn="taxa.summary.L.pdf", grid.ncol=2)
    plot.taxa.summary.L(taxa0=taxa, map0=map.000, x.var="Sample.Day.Since.First.Sample", outputfn="taxa.summary.IMP000.pdf")

    plot.taxa.summary(map0=map[cs,], otu=taxa, fn="taxa.summary.pdf")    

    # plot self samples
    ssamples <- map[map$Subject.ID == "IMP.000", "Sample.Order", drop=F]
    ssamples.ordered <- rownames(ssamples[ order(ssamples[,1]), , drop=F]) 
    day.labels <- as.character(sort(ssamples[,1])-3)
    day.labels[(length(day.labels)-3):length(day.labels)] <- paste0("P",1:4)
    plot.taxa.summary(map0=map[ssamples.ordered,], otu=taxa, fn="taxa.summary.IMP000.pdf", sample.order=ssamples.ordered,
                            x.labels = day.labels)
    
### PCOA

    plot.pcoa.long(map, samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), dm=uwuf_dm, plot.title="UnWeighted Unifrac - L")
    plot.pcoa.long(map, samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), dm=wuf_dm, plot.title="Weighted Unifrac - L")

    plot.pcoa.long(map, samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), dm=uwuf_dm, plot.title="UnWeighted Unifrac - L.outlined", convex.hull=TRUE)
    plot.pcoa.long(map, samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), dm=wuf_dm, plot.title="Weighted Unifrac - L.outlined", convex.hull=TRUE)
    plot.pcoa.long(map, samples=rownames(map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000",]), dm=bc_dm, plot.title="Bray-Curtis - L.outlined", convex.hull=TRUE)


    # plot IMP.000 by pre, traveling, post phases 
    #plot.pcoa.long(map, rownames(map[map$Subject.ID == "IMP.000",]), bc_dm, "BC - IMP.000")

    plot.pcoa(map[cs,], dm=wuf_dm, plot.title="Weighted Unifrac")
    plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="Unweighted Unifrac")
    plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="Unweighted Unifrac", axis1=1, axis2=3) # try pc1 and pc3
    plot.pcoa(map[cs,], dm=bc_dm, plot.title="Bray Curtis")

    # RDA: constrain ordination with selected env variables using Unweighted Unifrac DM
    plot.constrained.ordination(map[cs,], dm0=uwuf_dm, plot.title="Unweighted Unifrac", env.vars=c("Years.in.US","BMI","Age"))    
    plot.constrained.ordination(map[cs,], dm0=uwuf_dm, plot.title="Unweighted Unifrac - Full Model")    
    
    #plot.pcoa.by(map0=map[cs,], wuf_dm, fn="KCK.fmt.pcoa.pdf", ethnicity="Hmong", color.by="Years.in.US", and.by="BMI.Class", label.samples=paste0("TFSCS0", 20:29))

### Bacteroides-Prevotella
    plot.b.p.barplot(map[c(karenthai,karen_firstgen_cs),], taxa, bins=0:10, fn="b.p.barplot.karen.pdf")
    plot.b.p.barplot(map[cs,], taxa, bins=seq(0,45,5), fn="b.p.barplot.pdf")

    # CS samples only
    plot.b.p.ratio(map[firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.pdf")
    plot.b.p.ratio(map[hmong_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.hmong.pdf")
    plot.b.p.ratio(map[karen_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.karen.pdf")

    # BP ratio over time, but colored by BMI
    plot.b.p.ratio.x.bmi(map[hmong_firstgen_cs,,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.hmong.x.BMI.pdf")
    plot.b.p.ratio.x.bmi(map[karen_firstgen_cs,,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.karen.x.BMI.pdf")
    
    # ALL samples 
    plot.b.p.ratio.all(map, taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.pdf", g1=firstgen_cs, g2=c(hmongthai,karenthai),g3=(hmong_secondgen_cs))

    # plot longitudinal samples only
    plot.b.p.ratio(map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.L.pdf", longitudinal=T)

### Prediction 
    prediction <- predict.years(map[map$Ethnicity=="Karen" & (is.na(map$Sample.Order) | map$Sample.Order==1),], taxa)

### Intra-inter group variabilities
    bc_ret <- plot.within.group.distances(map0=map[cs,], dm=bc_dm, fn="within.group.bc.pdf", ylab="Bray-Curtis distance")
    wuf_ret <- plot.within.group.distances(map0=map[cs,], dm=wuf_dm, fn="within.group.wuf.pdf", ylab="Weighted Unifrac distance")
    uwuf_ret <- plot.within.group.distances(map0=map[cs,], dm=uwuf_dm, fn="within.group.uwuf.pdf", ylab="Unweighted Unifrac distance")

    plot.between.group.distances(map0=map[cs,], dm=bc_dm, fn="between.group.bc.pdf", ylab="Bray-Curtis distance")
    plot.between.group.distances(map0=map[cs,], dm=wuf_dm, fn="between.group.wuf.pdf", ylab="Weighted Unifrac distance")

### Relative Distance - Longitudinal Plots

    # by days since arrival
    mains <- c("Unweighted Unifrac","Weighted Unifrac", "Bray-Curtis")
    dms <- list(uwuf_dm,wuf_dm,bc_dm)
    xlab<-"Days Since US Arrival"
    # by days to first sample        
    multiplot.relative.L(mains=mains, dms=dms, map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Dissimilarity to First Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Karen_L_to_Day0.pdf")
    
    # by days to last sample    
    multiplot.relative.L(mains=mains, dms=dms, map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Dissimilarity to Last Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=6, outputfn="Karen_L_to_Day_M6.pdf")

    # days to first sample - day-to-day
    multiplot.relative.L(mains=mains, dms=dms, map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Dissimilarity to Previous Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Karen_L_to_Previous.pdf", to.previous=T)
    
    # L to KarenThai
    multiplot.relative.L(mains=mains, dms=dms, map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Distance to Karen in Thailand", 
                        xlab=xlab, x.var="Sample.Day.Since.Arrival", ref_samples=karenthai, outputfn="Karen_L_to_KarenThai.pdf")
    # L to 2ndGenHmong
    multiplot.relative.L(mains=mains, dms=dms, map0 = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Distance to 2nd-Generation", 
                        xlab=xlab, x.var="Sample.Day.Since.Arrival", ref_samples=hmong_secondgen_cs, outputfn="Karen_L_to_USborn.pdf")

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

### Relative Distance - Cross-Sectional Plots

    xlab="Years in US"
    x.var="Years.in.US"

    # all 1st gen to all thai
    multiplot.relative.CS(mains, dms, map, ylab="Distance to all Thai", xlab="BMI", x.var="BMI", query_samples=c(hmong_firstgen_cs, karen_firstgen_cs), 
                        ref_samples=c(karenthai,hmongthai), outputfn="BMI - allfirstgen_to_Thai.pdf")

    # all 1st gen to all thai
    multiplot.relative.CS(mains, dms, map, ylab="Distance to all Thai", xlab=xlab, x.var=x.var, query_samples=c(hmong_firstgen_cs, karen_firstgen_cs), 
                        ref_samples=c(karenthai,hmongthai), outputfn="allfirstgen_to_Thai.pdf")
    # 1stKaren KarenThai
    multiplot.relative.CS(mains, dms, map, ylab="Distance to Karen Thai", xlab=xlab, x.var=x.var, query_samples=karen_firstgen_cs, 
                        ref_samples=karenthai, outputfn="Karen_CS_to_Thai.pdf")
    # 1stKaren to 2ndGenHmong
    multiplot.relative.CS(mains, dms, map, ylab="Distance to 2nd-Generation Hmong", xlab=xlab, x.var=x.var, query_samples=karen_firstgen_cs, 
                        ref_samples=hmong_secondgen_cs, outputfn="Karen_CS_to_USbornHmong.pdf")
    # 1stHmong to HmongThai
    multiplot.relative.CS(mains, dms, map, ylab="Distance to Hmong Thai", xlab=xlab, x.var=x.var, query_samples=hmong_firstgen_cs, 
                        ref_samples=hmongthai, outputfn="Hmong_CS_to_Thai.pdf")
    # 1stHmong to 2ndHmong
    multiplot.relative.CS(mains, dms, map, ylab="Distance to US-born Hmong", xlab=xlab, x.var=x.var, query_samples=hmong_firstgen_cs, 
                        ref_samples=hmong_secondgen_cs, outputfn="Hmong_CS_to_USbornHmong.pdf")

### MB americanization (relative distance to reference group) - Lean vs. Obese only
    # 1. Hmong1st and Hmong2nd to HmongThai    
    # 2. all US to all Thai
    # 3. Karen1st to KarenThai
    query_samples_list <- list(c(hmong_secondgen_cs, hmong_firstgen_cs), c(hmong_secondgen_cs,hmong_firstgen_cs, karen_firstgen_cs), karen_firstgen_cs)
    ref_samples_list <- list(hmongthai, c(hmongthai, karenthai), karenthai)
    fn_ext <- c("Hmong", "All", "Karen")
    for(i in 1:length(query_samples_list))
    {
        map_temp <- map[query_samples_list[[i]],]
        map_temp <- map_temp[map_temp$BMI.Class %in% c("Lean","Obese"),]
        new_query_samples <- rownames(map_temp)
        # calculate the relative distances to reference groups as response variable
        rel.dists <- lapply(dms, function(xx) get.relative.distance(new_query_samples, ref_samples_list[[i]], xx))    
        multiplot.boxplot.by.group.x.bmi(map_temp, rel.dists, rep("Distance to Thai Reference",3), mains, outputfn = paste0("boxplot-BMI-x-RelativeDistance-",fn_ext[i],".pdf"))
    }

### MB americanization (BP Ratio) - Lean vs. Obese only

    bp <- log10(get.taxa.ratio(taxa, bug1=bacteroides, bug2=prevotella)) # calculate bp for everyone
    p <- plot.boxplot.by.group.x.bmi(map[cs,], bp[cs], "B-P Ratio", main="") # plot cross sectional only
    save_plot("boxplot-BMI-x-BPRatio.pdf", p, base_aspect_ratio = 1.3)

    p <- plot.boxplot.by.group.x.bmi(map[c(karenthai,karen_firstgen_cs),], bp[c(karenthai,karen_firstgen_cs)], "B-P Ratio", main="") # plot cross sectional only
    save_plot("boxplot-BMI-x-BPRatio-Karen.pdf", p, base_aspect_ratio = 1.3)

    p <- plot.boxplot.by.group.x.bmi(map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], bp[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs)], "B-P Ratio", main="") # plot cross sectional only
    save_plot("boxplot-BMI-x-BPRatio-Hmong.pdf", p, base_aspect_ratio = 1.3)


### plot differential taxa by BMI classes

    #iterate through all combos of looking at diff taxa in lean vs. obese by subgroup
    groups <- c("HmongThai","Hmong1st","Hmong2nd","KarenThai","Karen1st")
    taxatables <- list(taxa, taxa_L2, taxa_L6_rare)
    names(taxatables) <- c("L6", "L2", "L6_rare")
    combos <- expand.grid(names(taxatables),groups, stringsAsFactors=F)
    colnames(combos) <- c("taxa","group")
    combos[combos$group %in% c("KarenThai","HmongThai", "Hmong2nd"), "controls"] <- "Age"
    combos[combos$group %in% c("Karen1st","Hmong1st"), "controls"] <- c("Age,Years.in.US")

    sig.level=.25    
    lean.obese.cs.map <- map[map$BMI.Class %in% c("Lean", "Obese") & (is.na(map$Sample.Order) | map$Sample.Order==1),]
    
    for(i in 1:nrow(combos))
    {
        combo.name <- paste("BMI", paste(unlist(combos[i,c("group","taxa")]), collapse="."), sep=".")
        print(combo.name) # so we know which iteration we are in
        plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Sample.Group == combos$group[i],], taxatables[[combos$taxa[i]]], x.var="BMI.Class", 
            control.vars=unlist(strsplit(combos$controls[i],",")), outputfn.prepend=combo.name, sig.level=sig.level)
    }

    # create factor version of years in us so that data isn't omitted for NAs or 0s - to be used only when dealing with the entire dataset
    Years.in.US.Factor <- lean.obese.cs.map$Years.in.US
    Years.in.US.Factor[Years.in.US.Factor == 0] <- 50 # 2nd gen are 0, let's set them to something higher like 50??
    Years.in.US.Factor[is.na(Years.in.US.Factor)] <- 0 # Thai are NA, let's set them to 0
    
    Years.in.US.Factor <- factor(cut(Years.in.US.Factor,c(-1, .0001, seq(5,35,5), 41, 50)), ordered=T)
    
    lean.obese.cs.map[,"Years.in.US.Factor"] <- Years.in.US.Factor
    
    # try all the data, but add ethnicity as a control
    plot.diff.taxa(lean.obese.cs.map, taxa, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="BMI.all", sig.level=.10)

    plot.diff.taxa(lean.obese.cs.map, taxa_L6_rare, x.var="BMI.Class", 
        control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L6_rare", sig.level=.10)

    plot.diff.taxa(lean.obese.cs.map, taxa_L2, x.var="BMI.Class", 
        control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L2", sig.level=1)

### plot differential taxa by WHR

    plot.diff.taxa(lean.obese.cs.map, taxa, x.var="Waist.Height.Ratio", 
            control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="WHR.all", sig.level=.10)

### plot differential OFUs by BMI class
    output.table <- NULL
    for(i in 1:length(oFu_list))
    {
        ofu <- oFu_list[[i]]
        ret <- plot.diff.taxa(lean.obese.cs.map, ofu, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend=paste0("BMI.all.", names(oFu_list)[i]), sig.level=.25)
    
        if(nrow(ret) > 0)    
            output.table <- rbind(output.table, cbind(ret, similarity=names(oFu_list)[i], group="All Groups"))
    }
    # let's look at 1st gen immigrants only
    submap <- lean.obese.cs.map[lean.obese.cs.map$Sample.Group %in% c("Karen1st", "Hmong1st"),]
    for(i in 1:length(oFu_list))
    {
        ofu <- oFu_list[[i]]
        ret <- plot.diff.taxa(submap, ofu, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend=paste0("BMI.1stgen.", names(oFu_list)[i]), sig.level=.25)

        if(nrow(ret) > 0)
            output.table <- rbind(output.table, cbind(ret, similarity=names(oFu_list)[i], group="1stGen"))
    }
    write.table(output.table, file="IMP-OFU-lean-v-obese.txt", quote=F, row.names=F, sep="\t")


### plot differential taxa M1 vs M6 in longitudinal subjects
    ### nothing significant
    s1 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 1 & map$Subject.ID != "IMP.000",]))
    s2 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 6 & map$Subject.ID != "IMP.000",]))
    taxa0 <- taxa[c(s1,s2),]
    prevalences <- apply(taxa0, 2, function(bug.col) mean(bug.col > 0))
    taxa0 <- taxa0[, prevalences >= .10]
    ret <- collapse.by.correlation(taxa0, .95)
    taxa0 <- taxa0[, ret$reps]

    ret <- test.features.nonparametric.paired(taxa0, samples1=s1, samples2=s2, sig.level=.10)


### heatmap

    make.heatmap(otu, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.all.thaibaseline.pdf")
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.pdf")   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.karen.pdf") 

    make.heatmap(otu, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.all.thaibaseline.ordernumotus.pdf", order.by.num.otus=TRUE)
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.hmong.ordernumotus.pdf", order.by.num.otus=TRUE)   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=NULL, outputfn="heatmap.karen.ordernumotus.pdf", order.by.num.otus=TRUE) 


    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("HmongThai"), outputfn="heatmap.hmong.thaibaseline.pdf")   
    make.heatmap(otu, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("Hmong2nd"), outputfn="heatmap.hmong.2ndgenbaseline.pdf")   
    make.heatmap(otu, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=c("KarenThai"), outputfn="heatmap.karen.thaibaseline.pdf") 

    make.heatmap(taxa, map[cs,], .25, presence.absence=T, baseline.groups=c("HmongThai","KarenThai"), outputfn="heatmap.taxa.all.thaibaseline.pdf")
    make.heatmap(taxa, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("HmongThai"), outputfn="heatmap.taxa.hmong.thaibaseline.pdf")   
    make.heatmap(taxa, map[c(hmongthai,hmong_firstgen_cs, hmong_secondgen_cs),], .25, presence.absence=T, baseline.groups=c("Hmong2nd"), outputfn="heatmap.taxa.hmong.2ndgenbaseline.pdf")   
    make.heatmap(taxa, map[c(karenthai,karen_firstgen_cs),], .25, presence.absence=T, baseline.groups=c("KarenThai"), outputfn="heatmap.taxa.karen.thaibaseline.pdf") 
    
    