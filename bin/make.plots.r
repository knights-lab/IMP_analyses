setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

# load all source files, data, and prep any variables for reuse
source("bin/load.r")

### alpha diversity
    plot.alphadiv(map[cs,], alphadiv0[cs,], metric = "PD_whole_tree")
    plot.alphadiv(map[cs,], alphadiv0[cs,], metric = "shannon")
    plot.alphadiv(map[cs,], alphadiv0[cs,], metric = "simpson")
    
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
    plot.taxa.summary(map0=map[cs,], otu=taxa, fn="taxa.summary.pdf")    

    # plot self samples
    ssamples <- map[map$Subject.ID == "IMP.000", "Sample.Order", drop=F]
    ssamples.ordered <- rownames(ssamples[ order(ssamples[,1]), , drop=F]) 
    day.labels <- as.character(sort(ssamples[,1])-3)
    day.labels[(length(day.labels)-3):length(day.labels)] <- paste0("P",1:4)
    plot.taxa.summary(map0=map[ssamples.ordered,], otu=taxa, fn="taxa.summary.IMP000.pdf", sample.order=ssamples.ordered,
                            x.labels = day.labels)
    
### PCOA
    plot.pcoa.long(map_all, uwuf_dm, "UnWeighted Unifrac - L")
    plot.pcoa.long(map_all, wuf_dm, "Weighted Unifrac - L")

    plot.pcoa(map[cs,], wuf_dm, "Weighted Unifrac")
    plot.pcoa(map[cs,], uwuf_dm, "Unweighted Unifrac")
    plot.pcoa(map[cs,], bc_dm, "Bray Curtis")

plot.pcoa.by(map0=map[cs,], wuf_dm, fn="KCK.fmt.pcoa.pdf", ethnicity="Hmong", color.by="Years.in.US", and.by="BMI.Class", label.samples=paste0("TFSCS0", 20:29))

### Bacteroides-Prevotella
    plot.b.p.barplot(map[c(karenthai,karen_firstgen_cs),], taxa, bins=0:10, fn="b.p.barplot.karen.pdf")
    plot.b.p.barplot(map[cs,], taxa, bins=seq(0,45,5), fn="b.p.barplot.pdf")

    # CS samples only
    plot.b.p.ratio(map[firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.pdf")
    plot.b.p.ratio(map[hmong_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.hmong.pdf")
    plot.b.p.ratio(map[karen_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.karen.pdf")

    # ALL samples 
    plot.b.p.ratio.all(map, taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.pdf", g1=firstgen_cs, g2=c(hmongthai,karenthai),g3=(hmong_secondgen_cs))

    # plot longitudinal samples only
    plot.b.p.ratio(map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.L.pdf", longitudinal=T)

### Prediction 
    prediction <- predict.years(map[map$Ethnicity=="Karen" & (is.na(map$Sample.Order) | map$Sample.Order==1),], taxa)

### Intra-inter group variabilities
    bc_ret <- plot.within.group.distances(map0=map[cs,], bc_dm, fn="within.group.bc.pdf", ylab="Bray-Curtis distance")
    wuf_ret <- plot.within.group.distances(map0=map[cs,], wuf_dm, fn="within.group.wuf.pdf", ylab="Weighted Unifrac distance")
    uwuf_ret <- plot.within.group.distances(map0=map[cs,], uwuf_dm, fn="within.group.uwuf.pdf", ylab="Unweighted Unifrac distance")

    plot.between.group.distances(map0=map[cs,], bc_dm, fn="between.group.bc.pdf", ylab="Bray-Curtis distance")
    plot.between.group.distances(map0=map[cs,], wuf_dm, fn="between.group.wuf.pdf", ylab="Weighted Unifrac distance")

### Relative Distance - Longitudinal Plots

    # by days since arrival
    mains <- c("Unweighted Unifrac","Weighted Unifrac", "Bray-Curtis")
    dms <- list(uwuf_dm,wuf_dm,bc_dm)
    xlab<-"Days Since US Arrival"
    # by days to first sample        
    multiplot.relative.L(mains=mains, dms=dms, map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Dissimilarity to First Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Karen_L_to_Day0.pdf")
    
    # by days to last sample    
    multiplot.relative.L(mains=mains, dms=dms, map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Dissimilarity to Last Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=6, outputfn="Karen_L_to_Day_M6.pdf")

    # days to first sample - day-to-day
    multiplot.relative.L(mains=mains, dms=dms, map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Dissimilarity to Previous Sample", xlab=xlab, x.var="Sample.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Karen_L_to_Previous.pdf", to.previous=T)
    
    # L to KarenThai
    multiplot.relative.L(mains=mains, dms=dms, map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Distance to Karen in Thailand", 
                        xlab=xlab, x.var="Sample.Day.Since.Arrival", ref_samples=karenthai, outputfn="Karen_L_to_KarenThai.pdf")
    # L to 2ndGenHmong
    multiplot.relative.L(mains=mains, dms=dms, map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], ylab="Distance to 2nd-Generation", 
                        xlab=xlab, x.var="Sample.Day.Since.Arrival", ref_samples=hmong_secondgen_cs, outputfn="Karen_L_to_USborn.pdf")

    # IMP000 to Day 0
    # let's make an extra var for IMP000
    map.000 <- map[map$Subject.ID == "IMP.000",]
    map.000$Sample.Day.Since.First.Sample <- as.numeric(as.Date(map.000$Sample.Date, format="%m/%d/%y") - as.Date("08/01/16", format="%m/%d/%y")) # hard code first sample date
    map.000[,"travel.phase"] <- "Traveling"
    map.000[map.000$Sample.Day.Since.First.Sample <= 2,"travel.phase"] <- "Pre"
    map.000[map.000$Sample.Day.Since.First.Sample >= 28,"travel.phase"] <- "Post"
    map.000$travel.phase <- factor(map.000$travel.phase, levels=c("Pre", "Traveling", "Post")) 
    
    # relative to first day
    multiplot.relative.L(mains=mains, dms=dms, map=map.000, ylab="Distance to Self at Day 0", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref.sample.order = 1, outputfn="IMP000_to_Day0.pdf", override.cols.by = "travel.phase")
    # overall variability -- relative to previous sample
    multiplot.relative.L(mains=mains, dms=dms, map=map.000, ylab="Distance to Previous Sample", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref.sample.order = 1, outputfn="IMP000_to_Previous.pdf", to.previous=T, override.cols.by="travel.phase")                
    # relative to Hmong Thai
    multiplot.relative.L(mains=mains, dms=dms, map=map.000, ylab="Distance to Hmong Thai", 
                    xlab="Day", x.var="Sample.Day.Since.First.Sample", ref_samples=hmongthai, outputfn="IMP000_to_HmongThai.pdf", override.cols.by="travel.phase")
    # relative to 2nd Gen Hmong
    multiplot.relative.L(mains=mains, dms=dms, map=map.000, ylab="Distance to 2nd Gen Hmong", 
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
        map_temp <- map_temp[map_temp$BMI.Class %in% c("Normal","Obese"),]
        new_query_samples <- rownames(map_temp)
        # calculate the relative distances to reference groups as response variable
        rel.dists <- lapply(dms, function(xx) get.relative.distance(new_query_samples, ref_samples_list[[i]], xx))    
        multiplot.boxplot.by.group.x.bmi(map_temp, rel.dists, rep("Distance to Thai Reference",3), mains, outputfn = paste0("boxplot-BMI-x-RelativeDistance-",fn_ext[i],".pdf"))
    }

### MB americanization (BP Ratio) - Lean vs. Obese only





















########## ARCHIVED EXPLORATORY STUFF ############
    
# this doesnt work well
# heatmap nutrients + microbes
#    plot.nutrient.mb.heatmap(nutrients=nutrients, taxa=taxa, map0=map[cs,], fn="nutrients_mb_heatmap.all2.pdf")

# why isnt this working!?
#plot.nutrient.mb.heatmap(nutrients=nutrients, taxa=taxa, map0=map[hmong_firstgen_cs,], fn="nutrients_mb_heatmap.hmongthai.pdf")

# color by BMI
#     plot.pcoa.BMI(map[c(hmong_firstgen_cs, hmong_secondgen_cs),], pc, "Hmong", "pcoa_BMI_Hmong_US.pdf")
#     plot.pcoa.BMI(map[c(hmong_firstgen_cs),], pc, "Hmong", "pcoa_BMI_Hmong_1st.pdf")
#     plot.pcoa.by(map[c(hmong_firstgen_cs),], pc, "Hmong", "pcoa_YearsInUS_Hmong_1st.pdf")
        
# plot food and taxa summaries, then order samples by each other -- # this isn't very useful
#     food.sample.ids <- intersect(rownames(food_otu), rownames(taxa[cs,]))
#     ordered.taxa <- plot.taxa.summary(map0=map[food.sample.ids,], otu=taxa, fn="taxa.summary2.pdf")    
#     ordered.food <- plot.taxa.summary(map0=map[food.sample.ids,], otu=food_otu, fn="food.summary2.pdf")
#     temp <- plot.taxa.summary(map0=map[ordered.food,], otu=taxa, fn="taxa.summary.by.food.pdf", sample.order=ordered.food)    
#     temp <- plot.taxa.summary(map0=map[ordered.taxa,], otu=food_otu, fn="food.summary.by.taxa.pdf", sample.order=ordered.taxa)
