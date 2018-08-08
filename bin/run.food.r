# load food data
    food_wuf_dm <- read.table(paste(food_datadir, "wuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_uwuf_dm <- read.table(paste(food_datadir,"uwuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_bc_dm <- read.table(paste(food_datadir,"bc_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_otu_L3 <- read.table(paste(food_datadir,"food.otu_L3.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment="", skip=1)
    food_alpha_fn <- paste(food_datadir,"food.alpha.txt",sep="/")

    t_food_otu <- t(food_otu_L3)
    food_otu_L3 <- sweep(t_food_otu, 1, rowSums(t_food_otu), '/')
    food_alphadiv <- read.table(food_alpha_fn, sep="\t", quote="", row=1, head=T, comment.char="")
    food_alphadiv_L3 <- read.table(paste(food_datadir,"food.alpha_L3.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")

    food_otu_L5 <- t(read.table(paste(food_datadir,"food.otu_L5.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="", skip=1))


# nutrient file should already be mapped to sample ids
    nutrients <- read.table(nutrientsfn, sep="\t", header=T, check.names=F, as.is=T,row=1)
    nutrients[,"% of Calories from Total Sugars"] <- nutrients[,"Total Sugars in Grams"] * 4 / nutrients[,"Total Calories"] * 100
    nutrients[,"g Fiber per 1000 Calories"] <- (nutrients[,"Dietary Fiber in Grams",] / nutrients[,"Total Calories"]) * 1000


### ALPHA DIVERSITY
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(food_alphadiv[cs, alpha.metrics]), 
            ylabs=rep("",length(alpha.metrics)), mains=alpha.metrics, outputfn="food_boxplot.alphadiv.bmi.pdf", add.pvals=TRUE)

    p <- plot.scatter(x=food_alphadiv[cs,"shannon"], y=alphadiv[cs,"shannon"], xlab="Diet diversity", ylab="MB diversity", main="MB-Diet alpha diversity")
    save_plot("food-mb-alphadiv_shannon.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1 )

    p <- plot.scatter(x=food_alphadiv_L3[cs,"shannon"], y=alphadiv[cs,"shannon"], xlab="Diet diversity", ylab="MB diversity", main="MB-Diet alpha diversity")
    save_plot("food-mb-alphadiv_shannon_L3.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1 )


    map_L <- map[map$Sub.Study == "L" & map$Sample.Order %in% 0:6,] # there's one person with a 7th month sample, but we'll drop it
    
    ### food alpha diversity over time
    multiplot.alphadiv.L(map_L, food_alphadiv, alpha.metrics, outputfn="alphadiv.L.Month.pdf") # defaults to Sample Order
    multiplot.alphadiv.L(map_L, food_alphadiv, alpha.metrics, outputfn="alphadiv.L.first2.last2.pdf", num.clip.months=2)
    multiplot.alphadiv.L(map_L, food_alphadiv, alpha.metrics, outputfn="alphadiv.L.first.last.pdf", num.clip.months=1)

### Taxa summaries
    plot.taxa.summary(map0=map[cs,], otu=food_otu_L3, fn="food.summary.pdf")
    
    # plot ALL longitudinal food summaries
    plot.taxa.summary.L(taxa0=food_otu_L3, map0=map[map$Sub.Study=="L", ], outputfn="food.summary.L.pdf", grid.ncol=4)

    # look at individual summaries under a specific food category only
    food_otu <- read.table(paste(datadir,"../food.otu.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment="", strip.white=T)
    foods <- as.character(food_otu[,ncol(food_otu)])
    food_otu <- food_otu[,-ncol(food_otu)]
    t_food_otu <- t(food_otu)
    food_otu <- sweep(t_food_otu, 1, rowSums(t_food_otu), '/')
    foods_dm <- strsplit(foods, split=";")
    foods_dm <- data.frame(matrix(unlist(foods_dm), byrow=T, ncol=6),stringsAsFactors=FALSE, row.names = colnames(food_otu))
    foods_L1 <- c("L1_Meat_Poultry_Fish_and_Mixtures", "L1_Dry_Beans_Peas_Other_Legumes_Nuts_and_Seeds","L1_Grain_Product","L1_Fruits","L1_Vegetables","L1_Sugars_Sweets_and_Beverages")# doesn't make sense to look at all unique(foods_dm[,1])
    for(i in 1:length(foods_L1))
    {
        plot.taxa.summary.L(taxa0=food_otu[,rownames(foods_dm[foods_dm[,1] == foods_L1[i],])], map0=map[map$Sub.Study=="L" & map$Subject.ID != "IMP.000", ], 
                            x.var="Diet.Day.Since.Arrival", grid.ncol=4, outputfn=paste0("food.",foods_L1[i],".L.pdf"))

        plot.taxa.summary.L(taxa0=food_otu[,rownames(foods_dm[foods_dm[,1] == foods_L1[i],])], map0=map[map$Subject.ID %in% sequenced.subjects, ], 
                            x.var="Diet.Day.Since.Arrival", grid.ncol=2, outputfn=paste0("food.",foods_L1[i],".sequenced.L.pdf"))

    }

### PCOA
    #plot.pcoa(map[cs,], dm=food_euc_dm, plot.title="Food Euclidean")
    #plot.pcoa(map[cs,], dm=food_bc_dm, plot.title="Food Bray-Curtis")
    
    p <- plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", show.stats=F, flip.axis=1)
    save_plot("pcoa - Food Weighted Unifrac.pdf", p, base_aspect_ratio = 1)        

    p <- plot.pcoa(map[cs,], dm=food_uwuf_dm, plot.title="Food Unweighted Unifrac", show.stats=F)
    save_plot("pcoa - Food Unweighted Unifrac.pdf", p, base_aspect_ratio = 1)        

    plot.pcoa.long(map, dm=food_uwuf_dm, plot.title="Food UnWeighted Unifrac - L")
    plot.pcoa.long(map, dm=food_wuf_dm, plot.title="Food Weighted Unifrac - L")    
    
    
### Intra-inter group variabilities
    f_wuf_ret <- plot.within.group.distances(map0=map[cs,], food_wuf_dm, fn="within.group.food.wuf.pdf", ylab="Weighted Unifrac distance")
    f_uwuf_ret <- plot.within.group.distances(map0=map[cs,], food_uwuf_dm, fn="within.group.food.uwuf.pdf", ylab="Unweighted Unifrac distance")
    f_bc_ret <- plot.within.group.distances(map0=map[cs,], food_bc_dm, fn="within.group.food.bc.pdf", ylab="Bray-Curtis distance")

    plots <- list(f_wuf_ret$plot, wuf_ret$plot, f_uwuf_ret$plot, uwuf_ret$plot, f_bc_ret$plot, bc_ret$plot)
    # plot side by side MB and Food group distances
    multiplot <- plot_grid(plotlist=plots, ncol=2, nrow=3)
    
    # look at correlations of the variabilities 
    datalist <- list(f_wuf_ret$data, wuf_ret$data, f_uwuf_ret$data, uwuf_ret$data, f_bc_ret$data, bc_ret$data)
    cor.label <- NULL
    for(i in c(1,3,5))
    {
        valid <- intersect(rownames(datalist[[i]]), rownames(datalist[[i+1]]))
        cor.res <- cor.test(datalist[[i]][valid,"y"], datalist[[i+1]][valid,"y"])
        cor.label <- c(cor.label, substitute(paste(rho, " = ", estimate, ", P = ", pvalue),
                    list(estimate = signif(cor.res$estimate, 2), pvalue = signif(cor.res$p.value, 2))))
    }

    multiplot <- ggdraw(multiplot) + draw_label(cor.label[[1]], y=.98, hjust=.3) + draw_label(cor.label[[2]], y=.65, hjust=.3) + draw_label(cor.label[[3]], y=.32, hjust=.3) +
                draw_label(x=.05, y=.98, label="Food",fontface='bold',size=20,colour="blue3") + draw_label(x=.96, y=.98, label="MB",fontface='bold',size=20,colour="blue3") 
                    
    save_plot("food-x-mb.within.group.distances.pdf", multiplot,
          ncol = 2, nrow = 3, base_aspect_ratio = 1.3)
    
### NUTRIENTS
    # plot nutrients for ALL dietary data collected, by Sample Group only (NOT over Yrs in US)
        # find all of the outliers!
        #         outliers <- NULL
        #         for(nut in nutrient.vars)
        #         {
        #             outlier <- boxplot(nutrients[,nut])$out
        #             outliers <- rbind(outliers, 
        #                         data.frame(SuperTracker.DATE=nutrients[which(nutrients[,nut] %in% outlier), "SuperTracker.DATE"], 
        #                                 Sample.ID=rownames(nutrients)[which(nutrients[,nut] %in% outlier)], 
        #                                 Nutrient=nut, stringsAsFactors=F))
        # 
        #         }
        #         write.table(outliers[order(outliers$SuperTracker.DATE),], file="outliers.txt", sep="\t", quote=F, row.names=F)

    nutrient.vars <- c("Total Calories", "% of Calories from Total Sugars", "% of Calories from Added Sugars","% of Calories from Carbohydrate","% of Calories from Protein",
    "% of Calories from Saturated Fat", "% of Calories from Total Fat","g Fiber per 1000 Calories")
    
    p <- lapply(nutrient.vars, function(nut)  map.boxplot(y=nutrients[cs,nut], Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=F, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=10, x.size=9))
    save_plot("macronutrients.pdf", plot_grid(plotlist=p, ncol=4, nrow=2), base_aspect_ratio=3)
                                                    
    p <- lapply(nutrient.vars, function(nut)  map.boxplot(y=nutrients[cs,nut], Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=T, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=10, x.size=9))
    save_plot("macronutrients-pval.pdf", plot_grid(plotlist=p, ncol=4, nrow=2), base_aspect_ratio=3)

    p <- lapply(nutrient.vars, function(nut)  map.boxplot(y=log10(nutrients[cs,nut]+1), Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=F, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=10, x.size=9))
    save_plot("macronutrients-log.pdf", plot_grid(plotlist=p, ncol=4, nrow=2), base_aspect_ratio=3)

    p <- lapply(nutrient.vars, function(nut)  map.boxplot(y=log10(nutrients[cs,nut]+1), Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=T, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=10, x.size=9))
    save_plot("macronutrients-log-pval.pdf", plot_grid(plotlist=p, ncol=4, nrow=2), base_aspect_ratio=3)


### Relative Food Distances plots

    ### Longitudinal participants
    mains <- c("Unweighted Unifrac","Weighted Unifrac", "Bray-Curtis")
    food_dms <- list(food_uwuf_dm, food_wuf_dm, food_bc_dm)
    xlab<-"Days Since US Arrival"
    lsamples <- rownames(map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",])
    lsamples <- intersect(lsamples, rownames(diet_map)) # these are L samples with diet entries

    # FOOD relative distances of longitudinal - by days to first diet record
    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Dissimilarity to First Sample", xlab=xlab, x.var="Diet.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Diet - Karen_L_to_Day0.pdf")

    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Dissimilarity to First Diet", xlab="Month", x.var="Sample.Order",
                        ref.sample.order=1, outputfn="Diet - Karen_L_to_Month0.pdf")


    # FOOD relative distances of longitudinal - by days to last diet record
    # by days to last sample    
    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Dissimilarity to Last Sample", xlab=xlab, x.var="Diet.Day.Since.Arrival",
                        ref.sample.order=6, outputfn="Diet - Karen_L_to_Day_M6.pdf")

    # FOOD relative distances of longitudinal - by days to first diet record
    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Dissimilarity to First Sample", xlab=xlab, x.var="Diet.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Diet - Karen_L_to_Previous.pdf", to.previous=T)

    # L to KarenThai
    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Distance to Karen in Thailand", 
                        xlab=xlab, x.var="Diet.Day.Since.Arrival", ref_samples=karenthai, outputfn="Diet - Karen_L_to_KarenThai.pdf")
    # L to 2ndGenHmong
    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Distance to 2nd-Generation", 
                        xlab=xlab, x.var="Diet.Day.Since.Arrival", ref_samples=hmong_secondgen_cs, outputfn="Diet - Karen_L_to_USborn.pdf")

    ### CS participants
    xlab="Years in US"
    x.var="Years.in.US"
    
    # all 1st gen to all thai
    multiplot.relative.CS(mains, food_dms, map, ylab="Distance to all Thai", xlab=xlab, x.var=x.var, query_samples=c(hmong_firstgen_cs, karen_firstgen_cs), 
                        ref_samples=c(karenthai,hmongthai), outputfn="Diet - allfirstgen_to_Thai.pdf")
    # 1stKaren KarenThai
    multiplot.relative.CS(mains, food_dms, map, ylab="Distance to Karen Thai", xlab=xlab, x.var=x.var, query_samples=karen_firstgen_cs, 
                        ref_samples=karenthai, outputfn="Diet - Karen_CS_to_Thai.pdf")
    # 1stKaren to 2ndGenHmong
    multiplot.relative.CS(mains, food_dms, map, ylab="Distance to 2nd-Generation Hmong", xlab=xlab, x.var=x.var, query_samples=karen_firstgen_cs, 
                        ref_samples=hmong_secondgen_cs, outputfn="Diet - Karen_CS_to_USbornHmong.pdf")
    # 1stHmong to HmongThai
    multiplot.relative.CS(mains, food_dms, map, ylab="Distance to Hmong Thai", xlab=xlab, x.var=x.var, query_samples=hmong_firstgen_cs, 
                        ref_samples=hmongthai, outputfn="Diet - Hmong_CS_to_Thai.pdf")
    # 1stHmong to 2ndHmong
    multiplot.relative.CS(mains, food_dms, map, ylab="Distance to US-born Hmong", xlab=xlab, x.var=x.var, query_samples=hmong_firstgen_cs, 
                        ref_samples=hmong_secondgen_cs, outputfn="Diet - Hmong_CS_to_USbornHmong.pdf")


### FOOD AND MB correlations

    ### Procrustes: Food and Microbiome 
        draw.procrustes(map0=map[cs,], plot.title="procrustes-mb.uwuf-food.uwuf", max.axis=300, dm1=food_uwuf_dm, dm2=uwuf_dm, dm1.name="Food", dm2.name="MB")
        draw.procrustes(map0=map[cs,], plot.title="procrustes-mb.wuf-food.uwuf", max.axis=300, dm1=food_uwuf_dm, dm2=wuf_dm, dm1.name="Food", dm2.name="MB")
        draw.procrustes(map0=map[cs,], plot.title="procrustes-mb.wuf-food.wuf", max.axis=300, dm1=food_wuf_dm, dm2=wuf_dm, dm1.name="Food", dm2.name="MB")

    ### Scatterplot of Food and MB distances - this is ugly and not very useful - even after permuting it still looks about the same
        plot.two.dms(map[cs,], plot.title="scatter-uwuf.mb-vs-uwuf.food", dm1=food_uwuf_dm, dm2=uwuf_dm, xlab="Food", ylab="Microbiome")
        plot.two.dms(map[cs,], plot.title="scatter-wuf.mb-vs-wuf.food", dm1=food_wuf_dm, dm2=wuf_dm, xlab="Food", ylab="Microbiome")
        
    ### Relative Food AND MB Distances - this doesn't do as well
        # this gets the relative MB AND Food distance of each sample in a query group to to a reference group, and tests if they're correlated
        rel.distance <- get.relative.distance(query_samples=cs[!(cs %in% controls)], ref_samples=controls, dm=wuf_dm)
        food.rel.distance <- get.relative.distance(query_samples=cs[!(cs %in% controls)], ref_samples=controls, dm=food_wuf_dm)
        p <- plot.scatter(x=food.rel.distance, y=rel.distance, xlab="Diet distance", ylab="MB distance", main="Distance to Controls")
        save_plot("food-v-mb-distance-to-controls_wuf.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1 )

        rel.distance <- get.relative.distance(query_samples=cs[!(cs %in% controls)], ref_samples=controls, dm=uwuf_dm)
        food.rel.distance <- get.relative.distance(query_samples=cs[!(cs %in% controls)], ref_samples=controls, dm=food_uwuf_dm)
        p <- plot.scatter(x=food.rel.distance, y=rel.distance, xlab="Diet distance", ylab="MB distance", main="Distance to Controls")
        save_plot("food-v-mb-distance-to-controls_uwuf.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1 )
    
    ### Food PC1 AND Relative MB Distances - this does better
        plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", save.pc=TRUE, show.stats=FALSE)
        foodpc <- read.table("Food Weighted Unifrac-PC.txt", sep="\t", row=1)
        rel.distance <- get.relative.distance(query_samples=cs[!(cs %in% controls)], ref_samples=controls, dm=wuf_dm)    
        p <- plot.scatter(x=foodpc[cs[!(cs %in% controls)], 1], y=rel.distance, xlab="Diet PC1", ylab="MB distance", main="Distance to Controls")
        save_plot("foodPC1-v-mb-distance-to-controls_wuf.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1 )

    ### MB and Food stability - not super interesting
        subs <- unique(map_L$Subject.ID)
        food.dists <- list()
        mb.dists <- list()
        for(subject in subs)
        {
            food.dists[[subject]] <- mean(get.within.dist(dm=food_uwuf_dm, samples=rownames(map_L[map_L$Subject.ID == subject,])))
            mb.dists[[subject]] <- mean(get.within.dist(dm=uwuf_dm, samples=rownames(map_L[map_L$Subject.ID == subject,])))
        }

        ggdata <- data.frame(food=unlist(food.dists), mb=unlist(mb.dists))
        ggdata$food.stab <- 1-ggdata$food
        ggdata$mb.stab <- 1-ggdata$mb
        # 1/(1+dist) for WUF?
    
        p <- plot.scatter(x=ggdata$food.stab, y=ggdata$mb.stab, xlab="Diet stability", ylab="MB stability", main="UWUF")
        save_plot("food-mb-stability_L.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1 )


    ### see what L3 foods are correlated with first 5 food PCs
        plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", save.pc=TRUE, show.stats=FALSE, axis2=5)
        foodpc <- read.table("Food Weighted Unifrac-PC.txt", sep="\t", row=1)
        
        foodotu <- food_otu_L3[rownames(foodpc),]
        prevalences <- apply(foodotu, 2, function(bug.col) mean(bug.col > 0))
 		foodotu <- foodotu[, prevalences >= .05]

        food.corr <- list()
        for(i in 1:ncol(foodpc))
        {
            pval <- NULL
            rho <- NULL
            for(j in 1:ncol(foodotu))
            {
                ret <- cor.test(foodpc[,i], foodotu[,j], method="spear", exact=F)
                pval[j] <- ret$p.value
                rho[j] <- ret$estimate
            }
            adj.pval <- p.adjust(pval, method="fdr")
            food.corr[[i]] <- data.frame(pval, adj.pval, rho, foodname=colnames(foodotu))
        }
        pc1 <- food.corr[[1]][food.corr[[1]]$adj.pval < .05,]
        pc1 <- pc1[order(pc1$adj.pval),]
        write.table(pc1,file="pc1-foods.sig.txt", sep="\t", quote=F, row.names=F)
        pc2 <- food.corr[[2]][food.corr[[2]]$adj.pval < .05,]
        pc2 <- pc2[order(pc2$adj.pval),]
        write.table(pc2,file="pc2-foods.sig.txt", sep="\t", quote=F, row.names=F)        

############ ARCHIVED EXPLORATORY STUFF ###################

# plot nutrients over years in the US for sequenced samples only
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen"), nutrients0=nutrients)
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen"), nutrients0=nutrients, independent.var="Fraction.Life.in.US")

# plot food groups over years in the US
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen"), nutrients0=foodgroups, nutrient.vars=colnames(foodgroups))
#     plot.nutrients(map0=map[cs,], ethnicity=c("Karen"), nutrients0=foodgroups, bins=0:10,nutrient.vars=colnames(foodgroups))
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong"), nutrients0=foodgroups, bins=seq(0,45,5),nutrient.vars=colnames(foodgroups))

#     ### plot PCOA and color by nutrients
#     plot.pcoa.nutrient(map0=map[cs,], pc0=pc, nutrients0=nutrients, ethnicity=c("Hmong","Karen"), fn="pcoa_nutrients.pdf")
#     plot.pcoa.nutrient(map0=map[cs,], pc0=pc, nutrients0=nutrients, ethnicity="Hmong", fn="pcoa_nutrients_hmong.pdf")
#     plot.pcoa.nutrient(map0=map[cs,], pc0=pc, nutrients0=nutrients, ethnicity="Karen", fn="pcoa_nutrients_karen.pdf")
# 
#     ### plot PCOA and subset groups of interest before coloring nutrients
#     plot.pcoa.nutrient(map0=map[c(hmong_firstgen_cs, hmong_secondgen_cs),], pc0=pc, nutrients0=nutrients, ethnicity="Hmong", fn="pcoa_nutrients_hmong_US.pdf")
#     plot.pcoa.nutrient(map0=map[hmong_firstgen_cs,], pc0=pc, nutrients0=nutrients, ethnicity="Hmong", fn="pcoa_nutrients_hmong_1st.pdf")
#     plot.pcoa.nutrient(map0=map[hmongthai,], pc0=pc, nutrients0=nutrients, ethnicity="Hmong", fn="pcoa_nutrients_hmong_thai.pdf")
#     plot.pcoa.nutrient(map0=map[karen_firstgen_cs,], pc0=pc, nutrients0=nutrients, ethnicity="Karen", fn="pcoa_nutrients_karen_US.pdf")
#     plot.pcoa.nutrient(map0=map[karenthai,], pc0=pc, nutrients0=nutrients, ethnicity="Karen", fn="pcoa_nutrients_karen_thai.pdf")

#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen","Caucasian"), nutrients0=nutrients, nutrient.vars=nutrient.vars, independent.var="Sample.Group")
#     plot.nutrients(map0=map[cs,], ethnicity=c("Karen"), nutrients0=nutrients, bins=0:10,nutrient.vars=nutrient.vars, independent.var="Sample.Group")
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong"), nutrients0=nutrients, bins=seq(0,45,5),nutrient.vars=nutrient.vars, independent.var="Sample.Group")
# 
#     ### plot nutrients for all data collected -- over Decades in US subgrouped by BMI.class
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen"), nutrients0=nutrients, independent.var="BMI.Class", nutrient.vars=nutrient.vars)
#     plot.nutrients(map0=map[cs,], ethnicity=c("Karen"), nutrients0=nutrients, independent.var="BMI.Class", bins=0:10, nutrient.vars=nutrient.vars)
#     plot.nutrients(map0=map[cs,], ethnicity=c("Hmong"), nutrients0=nutrients, independent.var="BMI.Class", bins=seq(0,45,5), nutrient.vars=nutrient.vars)
