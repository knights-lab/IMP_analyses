# load food data
    food_wuf_dm <- read.table(paste(food_datadir, "wuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_uwuf_dm <- read.table(paste(food_datadir,"uwuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_euc_dm <- read.table(paste(food_datadir,"euc_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_bc_dm <- read.table(paste(food_datadir,"bc_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_otu_L3 <- read.table(paste(food_datadir,"food.otu_L3.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment="", skip=1)
    t_food_otu <- t(food_otu_L3)
    food_otu_L3 <- sweep(t_food_otu, 1, rowSums(t_food_otu), '/')


# nutrient file should already be mapped to sample ids
    nutrients <- read.table(nutrientsfn, sep="\t", header=T, check.names=F, as.is=T,row=1)
    nutrients[,"% of Calories from Total Sugars"] <- nutrients[,"Total Sugars in Grams"] * 4 / nutrients[,"Total Calories"]
    nutrients[,"g Fiber per 1000 Calories"] <- (nutrients[,"Dietary Fiber in Grams",] / nutrients[,"Total Calories"]) * 1000


### ALPHA DIVERSITY
    # where is the food alpha file?!
    # plot.alphadiv(map[cs,], food_alpha[cs,], metric = "PD_whole_tree")

#     hmong_and_controls <- c(controls,hmong_secondgen_cs, hmong_firstgen_cs, hmongthai)
#     karen_and_controls <- c(controls, karen_firstgen_cs, karenthai)
#     plot.alphadiv(map[hmong_and_controls,], food_alpha[hmong_and_controls,], metric = "PD_whole_tree")
#     plot.alphadiv(map[karen_and_controls,], food_alpha[karen_and_controls,], metric = "PD_whole_tree")
# 
#     plot.alphadiv(map[cs,], food_alpha[cs,], metric = "PD_whole_tree")    
#     plot.alphadiv(map[cs,], food_alpha[cs,], metric = "shannon")
#     plot.alphadiv(map[cs,], food_alpha[cs,], metric = "simpson")

#look at consumption of rice!
#L1_Grain_Product;L2_Pastas_cooked_cereals_rice;L3_Cooked_cereals_rice

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
    
    # let's only look at meat
#     food_otu <- read.table(paste(datadir,"../food.otu.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment="", strip.white=T)
#     foods <- as.character(food_otu[,ncol(food_otu)])
#     food_otu <- food_otu[,-ncol(food_otu)]
#     food_otu <- t(food_otu)
#     foods_dm <- strsplit(foods, split=";")
#     foods_dm <- data.frame(matrix(unlist(foods_dm), byrow=T, ncol=6),stringsAsFactors=FALSE, row.names = colnames(food_otu))
# 
#     meatnames <- rownames(foods_dm[foods_dm$X1 == "L1_Meat_Poultry_Fish_and_Mixtures" & foods_dm$X2 != "L2_Meatpoultry_fish_with_nonmeat",])
#     meat_otu <- food_otu[,meatnames]
#     valid_samples <- intersect(cs, rownames(meat_otu))
#     temp_map <- map[valid_samples,]
#     plot(map[valid_samples, "Years.in.US"], rowSums(meat_otu[valid_samples,]))
# 
#     plot(temp_map[temp_map$Years.in.US != 0, "Years.in.US"], rowSums(meat_otu[rownames(temp_map[temp_map$Years.in.US != 0,]),]))


### PCOA
    plot.pcoa(map[cs,], dm=food_euc_dm, plot.title="Food Euclidean")
    plot.pcoa(map[cs,], dm=food_bc_dm, plot.title="Food Bray-Curtis")
    plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac")
    plot.pcoa(map[cs,], dm=food_uwuf_dm, plot.title="Food Unweighted Unifrac")
    
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
    nutrient.vars <- c("Total Calories", "% of Calories from Total Sugars", "% of Calories from Added Sugars","% of Calories from Carbohydrate","% of Calories from Protein",
    "% of Calories from Saturated Fat", "% of Calories from Total Fat","g Fiber per 1000 Calories")
    plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen","Caucasian"), nutrients0=nutrients, nutrient.vars=nutrient.vars, independent.var="Sample.Group")
    plot.nutrients(map0=map[cs,], ethnicity=c("Karen"), nutrients0=nutrients, bins=0:10,nutrient.vars=nutrient.vars, independent.var="Sample.Group")
    plot.nutrients(map0=map[cs,], ethnicity=c("Hmong"), nutrients0=nutrients, bins=seq(0,45,5),nutrient.vars=nutrient.vars, independent.var="Sample.Group")

    ### plot nutrients for all data collected -- over Decades in US subgrouped by BMI.class
    plot.nutrients(map0=map[cs,], ethnicity=c("Hmong","Karen"), nutrients0=nutrients, independent.var="BMI.Class", nutrient.vars=nutrient.vars)
    plot.nutrients(map0=map[cs,], ethnicity=c("Karen"), nutrients0=nutrients, independent.var="BMI.Class", bins=0:10, nutrient.vars=nutrient.vars)
    plot.nutrients(map0=map[cs,], ethnicity=c("Hmong"), nutrients0=nutrients, independent.var="BMI.Class", bins=seq(0,45,5), nutrient.vars=nutrient.vars)

### Relative Food Distances plots

    ### Longitudinal participants
    mains <- c("Unweighted Unifrac","Weighted Unifrac", "Bray-Curtis")
    food_dms <- list(food_uwuf_dm, food_wuf_dm, food_bc_dm)
    xlab<-"Days Since US Arrival"
    diet_map <- nutrients[,1:7]
    lsamples <- rownames(map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",])
    lsamples <- intersect(lsamples, rownames(diet_map)) # these are L samples with diet entries

    # FOOD relative distances of longitudinal - by days to first diet record
    multiplot.relative.L(mains=mains, dms=food_dms, map = map[lsamples,], ylab="Dissimilarity to First Sample", xlab=xlab, x.var="Diet.Day.Since.Arrival",
                        ref.sample.order=1, outputfn="Diet - Karen_L_to_Day0.pdf")
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



### Relative Food AND MB Distances -- OLD (don't do these anymore - won't be as good)
    # this gets the relative MB AND Food distance of each sample in a query group to to a reference group, and tests if they're correlated
#     plot.relative.distance.with.food(dm=wuf_dm, food_dm=food_wuf_dm, query_samples=c(hmong_firstgen_cs,hmong_secondgen_cs), 
#     food_ref_samples=controls, mb_ref_samples=c(hmongthai))
# 
#     plot.relative.distance.with.food(dm=uwuf_dm, food_dm=food_uwuf_dm, query_samples=c(hmong_firstgen_cs,hmong_secondgen_cs, karen_firstgen_cs), 
#     food_ref_samples=c(hmongthai, karenthai), mb_ref_samples=c(hmongthai, karenthai))
# 
#     plot.relative.distance.with.food(dm=uwuf_dm, food_dm=food_uwuf_dm, query_samples=c(karen_firstgen_cs), 
#     food_ref_samples=karenthai, mb_ref_samples=karenthai)
# 
#     plot.relative.distance.with.food(dm=uwuf_dm, food_dm=food_uwuf_dm, query_samples=c(hmongthai, hmong_firstgen_cs), 
#     food_ref_samples=hmong_secondgen_cs, mb_ref_samples=hmong_secondgen_cs)
# 
#     plot.relative.distance.with.food(dm=wuf_dm, food_dm=food_wuf_dm, query_samples=c(hmong_firstgen_cs,hmong_secondgen_cs), 
#     food_ref_samples=hmongthai, mb_ref_samples=hmongthai)
# 
#     # significant
#     plot.relative.distance.with.food(dm=uwuf_dm, food_dm=food_uwuf_dm, query_samples=c(hmong_firstgen_cs,hmong_secondgen_cs), 
#     food_ref_samples=hmongthai, mb_ref_samples=hmongthai)

### Food PC1 AND Relative MB Distances
    plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", save.pc=TRUE)
    foodpc <- read.table("Food Weighted Unifrac-PC.txt", sep="\t", row=1)
    plot.relative.distance.with.food.pc(dm=uwuf_dm, food_pc=foodpc, query_samples=c(hmong_firstgen_cs,hmong_secondgen_cs, karen_firstgen_cs), 
                                      mb_ref_samples=c(hmongthai, karenthai))
    plot.relative.distance.with.food.pc(dm=uwuf_dm, food_pc=foodpc, query_samples=c(hmong_firstgen_cs,hmong_secondgen_cs), 
                                        mb_ref_samples=c(hmongthai))
    plot.relative.distance.with.food.pc(dm=uwuf_dm, food_pc=foodpc, query_samples=c(karen_firstgen_cs), 
                                        mb_ref_samples=c(karenthai))
    

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

