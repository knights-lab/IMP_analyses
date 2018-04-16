datadirname <- "denovo.clr"
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("bin/load.r")

# ----- CROSS-SECTIONAL PLOTS ----- 

### PCOA

    # let's flip the axis so that Years in US goes left to right, BMI/Age goes bottom to up
    plot.pcoa(map[cs,], otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR", show.stats=F, flip.axis=1)    

    plot.pcoa.by(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7[cs,], fn="pcoa - CLR - AgeArrival.pdf", color.by="Age.at.Arrival")
    plot.pcoa.by(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7[cs,], fn="pcoa - CLR - Years.in.US.pdf", color.by="Years.in.US")


    # this isn't super useful but helps get the point across -- PCOA with Environment Vectors (simply maps the vectors onto an unconstrained ordination)
    plot.pcoa(map[cs,], otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - PCOA with Vectors", env.vars=c("BMI","Years.in.US","Age"), flip.axis=2)    

    # RDA: constrain ordination with selected env variables
    # here we want to look at the % relative variation explained by the constrained model VS the unconstrained model
    plot.constrained.ordination(map[cs,], otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Constrained", env.vars=c("Years.in.US","BMI","Age"))    
    plot.constrained.ordination(map[cs,], otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Unconstrained")    

    # RDA: with 1st gen only
    #plot.constrained.ordination(map[firstgen_cs,], otu0=taxa_L7[firstgen_cs,], method="euclidean", plot.title="Euclidean - CLR - 1stgen", env.vars=c("Years.in.US","BMI","Age"))    
    #plot.constrained.ordination(map[firstgen_cs,], otu0=taxa_L7[firstgen_cs,], method="euclidean", plot.title="Euclidean - CLR - 1stgen - Unconstrained")    

    # try plotting constrained ordination with 1st 5 diet PCs as included vars
    
    # generate food PCs here (weighted unifrac)
    plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", save.pc=TRUE, axis1 = 1, axis2 = 5) # put axis 5 so we get all 5 PCs back
    # load in food PCs
    food.pc <- read.table("Food Weighted Unifrac-PC.txt", sep="\t", skip=1, row=1)
    colnames(food.pc) <- paste0("Food.PC", 1:5)
    
    # find correlation between Food and MB distances
    valid <- intersect(rownames(food_wuf_dm),rownames(taxa_L7))
    food_wuf_dm <- food_wuf_dm[valid,valid]
    food_wuf_dist <- as.dist(food_wuf_dm)
    euc_dist <- as.matrix(vegdist(taxa_L7[valid,], method="euclidean"))
    mantel(euc_dist, food_wuf_dist)
    
    # redo constrained ordination with just metadata vars of interest
    plot.constrained.ordination(data.frame(map[cs,], food.pc[cs,]), otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Constrained Vars", env.vars=c("Years.in.US","BMI","Age", "Ethnicity"))    
    # redo constrained ordination with metadata vars AND food PCs
    plot.constrained.ordination(data.frame(map[cs,], food.pc[cs,]), otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Constrained Vars and Food PCs", env.vars=c("Years.in.US","BMI", "Ethnicity","Age",colnames(food.pc)))    
    # redo constrained ordination with food PCs only
    plot.constrained.ordination(data.frame(map[cs,], food.pc[cs,]), otu0=taxa_L7[cs,], method="euclidean", plot.title="Euclidean - CLR - Food PCs", env.vars=colnames(food.pc))    
    
### Intra-inter group variabilities
    plot.within.group.distances(map0=map[cs,], otu0=taxa_L7[cs,], fn="within.group.euc.pdf", ylab="Euclidean distance")
    plot.between.group.distances(map0=map[cs,], otu0=taxa_L7[cs,], fn="between.group.euc.pdf", ylab="Euclidean distance")

### Relative Distance - Cross-Sectional Plots
    xlab="Years in US"
    x.var="Years.in.US"

    # all 1st gen to all thai
    plot.relative.CS(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7, main="Euclidean", ref_samples=c(karenthai,hmongthai), xlab=xlab, x.var=x.var, ylab="Distance to all Thai", outputfn="allfirstgen_to_Thai.pdf")

    # all to controls
    plot.relative.CS(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7, main="Euclidean", ref_samples=c(controls), xlab=xlab, x.var=x.var, ylab="Distance to Controls", outputfn="allfirstgen_to_Controls.pdf")

    # 1stKaren KarenThai
    plot.relative.CS(map0=map[karen_firstgen_cs,], otu0=taxa_L7, main="Euclidean", ref_samples=karenthai, xlab=xlab, x.var=x.var, ylab="Distance to Karen Thai", outputfn="Karen_CS_to_Thai.pdf")

    # 1stKaren to 2ndGenHmong
    plot.relative.CS(map0=map[karen_firstgen_cs,], otu0=taxa_L7, main="Euclidean", ref_samples=hmong_secondgen_cs, xlab=xlab, x.var=x.var, ylab="Distance to 2nd-gen", outputfn="Karen_CS_to_USborn.pdf")

    # 1stHmong to HmongThai
    plot.relative.CS(map0=map[hmong_firstgen_cs,], otu0=taxa_L7, main="Euclidean", ref_samples=hmongthai, xlab=xlab, x.var=x.var, ylab="Distance to Hmong Thai", outputfn="Hmong_CS_to_Thai.pdf")

    # 1stHmong to 2ndHmong
    plot.relative.CS(map0=map[hmong_firstgen_cs,], otu0=taxa_L7, main="Euclidean", ref_samples=hmong_secondgen_cs, xlab=xlab, x.var=x.var, ylab="Distance to US born Hmong", outputfn="Hmong_CS_to_USbornHmong.pdf")

    # 1stHmong to Controls
    plot.relative.CS(map0=map[hmong_firstgen_cs,], otu0=taxa_L7, main="Euclidean", ref_samples=controls, xlab=xlab, x.var=x.var, ylab="Distance to Controls", outputfn="Hmong_CS_to_Controls.pdf")

    plot.relative.CS(map0=map[karen_firstgen_cs,], otu0=taxa_L7, main="Euclidean", ref_samples=controls, xlab=xlab, x.var=x.var, ylab="Distance to Controls", outputfn="Karen_CS_to_Controls.pdf")

    # TODO fix this to work for euc distance 
    # Lean vs. Obese only
    # Hypothesis: are obese individuals more different from their groups?
    # 1. Hmong1st and Hmong2nd to HmongThai    
    # 2. all US to all Thai
    # 3. Karen1st to KarenThai
    mains <- c("Unweighted Unifrac","Weighted Unifrac", "Bray-Curtis")
    dms <- list(bc_dm)
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

### Differential Taxa
    # Lean vs Obese
    # create factor version of years in us so that data isn't omitted for NAs or 0s - to be used only when dealing with the entire dataset
    lean.obese.cs.map <- map[map$BMI.Class %in% c("Lean", "Obese") & (is.na(map$Sample.Order) | map$Sample.Order==1),]
    Years.in.US.Factor <- lean.obese.cs.map$Years.in.US
    Years.in.US.Factor[Years.in.US.Factor == 0] <- 50 # 2nd gen are 0, let's set them to something higher like 50??
    Years.in.US.Factor[is.na(Years.in.US.Factor)] <- 0 # Thai are NA, let's set them to 0
    Years.in.US.Factor <- factor(cut(Years.in.US.Factor,c(-1, .0001, seq(5,35,5), 41, 50)), ordered=T)
    lean.obese.cs.map[,"Years.in.US.Factor"] <- Years.in.US.Factor

    # Add ethnicity as a control
    plot.diff.taxa(lean.obese.cs.map, taxa_L6, x.var="BMI.Class", 
        control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L6", sig.level=.10, do.sqrt=FALSE, do.filter=TRUE)

    plot.diff.taxa(lean.obese.cs.map, taxa_L2, x.var="BMI.Class", 
        control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend="BMI.all.L2", sig.level=.1, do.sqrt=FALSE, do.filter=TRUE)

    plot.diff.taxa(lean.obese.cs.map, taxa_L6, x.var="Waist.Height.Ratio", 
        control.vars=c("Age","Years.in.US.Factor","Ethnicity"), outputfn.prepend="WHR.all.L6", sig.level=.10, do.sqrt=FALSE, do.filter=TRUE)

    ###### Not sure how useful this is?????? check code and params
        #iterate through all combos of looking at diff taxa in lean vs. obese by subgroup
        groups <- c("HmongThai","Hmong1st","Hmong2nd","KarenThai","Karen1st")
        taxatables <- list(taxa_L6, taxa_L7)
        names(taxatables) <- c("L6", "L7")
        combos <- expand.grid(names(taxatables),groups, stringsAsFactors=F)
        colnames(combos) <- c("taxa","group")
        combos[combos$group %in% c("KarenThai","HmongThai", "Hmong2nd"), "controls"] <- "Age"
        combos[combos$group %in% c("Karen1st","Hmong1st"), "controls"] <- c("Age,Years.in.US")
        for(i in 1:nrow(combos))
        {
            combo.name <- paste("BMI", paste(unlist(combos[i,c("group","taxa")]), collapse="."), sep=".")
            print(combo.name) # so we know which iteration we are in
            plot.diff.taxa(lean.obese.cs.map[lean.obese.cs.map$Sample.Group == combos$group[i],], taxatables[[combos$taxa[i]]], x.var="BMI.Class", 
                control.vars=unlist(strsplit(combos$controls[i],",")), outputfn.prepend=combo.name, sig.level=sig.level)
        }

### Heatmap
    make.heatmap.traditional(otu0=taxa, map0=map[c(hmongthai, hmong_secondgen_cs),], outputfn="heatmap.trad.taxa.hmong.pdf", filter.mode="clr", min.prev=.25, save.pvals=TRUE)
    pvals <- read.table("heatmap.trad.taxa.hmong.toppvals.txt", sep="\t", colClasses="character")[,1]
    make.heatmap.traditional(otu0=taxa[,pvals], map0=map[c(hmongthai, hmong_firstgen_cs, hmong_secondgen_cs),], outputfn="heatmap.trad.taxa.hmong.all.pdf", filter.mode="clr")

# ------ Longitudinal CLR ----

### PCOA
    plot.pcoa.long(map, samples=rownames(map[map$Sub.Study=="L",]), otu0 = taxa_L7, method="euclidean", plot.title="Euclidean - CLR - L.outlined", convex.hull=TRUE, flip.axis=1)


### Relative Distance - Longitudinal Plots
    map_L <- map[map$Sub.Study == "L",]
    # IMP.049 and IMP.050 were recruited at 2 months and 3 months respectively - these might be throwing things off
    map_L <- map_L[!(map_L$Subject.ID %in%  c("IMP.049", "IMP.050")),]

    # by sample order seems to work better
        # L to self
        xlab="Sample Number"
        x.var="Sample.Order"
        main="Euclidean"
        p1 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main=main, ref.sample.order=1, xlab=xlab, 
                        ylab="Dissimilarity to First Sample", x.var=x.var) + theme(legend.position='none')
        p2 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main=main, ref.sample.order=6, xlab=xlab, 
                        ylab="Dissimilarity to Last Sample", x.var=x.var) + theme(legend.position='none')
        p3 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main=main, ref.sample.order=1, to.previous=T, xlab=xlab, 
                        ylab="Dissimilarity to Previous Sample", x.var=x.var)
        save_plot("Karen_L_to_Selfs_Sample.Order.pdf", plot_grid(p1,p2,p3, nrow=1), base_aspect_ratio = 3)
    

        # L to KarenThai
        p4 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main="Euclidean", ref_samples=karenthai, xlab=xlab, 
                        ylab="Distance to Karen in Thailand", x.var=x.var) + theme(legend.position='none')
        # L to Controls
        p5 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main="Euclidean", ref_samples=controls, xlab=xlab, 
                        ylab="Distance to Controls", x.var=x.var)
        save_plot("Karen_L_to_Groups_Sample.Order.pdf", plot_grid(p4,p5, nrow=1), base_aspect_ratio = 1.3)


    # by sample day since arrival
        # L to self
        xlab="Days Since US Arrival"
        x.var="Sample.Day.Since.Arrival"
        main="Euclidean"
        p1 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main=main, ref.sample.order=1, xlab=xlab, 
                        ylab="Dissimilarity to First Sample", x.var=x.var) + theme(legend.position='none')
        p2 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main=main, ref.sample.order=6, xlab=xlab, 
                        ylab="Dissimilarity to Last Sample", x.var=x.var) + theme(legend.position='none')
        p3 <- plot.relative.L(map0 = map_L, otu0=taxa_L7, main=main, ref.sample.order=1, to.previous=T, xlab=xlab, 
                        ylab="Dissimilarity to Previous Sample", x.var=x.var)
        save_plot("Karen_L_to_Selfs_Sample.Day.pdf", list(p1,p2,p3), ncol = 3, nrow = 1, base_aspect_ratio = 1.3)
    

        # L to KarenThai
        p4 <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_L7, main="Euclidean", ref_samples=karenthai, xlab=xlab, 
                        ylab="Distance to Karen in Thailand", x.var=x.var) + theme(legend.position='none')
        # L to Controls
        p5 <- plot.relative.L(map0 = map[l.samples,], otu0=taxa_L7, main="Euclidean", ref_samples=controls, xlab=xlab, 
                        ylab="Distance to Controls", x.var=x.var)
        save_plot("Karen_L_to_Groups_Sample.Day.pdf", list(p4,p5), ncol = 2, nrow = 1, base_aspect_ratio = 1.3)

    # paired test of M1 vs M6 distances to: 
        # KarenThai
        m1 <- rownames(map)[map$Sample.Group=="Karen1st" & map$Sub.Study=="L" & map$Sample.Month == 1]    
        m6 <- rownames(map)[map$Sample.Group=="Karen1st" & map$Sub.Study=="L" & map$Sample.Month == 6]    
        ret <- prep.dm(map0=map[c(m1,m6,karenthai),], otu0=taxa_L7[c(m1,m6,karenthai),], dm=NULL, method="euclidean")
        rdist <- get.relative.distance(query_samples=c(m1,m6), ref_samples=karenthai, dm=ret$dm)
        
        d <- data.frame(y=rdist, month=as.factor(c(rep(1,length(m1)), rep(6,length(m6)))))
        p <- ggplot(d, aes(month, y)) + geom_quasirandom(dodge.width=.75) + geom_boxplot(alpha=0, colour="black") +
            ylab("Distance to Thai") + xlab("Month")  
        p <- ggdraw(p) + draw_figure_label(label=paste0("p=",signif(t.test(rdist[m1], rdist[m6], paired=T)$p.value,2)), size=8, position="top.right")
        save_plot("Karen_M1M6_to_Thai.pdf", p, base_aspect_ratio=1.3)

        # Controls
        ret <- prep.dm(map0=map[c(m1,m6,controls),], otu0=taxa_L7[c(m1,m6,controls),], dm=NULL, method="euclidean")
        rdist <- get.relative.distance(query_samples=c(m1,m6), ref_samples=controls, dm=ret$dm)        
        d <- data.frame(y=rdist, month=as.factor(c(rep(1,length(m1)), rep(6,length(m6)))))
        p <- ggplot(d, aes(month, y)) + geom_quasirandom(dodge.width=.75) + geom_boxplot(alpha=0, colour="black") +
            ylab("Distance to Controls") + xlab("Month")  
        p <- ggdraw(p) + draw_figure_label(label=paste0("p=",signif(t.test(rdist[m1], rdist[m6], paired=T)$p.value,2)), size=8, position="top.right")
        save_plot("Karen_M1M6_to_Controls.pdf", p, base_aspect_ratio=1.3)

    
### DIFF TAXA
    # diff genus between M1 vs M6
    s1 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 1,]))
    s2 <- sort(rownames(map[map$Sub.Study=="L" & map$Sample.Order == 6,]))
    taxa0 <- taxa[c(s1,s2),]
    prevalences <- apply(taxa0, 2, function(bug.col) mean(bug.col > 0))
    taxa0 <- taxa0[, prevalences >= .10]
    ret <- collapse.by.correlation(taxa0, .95)
    taxa0 <- taxa0[, ret$reps]

    ret <- test.features.paired(taxa0, samples1=s1, samples2=s2, sig.level=.10, parametric=TRUE)
    if(length(ret$features) > 0)
        for(i in 1:length(ret$features)){  
            ggdata <- data.frame(x=map[c(s1,s2),"Sample.Order"], y=taxa0[,ret$features[i]])
            p[[i]] <- ggplot(ggdata, aes(x,y)) + geom_boxplot()
        }
        
    