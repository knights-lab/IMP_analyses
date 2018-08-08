# change these next two lines to run on your local machine
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

datadirname <- "denovo"
source("bin/load.r")

plot.map<- function(database,center,...){
    Obj <- map(database,...,plot=F)
    coord <- cbind(Obj[[1]],Obj[[2]])

    # split up the coordinates
    id <- rle(!is.na(coord[,1]))
    id <- matrix(c(1,cumsum(id$lengths)),ncol=2,byrow=T)
    polygons <- apply(id,1,function(i){coord[i[1]:i[2],]})

    # split up polygons that differ too much
    polygons <- lapply(polygons,function(x){
        x[,1] <- x[,1] + center
        x[,1] <- ifelse(x[,1]>180,x[,1]-360,x[,1])
        if(sum(diff(x[,1])>300,na.rm=T) >0){
          id <- x[,1] < 0
          x <- rbind(x[id,],c(NA,NA),x[!id,])
       }
       x
    })
    # reconstruct the object
    polygons <- do.call(rbind,polygons)
    Obj[[1]] <- polygons[,1]
    Obj[[2]] <- polygons[,2]

    map(Obj,...)
}

save.split <- function(p0, fn, skip.legend=FALSE)
{
    if(!skip.legend)
    {
        leg <- get_legend(p0)
        p0 <- p0 + theme(legend.position='none')
        save_plot(paste0("output_manuscript/",fn,"-legend.pdf"), leg, useDingbats=FALSE, base_aspect_ratio = 1)
    }
    save_plot(paste0("output_manuscript/",fn,".pdf"), p0, useDingbats=FALSE, base_aspect_ratio = 1)
}

### STUDY DESIGN
pdf("worldmap.pdf", width=10); plot.map("world", center=170, col="gray",bg="white",border=NA,
        fill=TRUE,ylim=c(-60,90),mar=c(0,0,0,0)); dev.off()

### PCOA
    
    ### PCOA with boxplots
        p <- plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="", show.stats=F, save.pc=T)
        leg <- get_legend(p)
        p <- p + theme(legend.position='none', axis.ticks=element_blank(), axis.text=element_blank()) + 
            scale_x_continuous(position="top") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), 
            axis.title = element_text(margin = unit(c(0,0,0,0),"cm")))

        pc <- read.table("-PC.txt", sep="\t", skip=1, row=1)
        colnames(pc) <- paste0("PC",1:2)
    
        # add boxplots to the axis
        ggdata <- data.frame(pc[cs,], Sample.Group=map[cs,"Sample.Group"])

        cols <- get.group.colors()
        cols["Hmong2nd"] <- "black"
        fills <- get.group.colors()
        alphas <- get.group.alphas() 
        alphas[c("Karen1st","Hmong1st")] <- 0
        xlabels <- SAMPLE.GROUP.NAMES.SHORT[levels(ggdata$Sample.Group)]

        pc1 <- ggplot(ggdata, aes(x=Sample.Group, y=PC1)) + geom_boxplot(aes(color=Sample.Group, fill=Sample.Group, alpha=Sample.Group)) + coord_flip() + 
                scale_y_continuous(position = "top") + scale_x_discrete(labels=xlabels, position="top", limits=rev(levels(ggdata$Sample.Group))) +
                scale_fill_manual(name="Groups",values=fills) +
                scale_color_manual(name="Groups",values=cols) +
                scale_alpha_manual(name="Groups",values=alphas) + 
                theme(legend.position='none', axis.title=element_blank(), axis.text.y=element_blank()) +
                theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

        # top, right, bottom, left

        pc2 <- ggplot(ggdata, aes(x=Sample.Group, y=PC2)) + geom_boxplot(aes(color=Sample.Group, fill=Sample.Group, alpha=Sample.Group)) + 
                scale_x_discrete(labels=xlabels) +
                scale_fill_manual(name="Groups",values=fills) +
                scale_color_manual(name="Groups",values=cols) +
                scale_alpha_manual(name="Groups",values=alphas) + 
                theme(legend.position='none', axis.title=element_blank(), axis.text.x=element_blank()) + 
                 theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
                 
        # can't get the alignment of the boxplots to work out perfectly, paste these together manually
        prow1 <- plot_grid(p, pc2, ncol=2, nrow=1, rel_widths=c(4.25,1), align="h")
        save_plot("pcoa-boxplots1.pdf", prow1, useDingbats=FALSE, base_aspect_ratio = 1, base_height=5, base_width=6)        

        prow2 <- plot_grid(p, pc1, ncol=1, nrow=2, rel_heights=c(5,1), align="v")
        save_plot("pcoa-boxplots2.pdf", prow2, useDingbats=FALSE, base_aspect_ratio = 1, base_width=5, base_height=6)        
    
    
        pc1.blank <- ggplot(ggdata, aes(x=Sample.Group, y=PC1)) + geom_boxplot(color="white", fill="white") + 
                theme(legend.position='none', axis.title=element_blank(), axis.text=element_blank()) + 
                 theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), line = element_blank())

        pc2.blank <- ggplot(ggdata, aes(x=Sample.Group, y=PC2)) + geom_boxplot(color="white", fill="white") + 
                theme(legend.position='none', axis.title=element_blank(), axis.text=element_blank()) + 
                 theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), line = element_blank())

        # make sequential plots of PCOA for presentation purposes
        p1 <- plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="", show.stats=F, save.pc=F, hide.groups=c("Control", "Hmong2nd", "Hmong1st", "Karen1st"))
        p1 <- p1 + theme(legend.position='none', axis.ticks=element_blank(), axis.text=element_blank()) + 
            scale_x_continuous(position="top") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(margin = unit(c(0,0,0,0),"cm")))
        p1grid <- plot_grid(p1, pc2.blank, ncol=2, nrow=1, rel_widths=c(4.25,1), align="h")
        save_plot("pcoa-seq1.pdf", p1grid, useDingbats=FALSE, base_aspect_ratio = 1, base_height=5, base_width=6)        

        p2 <- plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="", show.stats=F, save.pc=F, hide.groups=c("Hmong2nd", "Hmong1st", "Karen1st"))
        p2 <- p2 + theme(legend.position='none', axis.ticks=element_blank(), axis.text=element_blank()) + 
            scale_x_continuous(position="top") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(margin = unit(c(0,0,0,0),"cm")))
        p2grid <- plot_grid(p2, pc2.blank, ncol=2, nrow=1, rel_widths=c(4.25,1), align="h")
        save_plot("pcoa-seq2.pdf", p2grid, useDingbats=FALSE, base_aspect_ratio = 1, base_height=5, base_width=6)        

        p3 <- plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="", show.stats=F, save.pc=F, hide.groups=c("Hmong1st", "Karen1st"))
        p3 <- p3 + theme(legend.position='none', axis.ticks=element_blank(), axis.text=element_blank()) + 
            scale_x_continuous(position="top") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.title = element_text(margin = unit(c(0,0,0,0),"cm")))
        p3grid <- plot_grid(p3, pc2.blank, ncol=2, nrow=1, rel_widths=c(4.25,1), align="h")
        save_plot("pcoa-seq3.pdf", p3grid, useDingbats=FALSE, base_aspect_ratio = 1, base_height=5, base_width=6)        
            

    ### PCOA with diversity, years
        p <- plot.pcoa.by(map0=cbind(map[cs,], alphadiv[cs, "PD_whole_tree",drop=F]), dm=uwuf_dm, fill.by="PD_whole_tree")
        save.split(p, "pcoa-PD_whole_tree")
        cor.test(map[c(hmong_firstgen_cs,karen_firstgen_cs),"Years.in.US"],alphadiv[c(hmong_firstgen_cs,karen_firstgen_cs), "PD_whole_tree"], method="spear")

        maptemp <- map[c(hmong_firstgen_cs, karen_firstgen_cs),] # first generation only
        p <- plot.pcoa.by(maptemp, dm=uwuf_dm, fill.by="Years.in.US", fill.color="#452D81", low.fill.color="#F4D968")
        save.split(p, "pcoa-years.in.us")

    ### Taxa heatmap    
    make.heatmap.binary(otu0=otu, map0=map[c(hmongthai,hmong_firstgen_cs),], min.prevalence=0.75, baseline.groups="HmongThai", show.colnames=T, sig.level=1, outputfn="heatmap.diffotu.HT-H1.binary.p75.pdf",
                        is.otu=T, taxamapfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/ALL_RUNS/denovo/rep_set.tax")


    ### alpha diversity
        map.obese.lean <- map[cs,]
        map.obese.lean <- map.obese.lean[map.obese.lean$BMI.Class %in% c("Lean","Obese"),]

        multiplot.boxplot.by.group.x.bmi(map00=map.obese.lean, y.list = as.list(alphadiv[rownames(map.obese.lean), c("shannon","PD_whole_tree")]), 
                    ylabs=rep("",2), 
                    mains=c( "Shannon index",paste0("Faith", "'", "s PD")),
                     outputfn="MANUSCRIPT/boxplot.alphadiv.bmi.pdf", add.pvals=T)

        # suppl figure
        k = map.obese.lean[map.obese.lean$Ethnicity=="Karen",]
        multiplot.boxplot.by.group.x.bmi(map00=k, y.list = as.list(alphadiv[rownames(k), c("shannon","PD_whole_tree")]), 
                    ylabs=rep("",2), 
                    mains=c( "shannon index",paste0("Faith", "'", "s PD")),
                     outputfn="MANUSCRIPT/boxplot.alphadiv.karen.obeselean.pdf", add.pvals=T)

        h = map.obese.lean[map.obese.lean$Ethnicity=="Hmong",]
        multiplot.boxplot.by.group.x.bmi(map00=h, y.list = as.list(alphadiv[rownames(h), c("shannon","PD_whole_tree")]), 
                    ylabs=rep("",2), 
                    mains=c( "shannon index",paste0("Faith", "'", "s PD")),
                     outputfn="MANUSCRIPT/boxplot.alphadiv.hmong.obeselean.pdf", add.pvals=T)

        # two-way anova lean vs obese across all groups
        ggdata <- cbind(map[cs,],alphadiv=alphadiv[cs, "PD_whole_tree"])
        ggdata <- ggdata[ggdata$BMI.Class %in% c("Lean","Obese"),]
        f <- as.formula("alphadiv ~ Sample.Group + BMI.Class")
        a <- Anova(mod=aov(f, data=ggdata), type = "III") # this for unbalanced classes!
        print(a)

        # alpha diversity x years, lines fitted by BMI class
        ggdata <- cbind(map[c(hmong_firstgen_cs,karen_firstgen_cs),],alphadiv=alphadiv[c(hmong_firstgen_cs,karen_firstgen_cs), "PD_whole_tree"])

        ggdata <- ggdata[ggdata$BMI.Class %in% c("Lean","Obese"),]
        ggdata$BMI.Class <- factor(ggdata$BMI.Class, levels=rev(levels(ggdata$BMI.Class)))


        p <- ggplot(data=ggdata,
                aes(y=Years.in.US,x=alphadiv,group=BMI.Class, color=BMI.Class)) + geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                xlab("Faith's PD") + ylab("Years in US")

        save.split(p, "alphadiv.bmi.firstgen", T)


        # all first
            this.map <- map[c(karen_firstgen_cs, hmong_firstgen_cs),]
            this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
            f <- as.formula("alphadiv ~ Years.in.US + BMI + Ethnicity")
            summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
            # Years.in.US    -0.17530    0.08543  -2.052  0.04112 *  
            # BMI            -0.16352    0.17041  -0.960  0.33811    
            # EthnicityKaren  6.45244    2.06078   3.131  0.00193 ** 
        
            p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv)) + 
                    geom_point(alpha=.5, aes(color=BMI.Class)) + geom_smooth(method="lm",alpha=.2, color="black") + 
                    scale_color_manual(name="", values=get.bmi.colors(alpha=1)) +
                    ylab("Faith's PD") + xlab("Years in US") + ggtitle("") + facet_grid(. ~ Ethnicity, scales="free")
            save_plot(paste0("MANUSCRIPT/","alphadiv.allfirstgen.oneline",".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 2)
        # hmong
            this.map <- map[hmong_firstgen_cs,]
            this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
            f <- as.formula("alphadiv ~ Years.in.US + BMI")
            summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
            # Years.in.US -0.18912    0.08483  -2.229   0.0275 *  
            # BMI         -0.05370    0.22491  -0.239   0.8116    
       
            # with interaction       
            # Years.in.US      0.37270    0.45592   0.817    0.415    
            # BMI              0.43315    0.44842   0.966    0.336    
            # Years.in.US:BMI -0.02105    0.01678  -1.254    0.212
 
            p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv)) + 
                    geom_point(alpha=.5, aes(color=BMI.Class)) + geom_smooth(method="lm",alpha=.2, color="black") + 
                    scale_color_manual(name="", values=get.bmi.colors(alpha=1)) +
                    ylab("Faith's PD") + xlab("Years in US") + ggtitle("Hmong") + ylim(25,100)
            save.split(p, "alphadiv.hmongfirstgen.oneline")

        #karen
            this.map <- map[karen_firstgen_cs,]
            this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
            f <- as.formula("alphadiv ~ Years.in.US + BMI")
            summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
            # Years.in.US  -0.1713     0.4524  -0.379    0.706    
            # BMI          -0.2793     0.2593  -1.077    0.283    
            p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv)) + 
                    geom_point(alpha=.5, aes(color=BMI.Class)) + geom_smooth(method="lm",alpha=.2, color="black") + 
                    scale_color_manual(name="", values=get.bmi.colors(alpha=1)) +
                    ylab("Faith's PD") + xlab("Years in US") + ggtitle("Karen")+ ylim(25,100)
            save.split(p, "alphadiv.karenfirstgen.oneline")


### euclidean relative distances
        datadirname <- "denovo.clr"
        LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
        setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
        source("bin/load.r")

        # all 1st gen to all thai
        xlab="Years in US"
        x.var="Years.in.US"
        p <- plot.relative.CS(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7, ref_samples=c(karenthai,hmongthai), xlab=xlab, x.var=x.var, main="", ylab="Similarity to Thai", show.stats=F, similarity=T)
        save.split(p, "first-to-thai",T)

        # all to controls
        p<-plot.relative.CS(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7, ref_samples=c(controls), xlab=xlab, x.var=x.var, main="", ylab="Similarity to US", show.stats=F, similarity=T)
        save.split(p, "first-to-controls",T)

        p <- plot.relative.CS.facet(map0=map[c(hmong_firstgen_cs, karen_firstgen_cs),], otu0=taxa_L7, ref_samples1=c(karenthai,hmongthai), ref_samples2=controls, xlab=xlab, x.var=x.var, similarity=T)
        save_plot("relative-distance.pdf",p, ncol=2, base_aspect_ratio=1)


### BP
    
    ### BP all
        plot.b.p.ratio.all(map[cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.pdf", do.insets=TRUE)
    
    ### BP - all as boxplots
    plot.b.p.ratio.all.boxplots(map[cs,], taxa,bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.boxplots.pdf")

    ### BP with Thai groups as one, no color
    plot.b.p.ratio.all.nocolor(map[cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.nocolor.pdf")

    ### BP with Thai groups as one, with color
    plot.b.p.ratio.all.onethai(map[cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.all.onethai.pdf")

    # supp
    plot.b.p.ratio.all(map[karen_firstgen_cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.karen1st.pdf", do.insets=F)

### FOOD
    # load food data
    food_wuf_dm <- read.table(paste(food_datadir, "wuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_uwuf_dm <- read.table(paste(food_datadir,"uwuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)

    # nutrient file should already be mapped to sample ids
    nutrients <- read.table(nutrientsfn, sep="\t", header=T, check.names=F, as.is=T,row=1)
    nutrients[,"% of Calories from Total Sugars"] <- nutrients[,"Total Sugars in Grams"] * 4 / nutrients[,"Total Calories"] * 100
    nutrients[,"g Fiber per 1000 Calories"] <- (nutrients[,"Dietary Fiber in Grams",] / nutrients[,"Total Calories"]) * 1000

    ### food pcoa
    p <- plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="", show.stats=F)
    save.split(p, "food-pcoa-wuf")

    p <- plot.pcoa(map[cs,], dm=food_uwuf_dm, plot.title="", show.stats=F)
    save.split(p, "food-pcoa-uwuf")

    
    ### RDA MB by Food PCs
    # generate food PCs here (weighted unifrac)
    plot.pcoa(map[cs,], dm=food_uwuf_dm, plot.title="Food Unweighted Unifrac", show.stats=FALSE, save.pc=TRUE, axis1 = 1, axis2 = 5) # put axis 5 so we get all 5 PCs back
    # load in food PCs
    food.pc <- read.table("Food Unweighted Unifrac-PC.txt", sep="\t", skip=1, row=1)
    colnames(food.pc) <- paste0("Food.PC", 1:5)
    # redo constrained ordination with food PCs only
    p <- plot.constrained.ordination(data.frame(map[cs,], food.pc[cs,]), dm0=uwuf_dm, plot.title="Unweighted Unifrac - Food PCs", env.vars=colnames(food.pc))    
    save.split(p0=p, fn="RDA-UWUF-x-foodPCs",T)


    nutrient.vars <- c("Total Calories", "% of Calories from Total Sugars", "% of Calories from Added Sugars","% of Calories from Carbohydrate","% of Calories from Protein",
    "% of Calories from Saturated Fat", "% of Calories from Total Fat","g Fiber per 1000 Calories")

    # two-way anova
    p.2way <- lapply(nutrient.vars, function(nut)  map.boxplot(y=nutrients[cs,nut], Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=F, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=12, x.size=9, show.x=T, group.vars.df=map[cs,c("Ethnicity","Birth.Continent","Resident.Continent")]) + 
                                                theme(plot.title = element_text(size=18)))

    p.pval.2way <- lapply(nutrient.vars, function(nut)  map.boxplot(y=nutrients[cs,nut], Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=T, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=12, x.size=9, show.x=T, group.vars.df=map[cs,c("Ethnicity","Birth.Continent","Resident.Continent")]) + 
                                                theme(plot.title = element_text(size=18)) )
    save_plot("MANUSCRIPT/macronutrients2.pdf", plot_grid(plotlist=p.2way[c(1,2,7,5)], nrow=2, ncol=2), nrow=2, ncol=2, base_aspect_ratio=1,  base_height=5)
    save_plot("MANUSCRIPT/macronutrients2.pval.pdf", plot_grid(plotlist=p.pval.2way[c(1,2,7,5)], nrow=2, ncol=2), nrow=2, ncol=2, base_aspect_ratio=1,  base_height=5)

    p <- lapply(nutrient.vars, function(nut)  map.boxplot(y=nutrients[cs,nut], Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=F, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=12, x.size=9, show.x=T, group.vars.df=NULL) + 
                                                theme(plot.title = element_text(size=18)))

    p.pval <- lapply(nutrient.vars, function(nut)  map.boxplot(y=nutrients[cs,nut], Group=map[cs,"Sample.Group"], main=nut, facet.var=NULL, alpha=.01, add.pval=T, 
                                                plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=12, x.size=9, show.x=T, group.vars.df=NULL) + 
                                                theme(plot.title = element_text(size=18)) )

    save_plot("MANUSCRIPT/macronutrients.pdf", plot_grid(plotlist=p[c(1,2,7,5)], nrow=2, ncol=2), nrow=2, ncol=2, base_aspect_ratio=1,  base_height=5)
    save_plot("MANUSCRIPT/macronutrients.pval.pdf", plot_grid(plotlist=p.pval[c(1,2,7,5)], nrow=2, ncol=2), nrow=2, ncol=2, base_aspect_ratio=1,  base_height=5)

    
### LONGITUDINAL
    map_L <- map[map$Sub.Study == "L",]
    # IMP.049 and IMP.050 were recruited at 2 months and 3 months respectively and samples came late - drop these out
    map_L <- map_L[!(map_L$Subject.ID %in%  c("IMP.049", "IMP.050")),]

    ## longitudinal with all other samples
        p <- plot.pcoa.long(map, samples=rownames(map_L), dm=uwuf_dm, plot.title="", return.pc=F)				
        save.split(p, "UnWeighted Unifrac - L")
    
        # longitudinal only
        p<-plot.pcoa.long(map, samples=rownames(map_L), dm=uwuf_dm[rownames(map_L),rownames(map_L)], plot.title="", return.pc=F)				
        save.split(p, "UnWeighted Unifrac - L-ONLY")

    # plot the PC barplots
        pc <- plot.pcoa.long(map_L, dm=uwuf_dm[rownames(map_L),rownames(map_L)], plot.title="", return.pc=T)				
        p <- plot.pcoa.long(map_L, dm=uwuf_dm, plot.title="", return.pc=F)				
        save.split(p, "pcoa.UWUF.L.samples", F)
        
        pc.start <- pc[pc$Type=="start",]
        pc.start <- pc.start[order(pc.start$subject),]
        pc.end <- pc[pc$Type=="end",]
        pc.end <- pc.end[order(pc.end$subject),]    
    
        d <- rbind(data.frame(delta = pc.end$x - pc.start$x, PC="PC1", subject=pc.start$subject),
                    data.frame(delta = pc.end$y - pc.start$y, PC="PC2", subject=pc.start$subject))
    
        # t.test is fine, data is normal here
        p1 <- paste0("P=", signif(t.test(d[d$PC=="PC1","delta"])$p.value, 2))
        p2 <- paste0("P=", signif(t.test(d[d$PC=="PC2","delta"])$p.value, 2))

        pbox <- ggplot(d, aes(x=PC,y=delta)) + geom_boxplot(lwd=1) + theme_bw() + ylim(c(-.2, .2)) +
                theme(axis.title=element_blank(), axis.text.x=element_text(size=20), axis.text.y=element_text(size=15),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_rect(colour = "black", fill=NA, size=2)) + annotate("text", x=1, y=.2, label=p1, size=8) + annotate("text", x=2, y=.2, label=p2, size=8)
        save.split(pbox, "pcoa.UWUF.Lsamples.boxplots", T)

    # plot with all samples 
        pc <- plot.pcoa.long(map, dm=uwuf_dm, plot.title="", return.pc=T)				
        p <- plot.pcoa.long(map, dm=uwuf_dm, plot.title="", return.pc=F)				
        save.split(p, "pcoa.UWUF.ALL.samples", F)

        pc.start <- pc[pc$Type=="start",]
        pc.start <- pc.start[order(pc.start$subject),]
        pc.end <- pc[pc$Type=="end",]
        pc.end <- pc.end[order(pc.end$subject),]    
    
        d <- rbind(data.frame(delta = pc.end$x - pc.start$x, PC="PC1", subject=pc.start$subject),
                    data.frame(delta = pc.end$y - pc.start$y, PC="PC2", subject=pc.start$subject))
    
        # t.test is fine, data is normal here
        p1 <- paste0("P=", signif(t.test(d[d$PC=="PC1","delta"])$p.value, 2))
        p2 <- paste0("P=", signif(t.test(d[d$PC=="PC2","delta"])$p.value, 2))

        pbox <- ggplot(d, aes(x=PC,y=delta)) + geom_boxplot(lwd=1) + theme_bw() + ylim(c(-.2, .2)) +
                theme(axis.title=element_blank(), axis.text.x=element_text(size=20), axis.text.y=element_text(size=15),
                                    panel.grid.major = element_blank(), 
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_rect(colour = "black", fill=NA, size=2)) + annotate("text", x=1, y=.2, label=p1, size=8) + annotate("text", x=2, y=.2, label=p2, size=8)
        save.split(pbox, "pcoa.UWUF.ALL.samples.boxplots", T)


    # taxa bar plots of n=6 subjects
        camp.subjects <- map_L[map_L$Recruitment.Location=="Maela Camp","Subject.ID"]
        plot.taxa.summary.L(taxa0=taxa, map=map_L[map_L$Subject.ID %in% camp.subjects,], outputfn="taxa.L.pdf", 
                        max.taxa=15, x.var="Sample.Day.Since.Arrival", grid.ncol=3)

        plot.taxa.summary.L(taxa0=food_otu_L3, map=map_L[map_L$Subject.ID %in% camp.subjects,], outputfn="foodtaxa.L.pdf", 
                        max.taxa=15, x.var="Sample.Day.Since.Arrival", grid.ncol=3)

    # plot BMI change
        map_L <- map[map$Sub.Study == "L",]
        mapL_M6 <- map_L[map_L$Sample.Month==6,"BMI.M6",drop=F]
        colnames(mapL_M6) <- "BMI"
        mapL_M1 <- unique(map_L[order(map_L$Sample.Order),"BMI",drop=F])
        BMIL <- rbind(mapL_M6, mapL_M1)
        p <- plot.response.L(map_L[rownames(BMIL),], y=BMIL$BMI, outputfn=NULL, ylab="(kg/m2)", ggtitle="BMI", num.clip.months=1, show.stats=F)        
        save.split(p$p, "BMI.L")

    # Food Alpha Diversity - L
        map_L <- map[map$Sub.Study == "L" & map$Sample.Order %in% 0:6,] # there's one person with a 7th month sample, but we'll drop it
        p <- plot.response.L(map_L, y=food_alphadiv[rownames(map_L), "PD_whole_tree"], outputfn=NULL, ggtitle="Food Diversity", ylab="Faith's PD", num.clip.months=1, show.stats=F)
        save.split(p$p, "food.diversity")
        
    # calories - L
        mac.names <- c("Total Calories", "% of Calories from Carbohydrate", "% of Calories from Total Fat", "g Fiber per 1000 Calories", "% of Calories from Total Sugars", "% of Calories from Protein")

        p <- list()
        for(mac in mac.names)
        {
            p[[mac]] <- plot.response.L(map_L, y=nutrients[rownames(map_L), mac], outputfn=NULL, ylab="", ggtitle=mac, num.clip.months=1, show.stats=F)
        }
        save.split(p[["% of Calories from Protein"]]$p, "protein.L") 

        pvals <- lapply(p, function(xx) xx$pval)
        adj.pvals <- p.adjust(pvals, method="fdr")
        names(adj.pvals) <- mac.names
        pp <- list()
        for(mac in mac.names) 
        {
            pp[[mac]] <- p[[mac]]$p
            pp[[mac]] <- pp[[mac]] + ggtitle(paste0("P=",adj.pvals[mac])) + theme(plot.title = element_text(face="plain"))
        }
        # protein and carbs are sig after adjustment
        save.split(pp[["% of Calories from Protein"]]$p, "protein.L") 

    # BP ratio
        map_L <- map[map$Sub.Study == "L",]
        # IMP.049 and IMP.050 were recruited at 2 months and 3 months respectively - these might be throwing things off
        map_L <- map_L[!(map_L$Subject.ID %in%  c("IMP.049", "IMP.050")),]
        p <- plot.b.p.ratio.L(map_L, taxa, bug1=bacteroides, bug2=prevotella, outputfn="b.p.ratio.L.first.last.pdf", num.clip.months=1, show.stats=F)
        save.split(p$p, "b.p.ratio.L")


### Shotgun BP Coverage and Abundance bubble plot
    datadirname <- "shotgun"
    LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
    setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
    source("bin/load.r")

    # plot CS and Karen at the same time, but save as separate files - with sample names
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                rescale="none", outputfn="bubble.BP.cov50.per.person.xlab.pdf", max.features=NULL, min.coverage=.5, per="person", 
                baseheight=4.5, basewidth=13, show.samplenames=TRUE)

    # plot CS and Karen at the same time, but save as separate files - withOUT sample names
    b.p.heatmap(mapfile=mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bcov.unique.binary.txt", 
                rescale="none", outputfn="bubble.BP.cov50.per.person.pdf", max.features=NULL, min.coverage=.5, per="person", 
                baseheight=4, basewidth=13)
    
### MOUSE
    LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
    setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
    source(paste0(LIBDIR,"mouse.load.data.r")) 

    ## PCOA
        map_pcoa_nobase1 <- map_mb[map_mb$Week %in% c(8,10) & map_mb$Donor %in% c("TFSCS023","IMP.263"),]
        plot.mouse.pcoa(map_pcoa_nobase1, dm=wuf_dm, outputfn="pcoa.pair1.lastpoints.wuf.pdf", type="allpoints", verbose=TRUE)

    ## IMMUNE
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
        
        selected.cells <- c("TCRgd", "CD8aaTCRab", "CD4TCRab")
        
        # Significant in both LPL and IEL (per CD45+): TCRgd, CD8aaTCRab, CD4TCRab
        p.iel.pval <- make.immune.plots.manuscript(ggdata=map[map$timepoint=="endpoint" & map$Date != max(map$Date) & map$Cohoused==FALSE,], 
                immune=immune.live, immune_apc=immune_apc, 
                group.var="Group.End", grouping.vars=c("Donor.Type","Diet.Type"),
                add.pval=TRUE, outdir="immune-cells-all-live/", selected.cells=selected.cells, cell.pop="IEL")

        save_plot("IEL-pval.pdf", plot_grid(plotlist=p.iel.pval, align="h", ncol=4), ncol=4, base_aspect_ratio=1)

        p.iel <- make.immune.plots.manuscript(ggdata=map[map$timepoint=="endpoint" & map$Date != max(map$Date) & map$Cohoused==FALSE,], 
                immune=immune.live, immune_apc=immune_apc, 
                group.var="Group.End", grouping.vars=c("Donor.Type","Diet.Type"),
                add.pval=FALSE, outdir="immune-cells-all-live/", selected.cells=selected.cells, cell.pop="IEL")
        titles <- paste0(c("CD45+",selected.cells),"\n")
        for(i in 1:length(p.iel))
            p.iel[[i]] <- p.iel[[i]] + ggtitle(titles[i])

        save_plot("IEL.pdf", plot_grid(plotlist=p.iel, align="h", ncol=4), ncol=4, base_aspect_ratio=1)

    ## BODY COMPOSITION
        colnames(body_comp) <- gsub("\\.\\.g\\.", "", colnames(body_comp))
        body_comp <- read.table(file=body_comp_fn, sep="\t", header=T, as.is=T)    
        gg_bodycomp <- merge(body_comp[,c("Mouse.ID","Percent.Fat","Percent.Lean")], map[map$timepoint=="endpoint",], by="Mouse.ID")
        gg_bodycomp <- gg_bodycomp[gg_bodycomp$Cohoused==FALSE,]
        p_fat <- mouse.boxplot(y=gg_bodycomp$Percent.Fat, Group=gg_bodycomp$Group.End, 
                main="Percent Fat\n", facet.var=NULL, add.pval=FALSE,
                ylab="", strip.text.size=10,group.vars.df=gg_bodycomp[,c("Donor.Type","Diet.Type")], hide.box=F)
        save_plot("mouse-percent-fat.pdf", p_fat, base_aspect_ratio = 1)
           
        p_fat <- mouse.boxplot(y=gg_bodycomp$Percent.Fat, Group=gg_bodycomp$Group.End, 
                main="Percent Fat", facet.var=NULL, add.pval=TRUE,
                ylab="", strip.text.size=10,group.vars.df=gg_bodycomp[,c("Donor.Type","Diet.Type")], hide.box=FALSE)
        save_plot("mouse-percent-fat-pval.pdf", p_fat, base_aspect_ratio = 1)
                
    ## Feed Efficiency
        map4groups <- map[map$Cohoused==FALSE,]
        plot.feed.efficiency.wk8(map4groups, group="Group.End", add.pval=TRUE, outputfn="feed-efficiency-endpoint-pval.pdf")                
        plot.feed.efficiency.wk8(map4groups, group="Group.End", add.pval=FALSE, outputfn="feed-efficiency-endpoint.pdf")                
                
    ## Blood Glucose
        map_glucose <- map[map$timepoint=="endpoint" & map$Date != max(map$Date) & map$Cohoused==FALSE,]
        plot.glucose(map[!is.na(map$Fasting.Glucose) & map$Cohoused==F,], group="Group.End", add.pval=TRUE, outputfn="mouse-glucose-foldchange-pval.pdf")              
        plot.glucose(map[!is.na(map$Fasting.Glucose) & map$Cohoused==F,], group="Group.End", add.pval=FALSE, outputfn="mouse-glucose-foldchange.pdf")              
                
                
    ## Villus Crypt
        plot.vc(merge(villus_crypt, map[map$timepoint=="endpoint",], by="Mouse.ID"), add.pval=FALSE, outputfn="villus-crypt-ratio.pdf")
        plot.v(merge(villus_crypt, map[map$timepoint=="endpoint",], by="Mouse.ID"), add.pval=FALSE, outputfn="villus.pdf")        
       
                
                
                
                
                
                
                