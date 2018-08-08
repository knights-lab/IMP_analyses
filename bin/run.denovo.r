# change these next two lines to run on your local machine
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

datadirname <- "denovo"
source("bin/load.r")

# ----- CROSS-SECTIONAL PLOTS ----- 

### alpha diversity
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(alphadiv[cs, alpha.metrics]), ylabs=rep("",length(alpha.metrics)), mains=alpha.metrics, outputfn="boxplot.alphadiv.bmi.pdf")

    am <- c("PD_whole_tree","chao1","shannon")
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(alphadiv[cs, am]), ylabs=rep("",length(am)), mains=am, outputfn="boxplot.alphadiv.bmi.2.pdf")

    # add pvals to plot
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(alphadiv[cs, alpha.metrics]), ylabs=rep("",length(alpha.metrics)), mains=alpha.metrics, outputfn="boxplot.alphadiv.bmi.PVAL.pdf", add.pvals=TRUE)
 
    # try alpha diversity with taxa table
    alphadiv.taxa <- read.table(paste(datadir,"alpha.taxa.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")
    multiplot.boxplot.by.group.x.bmi(map00=map[cs,], y.list = as.list(alphadiv.taxa[cs, alpha.nonphylo]), ylabs=rep("",length(alpha.nonphylo)), mains=alpha.nonphylo, outputfn="boxplot.alphadiv_taxa.bmi.pdf")

    # look at diversity correlations with Age AND Years.in.US
    ph <- plot.XY(data.frame(map[hmong_firstgen_cs,], diversity=alphadiv[hmong_firstgen_cs,"PD_whole_tree"])[,c("Years.in.US", "BMI", "diversity")], z.symmetric=F)
    
    pk <- plot.XY(data.frame(map[karen_firstgen_cs,], diversity=alphadiv[karen_firstgen_cs,"PD_whole_tree"])[,c("Years.in.US", "BMI", "diversity")], z.symmetric=F)

    phk <- plot.XY(data.frame(map[c(hmong_firstgen_cs,karen_firstgen_cs),], diversity=alphadiv[c(hmong_firstgen_cs,karen_firstgen_cs),"PD_whole_tree"])[,c("Years.in.US", "diversity", "BMI")], z.symmetric=F)
    
    save_plot("diversity-x-Years.Age-Hmong1st.pdf", ph, base_aspect_ratio = 1.3)
    save_plot("diversity-x-Years.Age-Karen1st.pdf", pk, base_aspect_ratio = 1.3)
    
    # run some stats to look at associations in first-gen only
        # hmong < 10 years
        this.map <- map[hmong_firstgen_cs,]
        this.map <- this.map[this.map$Years.in.US < 10,]
        this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
        f <- as.formula("alphadiv ~ Years.in.US + BMI")
        summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
        # Years.in.US  -1.2323     0.6549  -1.882   0.0687 .  
        # BMI           0.3236     0.4687   0.690   0.4948    
        
        p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv)) + 
                geom_point(alpha=.5, aes(color=BMI.Class)) + geom_smooth(method="lm",alpha=.2, color="black") + 
                scale_color_manual(name="", values=get.bmi.colors(alpha=1)) +
                ylab("Faith's PD") + xlab("Years in US") + ggtitle("Hmong") + ylim(25,100)
        save.split(p, "alphadiv.hmongfirstgen10.oneline")

        this.map <- map[map$Sample.Group %in% c("Hmong1st","Karen1st") & map$BMI.Class %in% c("Lean","Obese"),]
        f <- as.formula("alphadiv ~ Years.in.US + BMI.Class + Ethnicity")
        summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
        # Years.in.US     -0.2078     0.1116  -1.861  0.06386 .  
        # BMI.ClassObese  -2.9582     2.1021  -1.407  0.16059    
        # EthnicityKaren   7.4231     2.6544   2.797  0.00556 ** 

        f <- as.formula("alphadiv ~ Years.in.US + BMI.Class + Ethnicity + Years.in.US:BMI.Class + Years.in.US:Ethnicity + BMI.Class:Ethnicity")
        summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
        #         Years.in.US                   -0.15478    0.14906  -1.038  0.30008    
        #         BMI.ClassObese                 2.07547    7.15747   0.290  0.77208    
        #         EthnicityKaren                 8.49567    3.17900   2.672  0.00803 ** 
        #         Years.in.US:BMI.ClassObese    -0.18263    0.25668  -0.712  0.47742    
        #         Years.in.US:EthnicityKaren    -0.05536    0.46681  -0.119  0.90570    
        #         BMI.ClassObese:EthnicityKaren -5.00159    7.23748  -0.691  0.49016  

        f <- as.formula("alphadiv ~ Years.in.US * BMI.Class * Ethnicity")
        summary(lm(f, data=cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])))
        #             Years.in.US                                -0.1343     0.1503  -0.894  0.37236    
        #             BMI.ClassObese                              3.6349     7.3106   0.497  0.61949    
        #             EthnicityKaren                              9.1238     3.2350   2.820  0.00518 ** 
        #             Years.in.US:BMI.ClassObese                 -0.2455     0.2636  -0.931  0.35256    
        #             Years.in.US:EthnicityKaren                 -0.3039     0.5240  -0.580  0.56251    
        #             BMI.ClassObese:EthnicityKaren              -8.9276     8.1562  -1.095  0.27476    
        #             Years.in.US:BMI.ClassObese:EthnicityKaren   1.2026     1.1527   1.043  0.29784          

        f <- as.formula("alphadiv ~ Years.in.US + Ethnicity")
        this.map <- map[map$Sample.Group %in% c("Hmong1st","Karen1st") & map$BMI.Class=="Obese",]
        this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
        summary(lm(f, data=this.map))
        # Years.in.US     -0.3399     0.1892  -1.797   0.0781 .  
        # EthnicityKaren   3.3027     5.8919   0.561   0.5775    
        p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                ylab("Faith's PD") + xlab("Years in US") + ggtitle("Obese")

        this.map <- map[map$Sample.Group %in% c("Hmong1st","Karen1st") & map$BMI.Class=="Lean",]
        this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
        summary(lm(f, data=this.map))
        # Years.in.US     -0.1593     0.1475  -1.080  0.28137    
        # EthnicityKaren   8.3557     3.0234   2.764  0.00625 ** 
        p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                ylab("Faith's PD") + xlab("Years in US") + ggtitle("Obese")
        save.split(p, "alphadiv.bmi.karenfirstgen", T)
    
    
        # Hmong Obese
            f <- as.formula("alphadiv ~ Years.in.US")
            this.map <- map[map$Sample.Group %in% c("Hmong1st") & map$BMI.Class=="Obese",]
            this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
            summary(lm(f, data=this.map))
            # Years.in.US  -0.3799     0.1786  -2.127   0.0421 *  
            p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                    geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                    scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                    ylab("Faith's PD") + xlab("Years in US") + ggtitle("Obese")
        # Hmong Overweight
            f <- as.formula("alphadiv ~ Years.in.US")
            this.map <- map[map$Sample.Group %in% c("Hmong1st") & map$BMI.Class=="Overweight",]
            this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
            summary(lm(f, data=this.map))
            #Years.in.US -0.09546    0.12295  -0.776    0.441    
            p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                    geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                    scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                    ylab("Faith's PD") + xlab("Years in US") + ggtitle("Obese")
        # Hmong Lean
                this.map <- map[map$Sample.Group %in% c("Hmong1st") & map$BMI.Class=="Lean",]
                this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
                summary(lm(f, data=this.map))
                # Years.in.US  -0.1343     0.1539  -0.873    0.387    
                p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                        geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                        scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                        ylab("Faith's PD") + xlab("Years in US") + ggtitle("Lean") 
        # Karen Obese
                f <- as.formula("alphadiv ~ Years.in.US")
                this.map <- map[map$Sample.Group %in% c("Karen1st") & map$BMI.Class=="Obese",]
                this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
                summary(lm(f, data=this.map))
                # Years.in.US   0.5188     0.9780   0.531    0.601    
                p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                        geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                        scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                        ylab("Faith's PD") + xlab("Years in US") + ggtitle("Obese")
        # Karen Lean
                this.map <- map[map$Sample.Group %in% c("Karen1st") & map$BMI.Class=="Lean",]
                this.map <- cbind(this.map, alphadiv=alphadiv[rownames(this.map), "PD_whole_tree"])
                summary(lm(f, data=this.map))
                # Years.in.US  -0.4382     0.5155   -0.85    0.397    
                p <- ggplot(data=this.map,aes(x=Years.in.US,y=alphadiv,color=BMI.Class)) + 
                        geom_point(alpha=.5) + geom_smooth(method="lm",alpha=.2) + 
                        scale_color_manual(name="", values=get.bmi.colors(alpha=1)) + theme(legend.position="none") +
                        ylab("Faith's PD") + xlab("Years in US") + ggtitle("Lean")
    
    
### body trends

    p <- plot.BMI.barplot(map[cs,], bins=seq(0,45,5), freq=F)
    save_plot("BMI.pdf", p+theme(legend.position="none"), base_aspect_ratio=1)

    p <- plot.BMI.barplot(map[c(karenthai, karen_firstgen_cs),], bins=seq(0,10,1), freq=F)    
    save_plot("BMI.Karen.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)

    p <- plot.BMI.barplot(map[c(hmongthai, hmong_firstgen_cs, hmong_secondgen_cs),], bins=seq(0,45,5),freq=F)
    save_plot("BMI.Hmong.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)

    p <- plot.BMI.barplot(map[c(karenthai, karen_firstgen_cs),], bins=seq(0,10,1), freq=T)    
    save_plot("BMI.Karen.hist.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)

    p <- plot.BMI.barplot(map[c(hmongthai, hmong_firstgen_cs, hmong_secondgen_cs),], bins=seq(0,45,5),freq=T)
    save_plot("BMI.Hmong.hist.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)

    leg <- get_legend(p)
    save_plot("BMI.legend.pdf", leg, base_aspect_ratio=1)

    # try overweight/obese only
    maph <- map[c(hmongthai, hmong_firstgen_cs, hmong_secondgen_cs),]
    p <- plot.BMI.barplot(maph[maph$BMI.Class %in% c("Overweight","Obese"),], bins=seq(0,45,5),freq=F,main="Hmong")
    save_plot("BMI.Hmong.OO.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)
    p <- plot.BMI.barplot(maph[maph$BMI.Class %in% c("Overweight","Obese"),], bins=seq(0,45,5),freq=T)
    save_plot("BMI.Hmong.OO.hist.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)

    mapk <- map[c(karenthai,karen_firstgen_cs),]
    p <- plot.BMI.barplot(mapk[mapk$BMI.Class %in% c("Overweight","Obese"),], bins=seq(0,10,1), freq=F,main="Karen")    
    save_plot("BMI.Karen.OO.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)
    p <- plot.BMI.barplot(mapk[mapk$BMI.Class %in% c("Overweight","Obese"),], bins=seq(0,10,1), freq=T)    
    save_plot("BMI.Karen.OO.hist.pdf", p + theme(legend.position="none"), base_aspect_ratio=1)


    
    plot.WHR.boxplot(map[c(karenthai,karen_firstgen_cs),], bins=0:10, fn="WHR_boxplot_Karen.pdf")   
    plot.WHR.boxplot(map[c(hmongthai,hmong_firstgen_cs,hmong_secondgen_cs),], bins=seq(0,45,5), fn="WHR_boxplot_Hmong.pdf")
    
### taxa summaries
    plot.taxa.summary(map0=map[cs,], otu=taxa, fn="taxa.summary.pdf")    

### PCOA 
    p <- plot.pcoa(map[cs,], dm=wuf_dm, plot.title="Weighted Unifrac", flip.axis=1:2)
    save_plot("pcoa - Weighted Unifrac.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

    p <- plot.pcoa(map[cs,], dm=bc_dm, plot.title="Bray Curtis")
    save_plot("pcoa - Bray Curtis.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

    p <- plot.pcoa(map[cs,], dm=uwuf_dm, plot.title="Unweighted Unifrac")
    save_plot("pcoa - Unweighted Unifrac.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
    
    
    # plot PCOA of firstgen only and color and fill by different metadata vars
    maptemp <- map[c(hmong_firstgen_cs, karen_firstgen_cs),]
    maptemp$Arrived.As <- ifelse(maptemp$Age.at.Arrival >= 18, "Adult", "Child")
    p <- plot.pcoa.by(map0=maptemp, dm=bc_dm, fn=, fill.by="Years.in.US", shape.by="Arrived.As")
    save_plot("pcoa-BC-Years-x-AgeArrived.Binary.pdf", p, base_aspect_ratio=1.3)
    
    # plot PCOA by alpha diversity
    p <- plot.pcoa.by(map0=cbind(map[cs,], alphadiv[cs, "PD_whole_tree",drop=F]), dm=bc_dm, fill.by="PD_whole_tree")
    save_plot("pcoa-BC-All-PD_whole_tree.pdf", p, base_aspect_ratio=1.3)

    pdiv <- plot.pcoa.by(map0=cbind(maptemp, alphadiv[rownames(maptemp), "PD_whole_tree",drop=F]), dm=bc_dm, fill.by="PD_whole_tree")
    pdiv.leg <- get_legend(pdiv)
    pdiv <- pdiv + theme(legend.position='none')
    pyrs <- plot.pcoa.by(map0=maptemp, dm=bc_dm, fill.by="Years.in.US", fill.color="#00441B")
    pyrs.leg <- get_legend(pyrs)
    pyrs <- pyrs + theme(legend.position='none')
    
    multiplot <- plot_grid(pdiv, pyrs, plot_grid(pdiv.leg, pyrs.leg, ncol=1, align="v"), ncol=3, rel_widths=c(1,1,.5))
    save_plot("pcoa-BC-PD_whole_tree,Yrs.in.US.pdf", multiplot, ncol = 2, nrow = 1, base_aspect_ratio = 1.3)

    dm <- uwuf_dm
    pdiv <- plot.pcoa.by(map0=cbind(maptemp, alphadiv[rownames(maptemp), "PD_whole_tree",drop=F]), dm=dm, fill.by="PD_whole_tree")
    pdiv.leg <- get_legend(pdiv)
    pdiv <- pdiv + theme(legend.position='none')
    pyrs <- plot.pcoa.by(map0=maptemp, dm=dm, fill.by="Years.in.US", fill.color="#00441B")
    pyrs.leg <- get_legend(pyrs)
    pyrs <- pyrs + theme(legend.position='none')
    
    multiplot <- plot_grid(pdiv, pyrs, plot_grid(pdiv.leg, pyrs.leg, ncol=1, align="v"), ncol=3, rel_widths=c(1,1,.5))
    save_plot("pcoa-UWUF-PD_whole_tree,Yrs.in.US.pdf", multiplot, ncol = 2, nrow = 1, base_aspect_ratio = 1.3)

    # RDA: constrain ordination with selected env variables using Unweighted Unifrac DM
    #plot.constrained.ordination(map[cs,], dm0=uwuf_dm, plot.title="Unweighted Unifrac", env.vars=c("Years.in.US","BMI","Age"))    
    #plot.constrained.ordination(map[cs,], dm0=uwuf_dm, plot.title="Unweighted Unifrac - Full Model")    

    # generate food PCs here (unweighted unifrac)
    plot.pcoa(map[cs,], dm=food_uwuf_dm, plot.title="Food Unweighted Unifrac", show.stats=FALSE, save.pc=TRUE, axis1 = 1, axis2 = 20) # put axis 5 so we get all 5 PCs back
    # load in food PCs
    food.pc <- read.table("Food Unweighted Unifrac-PC.txt", sep="\t", skip=1, row=1)
    colnames(food.pc) <- paste0("Food.PC", 1:ncol(food.pc))

    # redo constrained ordination with food PCs only
    p <- plot.constrained.ordination(data.frame(map[cs,], food.pc[cs,1:10]), dm0=uwuf_dm, plot.title="RDA - Unweighted Unifrac - Food PCs", env.vars=colnames(food.pc)[1:10])    
    save_plot("RDA-UWUF.MB-UWUF.Food.PC.pdf",p, base_aspect_ratio=1.3)

    # generate food PCs here (weighted unifrac)
    plot.pcoa(map[cs,], dm=food_wuf_dm, plot.title="Food Weighted Unifrac", show.stats=FALSE, save.pc=TRUE, axis1 = 1, axis2 = 5) # put axis 5 so we get all 5 PCs back
    # load in food PCs
    food.pc <- read.table("Food Weighted Unifrac-PC.txt", sep="\t", skip=1, row=1)
    colnames(food.pc) <- paste0("Food.PC", 1:ncol(food.pc))

    # redo constrained ordination with food PCs only
    p <- plot.constrained.ordination(data.frame(map[cs,], food.pc[cs,]), dm0=wuf_dm, plot.title="RDA - Weighted Unifrac - Food PCs", env.vars=colnames(food.pc))    
    save_plot("RDA-WUF.MB-WUF.Food.PC.pdf",p, base_aspect_ratio=1.3)
    
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

    make.heatmap.binary(otu0=otu, map0=map[c(hmongthai,hmong_firstgen_cs),], min.prevalence=0.75, baseline.groups="HmongThai", show.colnames=F, outputfn="heatmap.otu.hmongthai.binary.pdf")

    # let differentiated OTUs drive what's shown, but still apply a small prevalence filter


#     make.heatmap.binary(otu0=otu, map0=map[c(hmongthai,hmong_firstgen_cs),], min.prevalence=0.75, baseline.groups="HmongThai", show.colnames=T, sig.level=.01, outputfn="heatmap.diffotu.HT-H1.binary.p75.s01.pdf",
#                         is.otu=T, taxamapfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/ALL_RUNS/denovo/rep_set.tax")

    # don't filter based on differentiation, instead just collect the qvalues so we can refer to them later
    make.heatmap.binary(otu0=otu, map0=map[c(hmongthai,hmong_firstgen_cs),], min.prevalence=0.75, baseline.groups="HmongThai", show.colnames=T, sig.level=1, outputfn="heatmap.diffotu.HT-H1.binary.p75.pdf",
                        is.otu=T, taxamapfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/ALL_RUNS/denovo/rep_set.tax")

    # check when H1 is baseline
    #    make.heatmap.binary(otu0=otu, map0=map[c(hmongthai,hmong_firstgen_cs),], min.prevalence=0.75, baseline.groups="Hmong1st", show.colnames=F, sig.level=.10, outputfn="heatmap.diffotu.H1-HT.binary.pdf")


# ----- LONGITUDINAL PLOTS ----- 
    
    map_L <- map[map$Sub.Study == "L",]
    # IMP.049 and IMP.050 were recruited at 2 months and 3 months respectively - these might be throwing things off
    map_L <- map_L[!(map_L$Subject.ID %in%  c("IMP.049", "IMP.050")),]

### alpha diversity, one per subject
    #multiplot.alphadiv.L(map_L, alphadiv, alpha.metrics, outputfn="alphadiv.L.Num.pdf", x.var="Sample.Order", xlab="Sample Number") # by sample number
    multiplot.alphadiv.L(map_L, alphadiv, alpha.metrics, outputfn="alphadiv.L.Month.pdf") # defaults to Sample Order
    multiplot.alphadiv.L(map_L, alphadiv, alpha.metrics, outputfn="alphadiv.L.first.last.pdf", num.clip.months=1)
    multiplot.alphadiv.L(map_L, alphadiv, alpha.metrics, outputfn="alphadiv.L.first2.last2.pdf", num.clip.months=2)
    
### taxa summaries
    # stream plots are useful for longitudinal changes over time
    plot.taxa.summary.L(taxa0=taxa, map0=map[map$Sub.Study=="L", ], outputfn="taxa.summary.L.pdf", grid.ncol=2)

### PCOA - consider using CLR
    puw <- plot.pcoa.long(map, samples=rownames(map_L), dm=uwuf_dm, plot.title="UnWeighted Unifrac - L")
    pw <- plot.pcoa.long(map, samples=rownames(map_L), dm=wuf_dm, plot.title="Weighted Unifrac - L", flip.axis=1)

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
