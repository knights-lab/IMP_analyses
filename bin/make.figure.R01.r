# R01 stuff only

# figure 1
    p1a <- plot.pcoa(map[cs,], otu0=taxa_clr_L7[cs,], method="euclidean", plot.title="Microbiome Westernization", env.vars=c("Years.in.US"), flip.axis=2, show.stats=FALSE)    

    # unweighted doesn't look as nice
    #p1b <- plot.pcoa(map_all[cs_all,], dm=food_uwuf_dm, plot.title="Unweighted - Dietary Westernization", show.stats=FALSE, flip.axis=1)    

    # weighted
    p1b <- plot.pcoa(map_all[cs_all,], dm=food_wuf_dm, plot.title="Dietary Westernization", show.stats=FALSE, env.vars=c("Years.in.US"))
    
    leg <- get_legend(p1b)
    lhs <- plot_grid(p1a + theme(legend.position="none"), p1b + theme(legend.position="none"), labels=c("A","B"), ncol=1)
    p <- plot_grid(lhs, leg, rel_widths=c(3,1))
    # for some reason manual saving looks better
#    save_plot("Figure 1 - PCOA.pdf", p, useDingbats=FALSE, base_aspect_ratio = .95 )

    # get only samples we have microbiome samples sequenced for
    this.cs <- cs[cs != "IMP000.L.T30"]
    map.cs <- map[this.cs,]    
    ret <- prep.dm(map0=map.cs, otu0=taxa_clr_L7, method="euclidean")
    this.dm <- ret$dm
    mb.distance <- get.relative.distance(query_samples=this.cs, ref_samples=c(karenthai,hmongthai), dm=this.dm)
    map.cs$mb.dist <- mb.distance
    food.distance <- get.relative.distance(query_samples=this.cs, ref_samples=c(karenthai_all,hmongthai_all), dm=food_wuf_dm)
    map.cs$food.dist <- food.distance
    this.bp <- log10(taxa[this.cs, bacteroides]/taxa[this.cs, prevotella])
    
    plot(food.distance, this.bp)
    cor.test(this.bp, food.distance, method="spear")    
        # 	Spearman's rank correlation rho
        # 
        # data:  this.bp and food.distance
        # S = 3444200, p-value = 1.637e-10
        # alternative hypothesis: true rho is not equal to 0
        # sample estimates:
        #       rho 
        # 0.3512581 
    summary(lm(mb.dist ~ food.dist + BMI + Age + Years.in.US + Ethnicity, data=map.cs))
        # 
        # Call:
        # lm(formula = mb.dist ~ food.dist + BMI + Age + Years.in.US + 
        #     Ethnicity, data = map.cs)
        # 
        # Residuals:
        #     Min      1Q  Median      3Q     Max 
        # -8.1055 -2.4050 -0.6694  2.0593 15.9363 
        # 
        # Coefficients:
        #                Estimate Std. Error t value Pr(>|t|)    
        # (Intercept)    35.19867    1.69232  20.799  < 2e-16 ***
        # food.dist       0.23660    0.13462   1.758 0.079809 .  
        # BMI            -0.05430    0.04662  -1.165 0.245026    
        # Age             0.05696    0.01816   3.136 0.001876 ** 
        # Years.in.US     0.08854    0.01625   5.447 1.04e-07 ***
        # EthnicityKaren -1.84500    0.48604  -3.796 0.000177 ***
        # ---
        # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
        # 
        # Residual standard error: 3.719 on 311 degrees of freedom
        # Multiple R-squared:  0.279,	Adjusted R-squared:  0.2675 
        # F-statistic: 24.07 on 5 and 311 DF,  p-value: < 2.2e-16
    lm.food <- lm(mb.dist ~ food.dist + BMI + Age + Years.in.US + Ethnicity, data=map.cs)
    lm.nofood <- lm(mb.dist ~ BMI + Age + Years.in.US + Ethnicity, data=map.cs)
    anova(lm.food, lm.nofood)
        # Analysis of Variance Table
        # 
        # Model 1: mb.dist ~ food.dist + BMI + Age + Years.in.US + Ethnicity
        # Model 2: mb.dist ~ BMI + Age + Years.in.US + Ethnicity
        #   Res.Df    RSS Df Sum of Sq     F  Pr(>F)  
        # 1    311 4300.7                             
        # 2    312 4343.5 -1   -42.717 3.089 0.07981 .

    # none of these variables seem to really matter
#    summary(lm(mb.dist ~ food.dist + BMI + Age + Years.in.US + Ethnicity + Type.Birth.Location + Type.location.before.US +
 #               Alcohol.Use + Tobacco.Use + Breastfed, data=map.cs))

    summary(lm(food.dist ~ Sample.Group, data=map.cs)) #--> R2 33% 
    summary(lm(mb.dist ~ Sample.Group, data=map.cs)) # --> R2 19%
    
    plot.constrained.ordination(map.cs, otu0=taxa_clr_L7[this.cs,], method="euclidean", plot.title="CLR - CS by WUF Food Distance", env.vars=c("food.dist"))    

    plot.constrained.ordination(map.cs, otu0=taxa_clr_L7[this.cs,], method="euclidean", plot.title="CLR - CS by Yrs AND WUF Food Distance", env.vars=c("Years.in.US", "food.dist"))    

    plot.constrained.ordination(map.cs, otu0=taxa_clr_L7[this.cs,], plot.title="RDA constrained by multiple vars", env.vars=c("Fraction.Life.in.US", "food.dist", "BMI", "Age"))    
    

## Fig 2 
    #####  make sure to use RARE alpha tables
    p2a <- plot.boxplot.by.group.x.bmi(map00=lean.obese.cs.map, y=alphadiv0[rownames(lean.obese.cs.map), "PD_whole_tree"], ylab="Faith's Phylogenetic Diversity", main="Alpha Diversity", parametric=TRUE, add.stats=F)
    query_samples_list <- list(c(hmong_secondgen_cs, hmong_firstgen_cs), c(hmong_secondgen_cs,hmong_firstgen_cs, karen_firstgen_cs), karen_firstgen_cs)
    ref_samples_list <- list(hmongthai, c(hmongthai, karenthai), karenthai)
    fn_ext <- c("Hmong", "All", "Karen")
    i = 2
        map_temp <- map[query_samples_list[[i]],]
        map_temp <- map_temp[map_temp$BMI.Class %in% c("Lean","Obese"),]
        new_query_samples <- rownames(map_temp)
        # calculate the relative distances to reference groups as response variable
        all.samples <- c(new_query_samples, ref_samples_list[[i]])
        ddm <- vegdist(taxa_clr_L7[all.samples,], method="euc")
        dm <- as.matrix(ddm)
        rel.dist <- get.relative.distance(new_query_samples, ref_samples_list[[i]], dm)
    p2b <- plot.boxplot.by.group.x.bmi(map00=map_temp, y=rel.dist, ylab="Distance to Counterparts in Thailand", main="Microbiome Dissimilarity", add.stats=F)
    p <- plot_grid(p2a + theme(legend.position="none"), p2b, labels=c("A","B"), ncol=2)
    p <- ggdraw(p) + draw_label("Figure 2", x = 1, y = 0, hjust=1.5, vjust=-.5, size = 12)
    save_plot("Figure 2.pdf", p, base_aspect_ratio = 2)
      
## Fig 3 
    plot.b.p.ratio.all.2 <- function(map0, otu, bug1, bug2, outputfn, g1, g2, g3)
    {
        xvar <-  "Years.in.US" #"Fraction.Life.in.US" #
        map0 <- map0[c(g1,g2,g3),]
        otu <- otu[rownames(map0),]
    
        bp <- get.taxa.ratio(otu, bug1, bug2)
        
        pch <- get.pch(map0)

        d <- data.frame(x = map0[g1, xvar], y=-1*log10(bp[g1]))
        d2 <- data.frame(x = map0[g2, xvar], y=-1*log10(bp[g2]), group=map0[g2, "Sample.Group"])
        d3 <- data.frame(x = map0[g3, xvar], y=-1*log10(bp[g3]), group=map0[g3, "Sample.Group"])
        
        if(xvar == "Years.in.US"){
            d2 <-data.frame(x=rep(-2, length(g2)), y=-1*log10(bp[g2]), group=map0[g2, "Sample.Group"])
            d3 <-data.frame(x=rep(43, length(g3)), y=-1*log10(bp[g3]), group=map0[g3, "Sample.Group"])
        }

        dd <- rbind(d2, d3)
        dd$group <- factor(dd$group, levels=sort(as.character(unique(dd$group))))

        cols <- get.group.colors(as.character(unique(dd$group)), .5)

        main <- ggplot() + geom_point(data=d, aes(x=x,y=y), color=alpha("black",.3)) + 
                geom_smooth(data=d, aes(x=x,y=y), color="black", size=.8) #+ xlim(-2,43)
       
        main <- main + geom_point(data=dd, aes(x=x,y=y, color=group)) +
            scale_color_manual(values=cols) +
            ggtitle("Bacteroides-Prevotella Ratio") + # title
            ylab("log10(Bacteroides-Prevotella Ratio)") + xlab(xvar) + 
            geom_hline(yintercept=0, linetype="dashed", color="black", size=.5) +
            theme(legend.title=element_blank(), legend.text=element_text(size=10))

#        p <- ggdraw(main) + draw_label("Figure 3", x = 1, y = 0, hjust=1.5, vjust=-.5, size = 12)
        save_plot(outputfn, main, useDingbats=FALSE, base_aspect_ratio = 1.3 )
    }
    # flip the order of the bugs, and the function should multiply by -1 to flip it back
    plot.b.p.ratio.all.2(map, taxa, bug2=bacteroides, bug1=prevotella, outputfn="b.p.ratio.all.pdf", g1=firstgen_cs, g2=c(hmongthai,karenthai),g3=(hmong_secondgen_cs))
  
   
