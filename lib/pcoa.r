require(ggpubr) 
require(vegan)
require(faraway)
require(RColorBrewer)
require(cowplot)

# label.samples = vector of samples to label by name
# used for plotting samples for selecting Hmong Mouse Donors
plot.pcoa.by <- function(map0, dm, ethnicity, fn, color.by="Years.in.US", and.by="BMI.Class", label.samples=NULL)
{
    map0 <- map0[map0$Ethnicity %in% ethnicity,]

    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    dm <- dm[valid_samples, valid_samples]
    
    ddm <- as.dist(dm)
    # plot clusters
    # get PCoA scores first
    pc0 <- cmdscale(ddm,2)

    pch <- rep(19,nrow(map0))
    pch[map0$BMI.Class=="Obese"] <- 17 #mark obese as triangles
    pch[map0$BMI.Class=="Lean"] <- 15 #mark lean as squares
    
    d <- data.frame(x = pc0[,1], y = pc0[,2], 
                    color.by=map0[,color.by], group=map0$Sample.Group, label=rownames(map0))
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(color=color.by), size=2, pch=pch) + xlab("PC1") + ylab("PC2") +
                scale_color_gradient(low="yellow",high="red")

    if(!is.null(label.samples))
        p <- p + geom_text(aes(label=ifelse(label %in% label.samples, as.character(label), '')), hjust=0, vjust=0)
        
#        (pc0[label.samples,1], pc0[label.samples,2],labels=label.samples, cex=.5)
    
    save_plot(fn, p,
        base_aspect_ratio = 1.3
    )
}

# dm is optional and only for things that have been generated outside of R (unifrac)
# env.vars = specify continuous vars to plot additional environment variables 
# this function constrains the ordination by a list of environment variables passed in
plot.constrained.ordination <- function(map0, otu0, method="euclidean", plot.title, dm0=NULL, env.vars=NULL)
{   
    ret <- prep.dm(map0, otu0, dm0, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0

    # RDA: distance is euclidean (CLR transformed) and relationship with env.vars is expected to be linear
    if(is.null(dm0) & !is.null(env.vars))
        rrda <- rda(as.formula(paste("otu0 ~ ", paste(env.vars, collapse=" + "))), map0) # basic RDA with euc distance 
    else if(is.null(dm0) & is.null(env.vars))
        rrda <- rda(otu0) # basic RDA with euc distance - full model
    else if(!is.null(dm0) & !is.null(env.vars))
        rrda <- dbrda(as.formula(paste("ddm ~ ", paste(env.vars, collapse=" + "))), map0) # distance based RDA
    else if(!is.null(dm0) & is.null(env.vars))
        rrda <- dbrda(ddm) # distance based RDA - full model
    
    rrda.plot <- plot(rrda)
    
    pc <- scores(rrda)$sites
    
    d <- data.frame(x = pc[,1], y = pc[,2], group=map0$Sample.Group)
    # set the levels of Sample.Group so that it's the same every time
    d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))
    group.cols <- alpha(c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"),.8)
    if(length(unique(d$group)) > 5)
        group.cols <- c(alpha("black", .5), group.cols)

    percent_var <- signif(eigenvals(rrda)/sum(eigenvals(rrda)), 4)*100
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=group), size=2) +
        xlab(paste0("Axis 1 [",percent_var[1],"%]")) +
        ylab(paste0("Axis 2 [",percent_var[2],"%]")) + 
        scale_color_manual(values=group.cols) + #sets the color palette of the fill
        stat_ellipse(data=d, aes(colour=group), show.legend=F, type="t", level=.6)

    if(!is.null(env.vars))
    {
        print(anova(rrda, by="terms")) # significant terms here contribute most to constrained model
    
        # total variation explained by environment vars
        label <- paste0("Variation explained: ", signif(rrda$CCA$tot.chi/rrda$tot.chi, 4)*100, "%")
    
        vector.multiplier <- attributes(rrda.plot$biplot)$arrow.mul

        fit_env_df <- as.data.frame(rrda.plot$biplot*vector.multiplier)
        fit_env_df <- cbind(fit_env_df, Variable = rownames(rrda.plot$biplot))
        colnames(fit_env_df)[1:2] <- c("Dim1","Dim2")

        text_fit_env_df <- fit_env_df
        text_fit_env_df[,2] <- text_fit_env_df[,2]*1.2 # let's shift the text out a little
    
        p <- p + coord_fixed() + ## need aspect ratio of 1!
        geom_segment(data = fit_env_df,
                   aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
                   arrow = arrow(length = unit(0.25, "cm")), colour = alpha("black",.5)) +
        geom_text(data = text_fit_env_df, aes(x = Dim1, y = Dim2, label = Variable),
                size = 3)
        
        p <- ggdraw(p) + draw_figure_label(label, size=8, position="bottom.right")
    }
    save_plot(paste0("RDA - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}


# dm is optional and only for things that have been generated outside of R (unifrac)
# env.vars = specify continuous vars to plot additional environment variables 
plot.pcoa <- function(map0, otu0, method="euclidean", plot.title, axis1 = 1, axis2 = 2, dm=NULL, env.vars=NULL, flip.axis=NULL, show.stats=TRUE)
{   
    ret <- prep.dm(map0, otu0, dm, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0
    ppc <- cmdscale(ddm,max(axis1,axis2), eig=TRUE) # return eigenvals for calculating % var explained
    pc <- ppc$points
    
    if(!is.null(flip.axis))  pc[,flip.axis] <- -1*pc[,flip.axis] # sometimes we'll want to flip certain axes so orientations are the same
    
#     ### note, we need to run beta dispersion FIRST to make sure that adonis usage is allowed here (i.e. dispersions must NOT be different in order to use adonis)
#    a <- anova(betadisper(ddm, map0$Sample.Group))
#    print(a)

    # adonis on residence class as factors for ALL samples
    adonis.sample.group <- adonis(dm ~ map0$Sample.Group)
    adonis.sample.group.pval <- adonis.sample.group$aov.tab[1,"Pr(>F)"]     
    adonis.sample.group.R2 <- adonis.sample.group$aov.tab[1,"R2"]
    print(adonis.sample.group)
 
    anosim.result <- anosim(dm, grouping=map0$Sample.Group)
    estimate <- anosim.result$statistic
    pval <- anosim.result$signif
    print(anosim.result)

    # adonis on continuous years.in.us for 1st generation ONLY
   firstgen <- rownames(map0)[map0$Sample.Group %in% c("Hmong1st","Karen1st")]    
    adonis.years <- adonis(dm[firstgen,firstgen] ~ map0[firstgen,"Years.in.US"])
    print(adonis.years)
    adonis.years.pval <- adonis.years$aov.tab[1,"Pr(>F)"]     
    adonis.years.R2 <- adonis.years$aov.tab[1,"R2"]

    # see how PC1 is correlated with years in US for everyone but 2nd generation
#     minus2nd <- rownames(map0)[map0$Sample.Group != "Hmong2nd"]
#     x <- pc[minus2nd,axis1]
#     y <- map0[minus2nd,"Years.in.US"]
     cortest <- cor.test(pc[,axis1],map0[,"Years.in.US"], method="spear")
     print(cortest)
     #cor.label <- paste0("PC1 x Yrs.in.US Rho = ", signif(cortest$estimate, 2), ", P = ", signif(cortest$p.value, 2))
                     
    # let's draw 95% standard error ellipses around Thai and 2nd gen samples
#    mod <- betadisper(ddm, map0$Sample.Group, type="centroid")
#    plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
        
    # calculate variability explained
    percent_var <- signif(eigenvals(ppc)/sum(eigenvals(ppc)), 4)*100
    
    
    d <- data.frame(x = pc[,axis1], y = pc[,axis2], group=map0$Sample.Group)
    # set the levels of Sample.Group so that it's the same every time
    d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))

    group.cols <- alpha(c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"),.8)
    if(length(unique(d$group)) > 5)
        group.cols <- c(alpha("black", .5), group.cols)

    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=group), size=2) + xlab(paste0("PC",axis1, " [",percent_var[axis1],"%]")) + ylab(paste0("PC",axis2, " [",percent_var[axis2],"%]")) + ggtitle(plot.title) +
        scale_color_manual(values=group.cols) + #sets the color palette of the fill
        stat_ellipse(data=d, aes(colour=group), show.legend=F, type="t", level=.6)


    label<-NULL
    if(!is.null(env.vars))
    {
        # Note, this will fit Env.Var ~ PC1 and PC2 -- and report the Multiple-R-squared of this fit (permutation based)
        # This is NOT the same as correlating Env.Var directly with PC1 or PC2 only (above)
        
        # Fit environmental variables
        fit_env <- envfit(ord=pc, env=map0[,env.vars])

        # need to plot these (invisibly) first then grab scaling and vectors out
        # try plotting PCOA with env variables
        oplot <- ordiplot(pc) # plots points only, plot is necessary for scaling vectors below
        # myplot <- plot(fit_env, p.max = 0.05) # plot significant environmental vectors
        vector.multiplier <- ordiArrowMul(scores(fit_env, choices=c(1,2), display="vectors"))
        fit_env_df <- as.data.frame(scores(fit_env, choices=c(1,2), display="vectors")*vector.multiplier)
#        fit_env_df <- cbind(fit_env_df, Variable = rownames(fit_env_df))
        fit_env_df <- cbind(fit_env_df, Variable = env.vars)  

        text_fit_env_df <- fit_env_df
        text_fit_env_df[,1] <- text_fit_env_df[,1]*.7 # let's shift the text out a little
        text_fit_env_df[,2] <- text_fit_env_df[,2]*1.3
        
# manually place text
#         text_fit_env_df[1,1] <- text_fit_env_df[1,1]*.3 # BMI
#         text_fit_env_df[2,1] <- text_fit_env_df[2,1]*.7 # Years X
#         text_fit_env_df[2,2] <- text_fit_env_df[2,2]*1.3 # Years Y
#         text_fit_env_df[3,2] <- text_fit_env_df[3,2]*1.2 # Age
        
        p <- p + coord_fixed() + ## need aspect ratio of 1!
        geom_segment(data = fit_env_df,
                   aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
                   arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
        geom_text(data = text_fit_env_df, aes(x = Dim1, y = Dim2, label = Variable),
                size = 4, fontface = "bold")
        
        # let's print all of the environmental variable correlations 
        #        print(fit_env$vectors)
        # let's not include these stats on the figure itself
        #         stats_df <- cbind(variable=rownames(fit_env_df), R2=fit_env$vectors$r, P=fit_env$vectors$pvals)
        #         label <- paste(apply(stats_df, 1, function(xx) paste0(xx[1], " R2=", signif(as.numeric(xx[2]),2), " P=", signif(as.numeric(xx[3]),2))), collapse="\n")
        #         label <- paste0(label, "\n")
    }
    if(show.stats)
    {
        adonis.label <- paste0("1st-gen R2 = ", signif(adonis.years.R2, 2), ", P = ", signif(adonis.years.pval, 2), 
                                "\nGroup R2 = ", signif(adonis.sample.group.R2, 2), ", P = ", signif(adonis.sample.group.pval, 2))
        p <- ggdraw(p) + draw_figure_label(paste0(adonis.label), size=8, position="bottom.right")
    }
    
    
    save_plot(paste0("pcoa - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

    invisible(p)
}



# convex.hulls are shaded backgrounds around subject sample points
plot.pcoa.long <- function(map0, samples, otu0, method="euclidean", plot.title, dm=NULL, convex.hull=FALSE)
{
    ret <- prep.dm(map0, otu0, dm, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0
    pc <- cmdscale(ddm,2)
                
    group <- rep("NA", nrow(pc))
    names(group) <- rownames(pc)
    group[samples] <- map0[samples,"Subject.ID"]
    d <- data.frame(x = pc[,1], y = pc[,2], group=group, row.names=rownames(pc))

    cols <- alpha(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(map0[samples,"Subject.ID"]))),.8)    

    first.samples <- rownames(map0[which(rownames(map0) %in% samples & map0[,"Sample.Order"] == min(map0[samples,"Sample.Order"])),])
    last.samples <- rownames(map0[which(rownames(map0) %in% samples & map0[,"Sample.Order"] == max(map0[samples,"Sample.Order"])),])
    
    
    background.df <- d[!(rownames(d) %in% samples),] # all other CS samples
    if(length(cols) == 1) # single subject
    {    
        map.000 <- map0[samples,]
        map.000$Sample.Day.Since.First.Sample <- as.numeric(as.Date(map.000$Sample.Date, format="%m/%d/%y") - as.Date("08/01/16", format="%m/%d/%y")) # hard code first sample date
        map.000[,"travel.phase"] <- "Traveling"
        map.000[map.000$Sample.Day.Since.First.Sample <= 2,"travel.phase"] <- "Pre"
        map.000[map.000$Sample.Day.Since.First.Sample >= 28,"travel.phase"] <- "Post"
        map.000$travel.phase <- factor(map.000$travel.phase, levels=c("Pre", "Traveling", "Post")) 
        d2 <- data.frame(d[samples,], group2=map.000$travel.phase)

        p <- ggplot(data=d, aes(x, y)) +
            geom_point(data = background.df, colour=alpha('gray',.2), shape=1, size=2) +
            geom_point(data = d2, aes(colour=group2, shape=19), size=2) +
            scale_shape_identity() +
            xlab("PC1") + ylab("PC2") +
            ggtitle(plot.title) + scale_color_manual(values=alpha(c("red","gray","blue"),.7)) 

    }
    else{
        p <- ggplot(data=d, aes(x, y)) +
            geom_point(data = background.df, colour=alpha('gray',.2), shape=1, size=2) +
            geom_point(data = d[samples,], aes(colour=group, shape=20), size=2) +
            geom_point(data = d[last.samples,], mapping=aes(colour=group, shape=17), size=3) + # ends in triangle
            geom_point(data = d[first.samples,], mapping=aes(colour=group, shape=16), size=3) + # starts in large circle
            scale_shape_identity() +
            xlab("PC1") + ylab("PC2") +
            ggtitle(plot.title) + scale_color_manual(values=cols)

        if(convex.hull)
            p <- p + stat_chull(data=d[samples,], mapping=aes(colour=group), alpha = 0.1, geom = "polygon")
    }

    save_plot(paste0("pcoa - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}
