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
plot.pcoa <- function(map0, otu0, method="euclidean", plot.title, axis1 = 1, axis2 = 2, dm=NULL)
{   
    ret <- prep.dm(map0, otu0, dm, method, axis1, axis2)
    pc <- ret$pc
    dm <- ret$dm
    map0 <- ret$map0
 
#     ### note, we need to run beta dispersion FIRST to make sure that adonis usage is allowed here (i.e. dispersions must NOT be different in order to use adonis)
#    a <- anova(betadisper(ddm, map0$Sample.Group))
#    print(a)

    # adonis on residence class as factors for ALL samples
     adonis.residence <- adonis(dm ~ map0$Sample.Group)
     print(adonis.residence)
 
    anosim.result <- anosim(dm, grouping=map0$Sample.Group)
    estimate <- anosim.result$statistic
    pval <- anosim.result$signif
 
    # adonis on continuous years.in.us for 1st generation ONLY
    firstgen <- rownames(map0)[map0$Sample.Group %in% c("Hmong1st","Karen1st")]    
     adonis.years <- adonis(dm[firstgen,firstgen] ~ map0[firstgen,"Years.in.US"])
     print(adonis.years)
     adonis.years.pval <- adonis.years$aov.tab[1,"Pr(>F)"]     

 
    # see how PC1 is correlated with years in US for everyone but 2nd generation
    minus2nd <- rownames(map0)[map0$Sample.Group != "Hmong2nd"]
    x <- pc[minus2nd,axis1]
    y <- map0[minus2nd,"Years.in.US"]
    y[is.na(y)] <- 0
    cortest <- cor.test(x,y,method="spear")

    
    # let's draw 95% standard error ellipses around Thai and 2nd gen samples
#    mod <- betadisper(ddm, map0$Sample.Group, type="centroid")
#    plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
        
    d <- data.frame(x = pc[,axis1], y = pc[,axis2], group=map0$Sample.Group)
    # set the levels of Sample.Group so that it's the same every time
    d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))

#    group.cols <- c(alpha(wes_palette(n=5, name="Moonrise3"), .8))
#    group.cols <- alpha(brewer.pal(5,"YlGnBu"), .8)
#    group.cols[3] <- "#dfc27d"
    group.cols <- alpha(c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"),.8)
    if(length(unique(d$group)) > 5)
        group.cols <- c(alpha("black", .5), group.cols)
    
    
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=group), size=2) + xlab(paste0("PC",axis1)) + ylab(paste0("PC",axis2)) + ggtitle(plot.title) +
        scale_color_manual(values=group.cols) + #sets the color palette of the fill
        stat_ellipse(data=d, aes(colour=group), show.legend=F, type="t", level=.6)

    cor.label <- substitute(paste("PC1 x Yrs.in.US ", rho, " = ", estimate, ", P = ", pvalue),
                    list(estimate = signif(cortest$estimate, 2), pvalue = signif(cortest$p.value, 2)))

    perm.label <- substitute(paste("ANOSIM R2 = ", estimate ,", P = ",pvalue),
                    list(estimate = signif(estimate, 2), pvalue = signif(pval, 2)))

    p <- ggdraw(p) + draw_figure_label(paste(cor.label, perm.label, sep="\n"), size=8, position="bottom.right")
    
    save_plot(paste0("pcoa - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}



# convex.hulls are shaded backgrounds around subject sample points
plot.pcoa.long <- function(map0, samples, otu0, method="euclidean", plot.title, dm=NULL, convex.hull=FALSE)
{
    ret <- prep.dm(map0, otu0, dm, method)
    dm <- ret$dm
    pc <- ret$pc
    map0 <- ret$map0
                
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
            geom_point(data = d[last.samples,], mapping=aes(colour=group, shape=17), size=3) +
            geom_point(data = d[first.samples,], mapping=aes(colour=group, shape=16), size=3) +
            scale_shape_identity() +
            xlab("PC1") + ylab("PC2") +
            ggtitle(plot.title) + scale_color_manual(values=cols)

        if(convex.hull)
            p <- p + stat_chull(data=d[samples,], mapping=aes(colour=group), alpha = 0.1, geom = "polygon")
    }

    save_plot(paste0("pcoa - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}
