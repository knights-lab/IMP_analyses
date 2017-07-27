
require(vegan)
require(faraway)
require(RColorBrewer)
require(wesanderson) #some new color palettes
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
    pch[map0$BMI.Class=="Normal"] <- 15 #mark lean as squares
    
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

plot.pcoa.phyloseq <- function(map0, taxa0, otu0)
{
    OTU = otu_table(otu0, taxa_are_rows = TRUE)
    #TAX = tax_table(taxa0)
    physeq = phyloseq(OTU)
    sampledata = sample_data(map0)

    GP = physeq
    wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
    GP1 = prune_taxa(wh0, GP)
    GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

    GP.ord <- ordinate(physeq, "NMDS", "bray")
    p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
    print(p1)


}

plot.pcoa <- function(map0, dm, plot.title, otus=NULL)
{
    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    dm <- dm[valid_samples, valid_samples]
    
    ddm <- as.dist(dm)
    # plot clusters
    # get PCoA scores first
    pc <- cmdscale(ddm,2)
    
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
    x <- pc[minus2nd,1]
    y <- map0[minus2nd,"Years.in.US"]
    y[is.na(y)] <- 0
    cortest <- cor.test(x,y,method="spear")
#    plot(x,y)

    
    # let's draw 95% standard error ellipses around Thai and 2nd gen samples
#    mod <- betadisper(ddm, map0$Sample.Group, type="centroid")
#    plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
        
    d <- data.frame(x = pc[,1], y = pc[,2], group=map0$Sample.Group)
    # set the levels of Sample.Group so that it's the same every time
    d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))

#    group.cols <- c(alpha(wes_palette(n=5, name="Moonrise3"), .8))
#    group.cols <- alpha(brewer.pal(5,"YlGnBu"), .8)
#    group.cols[3] <- "#dfc27d"
    group.cols <- alpha(c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"),.8)
    if(length(unique(d$group)) > 5)
        group.cols <- c(alpha("black", .5), group.cols)
    
    
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=group), size=2) + xlab("PC1") + ylab("PC2") + ggtitle(plot.title) +
        scale_color_manual(values=group.cols) + #sets the color palette of the fill
        stat_ellipse(data=d, aes(colour=group), show.legend=F, type="t", level=.6)

    cor.label <- substitute(paste("MB Corr ", rho, " = ", estimate, ", P = ", pvalue),
                    list(estimate = signif(cortest$estimate, 2), pvalue = signif(cortest$p.value, 2)))

    perm.label <- substitute(paste("ANOSIM R2 = ", estimate ,", P = ",pvalue),
                    list(estimate = signif(estimate, 2), pvalue = signif(pval, 2)))

    p <- ggdraw(p) + draw_label(cor.label, .8, .05, size=8) + draw_label(perm.label, .8, .1,size=8)
    
    save_plot(paste0("pcoa - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}

plot.pcoa.long <- function(map0, dm, plot.title, otus=NULL)
{
    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    dm <- dm[valid_samples, valid_samples]
    
    ddm <- as.dist(dm)
    # plot clusters
    # get PCoA scores first
    pc <- cmdscale(ddm,2)
            
    ix.start <- which((map0$Sample.Order == 0) & map0$Subject.ID != "IMP.000")
    subs <- map0[(map0$Sample.Order == 0) & map0$Subject.ID != "IMP.000", "Subject.ID"]
    ix.start <- c(ix.start, which((map0$Sample.Order == 1) & map0$Subject.ID != "IMP.000" & !(map0$Subject.ID %in% subs)))
    ix.end <- which((map0$Sample.Order == 6) & map0$Subject.ID != "IMP.000")
    startend <- rep("NA", nrow(pc))
    startend[ix.start] <- "Begin"
    startend[ix.end] <- "End"    

    d <- data.frame(x = pc[,1], y = pc[,2], group=startend, subject=map0$Subject.ID)

    d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))

    group.cols <- c("red", "blue",alpha("gray",.5))
    
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=group), size=2) + xlab("PC1") + ylab("PC2") + ggtitle(plot.title) +
        scale_color_manual(values=group.cols) 
        
    p <- p + geom_line(data=d[d$group %in% c("Begin","End"),], aes(x,y,group=subject))
    
    save_plot(paste0("pcoa - ", plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}

