# plots taxa summaries and orders by 2 largest relative abundances within sample groups
# optional sample.order = predefined order of sample IDs
# x.labels = optional relabels of x axis 
plot.taxa.summary <- function(map0, otu, fn, sample.order=NULL, x.labels=NULL)
{
    valid.samples <- intersect(rownames(map0), rownames(otu))
    map0 <- map0[valid.samples,]
    otu <- otu[valid.samples,]
    
    map0$Sample.Group <- factor(map0$Sample.Group) # drop levels that arent present
    
    ordered.taxa <- names(sort(colSums(otu), decreasing=T))
    # reorder the colors here so that we get distinct colors next to each other 
    # but maintain the same color-taxa assignment
    cols <- gg_color_hue(ncol(otu))
    names(cols) <- colnames(otu)
    reordered.cols <- cols[ordered.taxa]
    
    # order entire OTU by samples by largest rel abundance of taxa    
#    ordered.otu <- otu[ order(-otu[,ordered.taxa[1]], -otu[,ordered.taxa[2]]), ]

    # order samples within each sample group
    ordered.rownames <- NULL
    for(i in 1:length(levels(map0$Sample.Group)))
    {
        group <- levels(map0$Sample.Group)[i]
        group.samples <- rownames(map0)[as.character(map0$Sample.Group) %in% group]
        this.otu <- otu[group.samples,]
        o.this.otu <- this.otu[ order(-this.otu[,ordered.taxa[1]], -this.otu[,ordered.taxa[2]]), ]
        ordered.rownames <- c(ordered.rownames, rownames(o.this.otu))

    }
    ordered.otu <- otu[ordered.rownames,]
    
    # if a sample order is given, override any reordering
    if(length(sample.order) > 0) ordered.otu <- otu[sample.order,]
    
    otu0 <- melt(ordered.otu, id.vars = 0, variable.name = "Taxa",
                value.name = "RelativeAbundance")
    colnames(otu0)[1:2] <- c("SampleID","Taxa")
    merged <- merge(otu0, map0, by.x="SampleID", by.y=0)

    # always order the taxa by their names in alphabetical order
    merged$Taxa <- factor(merged$Taxa, levels=rev(ordered.taxa))

    p <- ggplot(merged[order(merged$Taxa),], aes(x = SampleID, y = RelativeAbundance, fill = Taxa)) +
        geom_bar(stat = "identity", position="fill") +
        scale_x_discrete(labels = x.labels) + theme(legend.text=element_text(size=8)) +
       scale_fill_manual(breaks=ordered.taxa[1:20], values=reordered.cols) + # include only the first 10 taxa in legend
        ggtitle("Taxa Summary Plot")
        
    l <- get_legend(p)    
    
    pg <- plot_grid(p + theme(legend.position="none"))
    ggsave(plot=pg, fn, useDingbats=FALSE)
    ggsave(plot=l, paste0("legend_", fn), useDingbats=FALSE)

    invisible(ordered.rownames)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# cols = named vector of colors
plot.legend <- function(cols, outputfn)
{
  max.y <- length(cols)
  pdf(outputfn, width=10, height=7*(max.y/20))  
  plot(1:max.y, axes = 0, xlab="", ylab="", type="n")
  legend(1, max.y, col="white", legend=names(cols), pch=22, bty="n", pt.cex=3, pt.bg=cols)
  dev.off()
}


# multiplot of longitudinal stream plots
# pass in mapping file that already contains the longitudinal samples of interest
# x.var = timescale variable to plot by
plot.taxa.summary.L <- function(taxa0, map0, outputfn, max.taxa=15, x.var="Sample.Day.Since.Arrival", grid.ncol=1)
{
      valid.samples <- intersect(rownames(map0), rownames(taxa0))
      map0 <- map0[valid.samples,]
      taxa0 <- taxa0[valid.samples,]
  
    # manually normalize here
    #taxa0 <- sweep(taxa0, 1, rowSums(taxa0), '/')
    #prevalences <- apply(taxa0, 2, function(bug.col) mean(bug.col > 0))
    #taxa0 <- taxa0[, prevalences >= .10]

    subjects <- sort(unique(map0$Subject.ID)) # maintain order
   
    p <- NULL
    all.top <- NULL # list of all top taxa to show (for legend purposes)
    for(i in 1:length(subjects))
    {
        this.map0 <- map0[map0$Subject.ID == subjects[i],]
        this.samples <- rownames(this.map0[order(this.map0[,x.var]),]) # order samples by time
        this.top <- names(sort(colMeans(taxa0[this.samples,]), decreasing=T))[1:max.taxa] # order taxa by abundance
#print(subjects[i])
#print(taxa0[this.samples, this.top])

        dm <- melt(taxa0[this.samples, rev(this.top)])
        colnames(dm) <- c("sample.names", "taxa", "rel.abundance")
        d <- data.frame(dm, x=map0[as.character(dm[,1]), x.var])

        p[[i]] <- ggplot(d, aes(x=x, y=rel.abundance, fill=taxa)) +
                geom_area(colour="white", size=.05, alpha=.8) + xlab("") + ylab("") +
                theme(axis.text = element_text(size=8),  legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=7)) + 
                #xlim(min(max(map0[,x.var]),max(map0[,x.var])) +
                ggtitle(subjects[i]) +
                theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) # removes extra spacing between plots
    
        all.top[[i]] <- this.top
    }

    # add global taxa colors and remove legend 
    cols <- colorRampPalette(brewer.pal(11,"Set3"))(length(unique(unlist(all.top))))
    names(cols) <- unique(unlist(all.top))
    cols2 <- cols
    names(cols2) <- shorten.taxonomy(names(cols2))
    
    # manually plot and save legend separately
    plot.legend(cols2, paste0("legend_",outputfn))
    
    # add global colors and remove legend from plots
    for(i in 1:length(p))
        p[[i]] <- p[[i]] + scale_fill_manual(values=cols) + theme(legend.position='none')

    nrow <- floor(length(p) / grid.ncol)
    if(length(p) %% grid.ncol > 0)
      nrow <- nrow + 1
    
    if(length(p)>1)
        multiplot <- plot_grid(plotlist=p, ncol=grid.ncol, nrow=nrow)
    else # IMP.000 only
    {
        p[[1]] <- p[[1]] + geom_vline(xintercept=4, color="white", linetype="dashed") + geom_vline(xintercept=28, color="white", linetype="dashed")
        multiplot <- plot_grid(plotlist=p, ncol=1, nrow=1)
    }
    multiplot <- add_sub(multiplot, gsub("\\.", " ", x.var), vpadding=grid::unit(0,"lines"),y=6, x=0.5, vjust=4.5)
    save_plot(outputfn, multiplot, ncol = grid.ncol, nrow = nrow, base_aspect_ratio = 1.3)
}

