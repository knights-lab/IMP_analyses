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