# pick samples for deep shotgun

# load all source files, data, and prep any variables for reuse
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/bin/load.r")

plot.b.p.barplot.2 <- function(map0, otu, fn)
{
    bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
    prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"
    
    otu0 <- melt(otu[,c(bacteroides,prevotella)], id.vars = 0, variable.name = "Taxa",
                value.name = "RelativeAbundance")
    colnames(otu0)[1:2] <- c("SampleID","Taxa")
    merged <- merge(otu0, map0, by.x="SampleID", by.y=0)

    merged$Taxa <- factor(merged$Taxa, levels=c(prevotella,bacteroides)) 

    p <- ggplot(merged, aes(x = Override.Order, y = RelativeAbundance, fill = Taxa)) +
        geom_bar(stat = "summary", fun.y = "mean") + # summarizes the relative abundances across samples by averaging them
        scale_fill_manual(labels = c("Prevotella","Bacteroides"), values=c("#fcc5c0","#ae017e")) + # legend labels and colors
        ggtitle("Bacteroides-Prevotella") +  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
        theme(axis.text=element_text(size=7), axis.text.x = element_text(angle=90, hjust=0))
    
    save_plot(fn, p, useDingbats=FALSE, base_aspect_ratio = 1.3)

}

# exclude these, did NOT pass QC
# CS.171
# CS.082
# CS.083
# CS.010

# only CS.171 is an issue for our purposes 
#map[c("CS.171","CS.082","CS.083","CS.010"),]

dm <- wuf_dm[cs,cs] # best to use unifrac here
ddm <- as.dist(dm)
pc <- cmdscale(ddm,2)

map0 <- map[cs,]

bd <- betadisper(ddm, map0$Sample.Group)
centroid.dist <- bd$distances # distances of all samples to their group centroids!

# KCK full stool samples
samples <- paste0("TFSCS0", 20:29)

sort(centroid.dist[samples])
first5 <- names(sort(centroid.dist[samples])[1:5])

bp <- plot.b.p.ratio(map0[cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="temp.b.p.ratio.pdf")

ordered.samples <- names(sort(bp))

fmt.map <- data.frame(map0[ordered.samples,], Override.Order=ordered.samples)
fmt.map$Override.Order <- factor(fmt.map$Override.Order, levels=fmt.map$Override.Order)
plot.b.p.barplot.2(fmt.map, taxa[ordered.samples,], fn="fmt.barplot.bp.pdf")

d <- data.frame(x = pc[,1], y = pc[,2], group=map0$Sample.Group, distance=substring(as.character(centroid.dist[rownames(map0)]),1,4))
# set the levels of Sample.Group so that it's the same every time
d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))
group.cols <- alpha(c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"),.8)
p <- ggplot(data=d, aes(x, y)) + geom_point(colour=alpha("gray",.5), size=2) +
    scale_color_manual(values=group.cols) + #sets the color palette of the fill
    stat_ellipse(data=d, aes(colour=group), show.legend=T, type="t", level=.6)

p <- p + geom_text(data = d[samples, ], aes(label=samples) ,hjust=0, vjust=0, size=3)
p <- p + geom_point(data=d[hmongthai,], colour="red")
p <- p + geom_point(data=d[first5,], colour="blue")
hmongus <- c("CS.224", "CS.283", "CS.333", "CS.380")
p <- p + geom_point(data=d[hmongus,], colour="yellow")
p <- p + geom_text(data = d[hmongus, ], aes(label=hmongus) ,hjust=0, vjust=0, size=3)

hmong2 <- rownames(map0[map0$Sample.Group=="Hmong2nd",])
d2 <- data.frame(Age=map0$Age, BMI=map0$BMI, group=map0$Sample.Group, distance=substring(as.character(centroid.dist[rownames(map0)]),1,4), row.names=rownames(map0))
hmong2.d2 <- d2[hmong2,]
