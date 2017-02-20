require(reshape2) 
require(plyr)
require(ggplot2)
plot.b.p.ratio <- function(map0, otu)
{
    bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
    prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"


    residency <- vector(mode="character", length=nrow(map0))
    residency[map0$Years.in.US > 0 & map0$Years.in.US < 5] <- "0-5"
    residency[map0$Years.in.US >= 5 & map0$Years.in.US < 10] <- "5-10" 
    residency[map0$Years.in.US >= 10 & map0$Years.in.US < 20] <- "10-20"
    residency[map0$Years.in.US >= 20 & map0$Years.in.US < 30] <- "20-30"
    residency[map0$Years.in.US >= 30] <- "> 30"
    
    residency[map0$Years.in.US == 0] <- "US-born"
    residency[is.na(map0$Years.in.US)] <- "Thai-resident"
    
    residency <- factor(residency, levels=c("Thai-resident", "0-5","5-10","10-20","20-30","> 30","US-born")) 
    
    map0[,"residency"] <- residency
    
    otu0 <- melt(otu[,c(bacteroides,prevotella)], id.vars = 0, variable.name = "Taxa",
                value.name = "RelativeAbundance")
    colnames(otu0)[1:2] <- c("SampleID","Taxa")
    merged <- merge(otu0, map0, by.x="SampleID", by.y=0)


    p <- ggplot(merged[order(merged$Taxa),], aes(x = residency, y = RelativeAbundance, fill = Taxa)) +
        geom_bar(stat = "summary", fun.y = "mean") + # summarizes the relative abundances across samples by averaging them
        theme_bw() + # white background
        scale_fill_manual(labels = c("Bacteroides","Prevotella"), values=c("#ae017e", "#fcc5c0")) + # legend labels and colors
        ggtitle("Bacteroides and Prevotella over US residency") + # title
        theme(legend.key = element_blank()) # removes borders from legend
    
    ggsave(plot=p,"b.p.ratio.pdf", useDingbats=FALSE)

}