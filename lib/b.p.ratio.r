require(reshape2) 
require(plyr)
require(ggplot2)
require(cowplot)

# plot.b.p.ratio.L <- function(map0, otu, bug1, bug2, outputfn)
# {        
#     otu <- otu[rownames(map0),]
# 
#     # add .00001 to any that are 0 abundance
#     otu[otu[,bug1] == 0, bug1] <- .00001
#     otu[otu[,bug2] == 0, bug2] <- .00001
#     
#     bp <- otu[,bug1]/otu[,bug2]
#     names(bp) <- rownames(otu)
# 
#     lookup <- c(19,17) # point type
#     names(lookup) <- sort(unique(map0$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
#     pch <- lookup[as.character(map0$Ethnicity)] 
# 
#     d <- data.frame(x = map0$Sample.Order, y = log10(bp), id = map0$Subject.ID)
#     p <- ggplot(data=d, aes(x, y, color = id)) + geom_point() + geom_line() +
#         theme_bw() + # white background
#         ggtitle("Bacteroides-Prevotella Ratio - Longitudinal") + # title
#         ylab("log10(Bacteroides-Prevotella Ratio)") + xlab("Month in the US") +
#         scale_x_continuous(breaks=1:6) +
#         theme(legend.key = element_blank()) # removes borders from legend
#             
#     ggsave(plot=p,outputfn, useDingbats=FALSE)
# }


plot.b.p.ratio <- function(map0, otu, bug1, bug2, outputfn, longitudinal=F)
{        
    otu <- otu[rownames(map0),]

    # add .00001 to any that are 0 abundance
    otu[otu[,bug1] == 0, bug1] <- .00001
    otu[otu[,bug2] == 0, bug2] <- .00001
    
    bp <- otu[,bug1]/otu[,bug2]
    names(bp) <- rownames(otu)

    pch <- get.pch(map0)
        
    p <- NULL
    if(longitudinal)
    {
        d <- data.frame(x = map0$Sample.Order, y = log10(bp), id = map0$Subject.ID)
        p <- ggplot(data=d, aes(x, y, color = id)) + geom_point() + geom_line() + xlab("Month in the US") 

    }
    else
    {
        d <- data.frame(x = map0$Years.in.US, y = log10(bp))
        p <- ggplot(data=d, aes(x, y)) + geom_point(color = alpha("black", .3), shape=pch) + geom_smooth(color="#ae017e", size=.5) + xlab("Years in the US")
    }
    p <- p + 
        theme_bw() + # white background
        ggtitle("Bacteroides-Prevotella Ratio over US residency") + # title
        ylab("log10(Bacteroides-Prevotella Ratio)") +
        theme(legend.key = element_blank()) # removes borders from legend

    save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
    
    invisible(bp)
}


plot.b.p.ratio.all <- function(map0, otu, bug1, bug2, outputfn, g1, g2, g3)
{
    map0 <- map0[c(g1,g2,g3),]
    otu <- otu[rownames(map0),]

    # add .00001 to any that are 0 abundance
    otu[otu[,bug1] == 0, bug1] <- .00001
    otu[otu[,bug2] == 0, bug2] <- .00001
    
    bp <- otu[,bug1]/otu[,bug2]
    names(bp) <- rownames(otu)

    pch <- get.pch(map0)

    d <- data.frame(x = map0[g1, "Years.in.US"], y = log10(bp[g1]))
    d2 <-data.frame(x=rep(-2, length(g2)), y=log10(bp[g2]))
    d3 <-data.frame(x=rep(43, length(g3)), y=log10(bp[g3]))

    
    main <- ggplot() + geom_point(data=d, aes(x=x,y=y), color=alpha("black",.3)) + geom_smooth(data=d, aes(x=x,y=y), color="#ae017e", size=.5) +
        geom_point(data=d2, aes(x=x,y=y), color=alpha("#ae017e",.3), pch=6) +
        geom_point(data=d3, aes(x=x,y=y), color=alpha("#ae017e",.3), pch=2) +
        xlim(-2,43) + 
        ggtitle("Bacteroides-Prevotella Ratio over US residency") + # title
        ylab("log10(Bacteroides-Prevotella Ratio)") + xlab("Years in the US")
    
    final <- ggdraw(main) + draw_label("Thai", angle=45, colour="#ae017e", size = 8, x = 0.165, y = 0.105) + 
            draw_label("2ndGen", angle=45, colour="#ae017e", size = 8, x = 0.94, y = 0.105)

    
    save_plot(outputfn, final, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}



# plots taxa bar plots of B and P by blocks of time in US
plot.b.p.barplot <- function(map0, otu, bins=seq(0,10,2), fn)
{
    bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
    prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"

    map0$Years.in.US[is.na(map0$Years.in.US)] <- -2
    map0$Years.in.US[map0$Years.in.US==0] <- 45

    bin.labels <- paste0("<",bins[-1])
    map0$Decade <- cut(map0$Years.in.US, bins, labels=bin.labels)
    map0$Decade <- as.character(map0$Decade)
    map0[map0$Years.in.US==-2,"Decade"] <- "Thai"

    if(sum(map0$Years.in.US==45)==0) # if no US born
        map0$Decade <- factor(map0$Decade, levels=c("Thai", bin.labels)) 
    else {
        map0[map0$Years.in.US==45,"Decade"] <- "USborn"
        map0$Decade <- factor(map0$Decade, levels=c("Thai", bin.labels, "USborn")) 
    }
    
    otu0 <- melt(otu[,c(bacteroides,prevotella)], id.vars = 0, variable.name = "Taxa",
                value.name = "RelativeAbundance")
    colnames(otu0)[1:2] <- c("SampleID","Taxa")
    merged <- merge(otu0, map0, by.x="SampleID", by.y=0)

    merged$Taxa <- factor(merged$Taxa, levels=c(prevotella,bacteroides)) 


    p <- ggplot(merged[order(merged$Taxa),], aes(x = Decade, y = RelativeAbundance, fill = Taxa)) +
        geom_bar(stat = "summary", fun.y = "mean") + # summarizes the relative abundances across samples by averaging them
        scale_fill_manual(labels = c("Prevotella","Bacteroides"), values=c("#fcc5c0","#ae017e")) + # legend labels and colors
        ggtitle("Bacteroides and Prevotella over US residency") 
    
    save_plot(fn, p, useDingbats=FALSE, base_aspect_ratio = 1.3)


}