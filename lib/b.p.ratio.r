require(reshape2) 
require(plyr)
require(ggplot2)
require(cowplot)

plot.b.p.ratio.x.bmi <- function(map0, otu, bug1, bug2, outputfn=NULL)
{        
    otu <- otu[rownames(map0),]

    # add .00001 to any that are 0 abundance
    otu[otu[,bug1] == 0, bug1] <- .00001
    otu[otu[,bug2] == 0, bug2] <- .00001
    
    bp <- otu[,bug1]/otu[,bug2]
    names(bp) <- rownames(otu)
        
    p <- NULL
    d <- data.frame(x = map0$Years.in.US, y = log10(bp), group=map0$BMI.Class)
    p <- ggplot(data=d, aes(x, y, color=group)) + geom_point() + 
        geom_smooth(method = "nls", formula = y ~ a * x + b, se = F, method.args = list(start = list(a = 0.1, b = 0.1))) + 
        scale_color_manual(name = "BMI Class", values = alpha(c("#31625b","#fcd167","#c45565"),.5)) +
        xlab("Years in the US") +
        ggtitle("Bacteroides-Prevotella Ratio over US residency") + 
        ylab("log10(Bacteroides-Prevotella Ratio)") + theme(legend.position='none') # remove legend

    if(is.null(outputfn))
        return(p)
    else
        save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}

# num.clip.months = number of months to clip and average at the start and end
plot.b.p.ratio.L <- function(map0, otu0, bug1, bug2, outputfn, num.clip.months=NULL, show.stats)
{        
    otu0 <- otu0[rownames(map0),]
    
    bp <- get.log.taxa.ratio(otu0, bug1, bug2)

    p <- plot.response.L(map0=map0, y=bp, ggtitle="Bacteroides/Prevotella", ylab="log10(B/P)", num.clip.months=num.clip.months, show.stats=show.stats)
    return(p)
}

get.log.taxa.ratio <- function(otu, bug1, bug2)
{    
    if(sum(otu[,bug1] == 0) > 0){
        print(paste0("warning zero detected in ", bug1))
    }
    if(sum(otu[,bug2] == 0) > 0){
        print(paste0("warning zero detected in ", bug2))
    }

    bp <- log10(otu[,bug1]/otu[,bug2])
    return(bp)
    
    # reverse the trend - we shouldn't have to do this anymore!
    # flip the ratios to avoid dividing by 0 (prevotella only)
#     bp <- otu[,bug2]/otu[,bug1]
#     names(bp) <- rownames(otu)
#     # flip sign to get it back to bp
#     bp <- -1*log10(bp)
}

plot.bp.main <- function(ggdata, hide.y=TRUE)
{
    main <- ggplot(data=ggdata, aes(x=Years.in.US,y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5, stroke=.75, fill=NA) +
        geom_smooth(data=ggdata, aes(x=Years.in.US,y=y), color="black", size=.5) +
        scale_color_manual(name="Groups", values=cols) + 
        scale_alpha_manual(name="Groups", values=alphas) +
        scale_shape_manual(name="Groups", values=shapes) +
        ylim(ylim)  +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ylab("log10(B/P)") + xlab("Years in the US") + theme(legend.position='none') + 
        
    if(hide.y)# remove all y axis
        main <- main + theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size=10))

    return(main)
}

plot.b.p.ratio.all <- function(map0, otu0, bug1, bug2, outputfn, do.insets=FALSE)
{
    # don't include any samples that have completely zero of any because it screws up the log10 transform
    otu0 <- otu0[otu0[,bug1]!=0,]
    otu0 <- otu0[otu0[,bug2]!=0,]
    
    validsamples <- intersect(rownames(map0),rownames(otu0))
    otu0 <- otu0[validsamples,]
    map0 <- map0[validsamples,]
    
    bp <- get.log.taxa.ratio(otu0, bug1, bug2)

    d <- data.frame(map0, y = bp)
    ylim <- range(d$y)
    d.1st <- d[d$Years.in.US > 0 & d$Years.in.US < 50,]
    d.2nd.control <- d[d$Years.in.US %in% c(50,60),]
    d.2nd.control$Sample.Group <- factor(d.2nd.control$Sample.Group)
    d.thai <-  d[d$Years.in.US == 0,]
    d.thai$Sample.Group <- factor(d.thai$Sample.Group)
    
    groupnames <- as.character(unique(d$Sample.Group))

    cols <- get.group.colors(groups=groupnames) 
    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 

    main <- ggplot(data=d.1st, aes(x=Years.in.US,y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5, stroke=.75, fill=NA) +
        geom_smooth(data=d.1st, aes(x=Years.in.US,y=y), color="black", size=.5, method="lm") +
        scale_color_manual(name="Groups", values=cols) + 
        scale_alpha_manual(name="Groups", values=alphas) +
        scale_shape_manual(name="Groups", values=shapes) +
        ylim(ylim)  +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ylab("") + xlab("Years in the US") + theme(legend.position='none') + 
        # remove all y axis
        theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size=10))
    print(cor.test(d.1st$Years.in.US, d.1st$y, method="spear", exact=F))

    if(do.insets)
    {
        dh <- d.1st[d.1st$Ethnicity=="Hmong",]
        print(cor.test(dh$Years.in.US, dh$y, method="spear", exact=F))
        main.h <- ggplot(data=dh, aes(x=Years.in.US,y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5, stroke=.75, fill=NA) +
            geom_smooth(aes(x=Years.in.US,y=y), color="black", size=.5, method="lm") + ggtitle("Hmong") +
            scale_color_manual(name="Groups", values=cols) + 
            scale_alpha_manual(name="Groups", values=alphas) +
            scale_shape_manual(name="Groups", values=shapes) +
            ylab("log10(B/P)") + xlab("Years in the US") + theme(legend.position='none')

        dk <- d.1st[d.1st$Ethnicity=="Karen",]
        print(cor.test(dk$Years.in.US, dk$y, method="spear", exact=F))
        main.k <- ggplot(data=dk, aes(x=Years.in.US,y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5, stroke=.75, fill=NA) +
            geom_smooth(aes(x=Years.in.US,y=y), color="black", size=.5, method="lm") + ggtitle("Karen") +
            scale_color_manual(name="Groups", values=cols) + 
            scale_alpha_manual(name="Groups", values=alphas) +
            scale_shape_manual(name="Groups", values=shapes) +
            ylab("log10(B/P)") + xlab("Years in the US") + theme(legend.position='none')
    
        # generate these as insets
        save_plot("bp.hmong1st.pdf", main.h, base_aspect_ratio=1)
        save_plot("bp.karen1st.pdf", main.k, base_aspect_ratio=1)
    }
    left <- ggplot(d.thai, aes(Sample.Group, y, color=Sample.Group, alpha=Sample.Group)) + geom_quasirandom(dodge.width=.75) + geom_boxplot(alpha=0, colour="black") +  
            scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
            scale_alpha_manual(name="Groups", values=alphas) +
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("log10(B/P)") + xlab("") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),axis.text = element_text(size=10), axis.title.y = element_text(size=12)) + scale_x_discrete(labels=SAMPLE.GROUP.NAMES.SHORT)
    
    right <- ggplot(d.2nd.control, aes(Sample.Group, y, color=Sample.Group)) + 
            geom_quasirandom(dodge.width=.75, aes(color=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5) +
            geom_boxplot(alpha=0, colour="black") +  
            scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
            scale_shape_manual(name="Groups", values=shapes) +
            scale_alpha_manual(name="Groups", values=alphas) +
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("") + xlab("")  +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + scale_x_discrete(labels=SAMPLE.GROUP.NAMES.SHORT) +
            # remove all y axis
            theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.text.x = element_text(size=10))

    plots <- plot_grid(left, main, right, nrow=1, rel_widths=c(1.5,2,1), align="h") 
    save_plot(outputfn, plots, useDingbats=FALSE, base_aspect_ratio = 1)
}

plot.b.p.ratio.all.nocolor <- function(map0, otu0, bug1, bug2, outputfn)
{
    # don't include any samples that have completely zero of any because it screws up the log10 transform
    otu0 <- otu0[otu0[,bug1]!=0,]
    otu0 <- otu0[otu0[,bug2]!=0,]
    
    validsamples <- intersect(rownames(map0),rownames(otu0))
    otu0 <- otu0[validsamples,]
    map0 <- map0[validsamples,]
    
    bp <- get.log.taxa.ratio(otu0, bug1, bug2)

    d <- data.frame(map0, y = bp)
    ylim <- range(d$y)
    d.1st <- d[d$Years.in.US > 0 & d$Years.in.US < 50,]
    d.2nd.control <- d[d$Years.in.US %in% c(50,60),]
    d.2nd.control$Sample.Group <- factor(d.2nd.control$Sample.Group)
    d.thai <-  d[d$Years.in.US == 0,]
    d.thai$Sample.Group <- factor(rep(d.thai$Sample.Group[1],nrow(d.thai)))
    
    groupnames <- as.character(unique(d$Sample.Group))

    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 

    alphaval=.3

    main <- ggplot(data=d.1st, aes(x=Years.in.US,y)) + geom_point(alpha=alphaval, fill="black") +
        geom_smooth(data=d.1st, aes(x=Years.in.US,y=y), color="black", size=.5, method="lm") +
        ylim(ylim)  +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ylab("") + xlab("Years in the US") + theme(legend.position='none') + 
        # remove all y axis
        theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size=10))
    print(cor.test(d.1st$Years.in.US, d.1st$y, method="spear", exact=F))

    left <- ggplot(d.thai, aes(x=Sample.Group, y)) + geom_quasirandom(dodge.width=.75, alpha=alphaval, color="black") + geom_boxplot(alpha=0, colour="black") +  
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("log10(B/P)") + xlab("") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(size=10), axis.title.y = element_text(size=12)) + scale_x_discrete(labels="Thai")
    
    right <- ggplot(d.2nd.control, aes(Sample.Group, y)) + 
            geom_quasirandom(dodge.width=.75, color="black", alpha=alphaval) +
            geom_boxplot(alpha=0, colour="black") +  
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("") + xlab("")  +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + scale_x_discrete(labels=SAMPLE.GROUP.NAMES.SHORT) +
            # remove all y axis
            theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.text.x = element_text(size=10))

    plots <- plot_grid(left, main, right, nrow=1, rel_widths=c(1.5,2,1), align="h") 
    save_plot(outputfn, plots, useDingbats=FALSE, base_aspect_ratio = 1)
}

plot.b.p.ratio.all.onethai <- function(map0, otu0, bug1, bug2, outputfn)
{
    # don't include any samples that have completely zero of any because it screws up the log10 transform
    otu0 <- otu0[otu0[,bug1]!=0,]
    otu0 <- otu0[otu0[,bug2]!=0,]
    
    validsamples <- intersect(rownames(map0),rownames(otu0))
    otu0 <- otu0[validsamples,]
    map0 <- map0[validsamples,]
    
    bp <- get.log.taxa.ratio(otu0, bug1, bug2)

    d <- data.frame(map0, y = bp)
    ylim <- range(d$y)
    d.1st <- d[d$Years.in.US > 0 & d$Years.in.US < 50,]
    d.2nd.control <- d[d$Years.in.US %in% c(50,60),]
    d.2nd.control$Sample.Group <- factor(d.2nd.control$Sample.Group)
    d.thai <-  d[d$Years.in.US == 0,]
    d.thai$One.Sample.Group <- factor(rep(d.thai$Sample.Group[1],nrow(d.thai)))
    
    groupnames <- as.character(unique(d$Sample.Group))

    cols <- get.group.colors(groups=groupnames) 
    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 

    main <- ggplot(data=d.1st, aes(x=Years.in.US,y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5, stroke=.75, fill=NA) +
        geom_smooth(data=d.1st, aes(x=Years.in.US,y=y), color="black", size=.5, method="lm") +
        scale_color_manual(name="Groups", values=cols) + 
        scale_alpha_manual(name="Groups", values=alphas) +
        scale_shape_manual(name="Groups", values=shapes) +
        ylim(ylim)  +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ylab("") + xlab("Years in the US") + theme(legend.position='none') + 
        # remove all y axis
        theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size=10))

    print(cor.test(d.1st$Years.in.US, d.1st$y, method="spear", exact=F))

    left <- ggplot(d.thai, aes(x=One.Sample.Group, y, color=Sample.Group, alpha=Sample.Group, group=One.Sample.Group)) + geom_quasirandom(dodge.width=.75) + geom_boxplot(alpha=0, colour="black") +  
            scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
            scale_alpha_manual(name="Groups", values=alphas) +
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("log10(B/P)") + xlab("") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),axis.text = element_text(size=10), axis.title.y = element_text(size=12)) + scale_x_discrete(labels="HT-KT")

    
    right <- ggplot(d.2nd.control, aes(Sample.Group, y, color=Sample.Group)) + 
            geom_quasirandom(dodge.width=.75, aes(color=Sample.Group, shape=Sample.Group, alpha=Sample.Group), size=1.5) +
            geom_boxplot(alpha=0, colour="black") +  
            scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
            scale_shape_manual(name="Groups", values=shapes) +
            scale_alpha_manual(name="Groups", values=alphas) +
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("") + xlab("")  +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + scale_x_discrete(labels=SAMPLE.GROUP.NAMES.SHORT) +
            # remove all y axis
            theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.text.x = element_text(size=10))

    plots <- plot_grid(left, main, right, nrow=1, rel_widths=c(1.5,2,1), align="h") 
    save_plot(outputfn, plots, useDingbats=FALSE, base_aspect_ratio = 1)
}


plot.b.p.ratio.all.boxplots <- function(map0, otu0, bug1, bug2, outputfn)
{
    # don't include any samples that have completely zero of any because it screws up the log10 transform
    otu0 <- otu0[otu0[,bug1]!=0,]
    otu0 <- otu0[otu0[,bug2]!=0,]
    
    validsamples <- intersect(rownames(map0),rownames(otu0))
    otu0 <- otu0[validsamples,]
    map0 <- map0[validsamples,]
    
    bp <- get.log.taxa.ratio(otu0, bug1, bug2)

    d <- data.frame(map0, y = bp)
    
    groupnames <- as.character(unique(d$Sample.Group))

    cols <- get.group.colors(groups=groupnames) 
    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 

#     p <- ggplot(d, aes(Sample.Group, y, color=Sample.Group, shape=Sample.Group, size=Sample.Group, alpha=Sample.Group, group=Sample.Group)) + geom_quasirandom(dodge.width=.75) + 
#             geom_boxplot(alpha=0, colour="black", size=.5) +  
#             scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
#             scale_alpha_manual(name="Groups", values=alphas) +
#             scale_shape_manual(name="Groups", values=shapes) +
#             scale_size_manual(name="Groups", values=sizes) +
#             guides(fill=FALSE) + theme(legend.position='none') +
#             ylab("log10(B/P)") + xlab("") +
#             theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), axis.text = element_text(size=10), axis.title.y = element_text(size=12)) + scale_x_discrete(labels=SAMPLE.GROUP.NAMES.SHORT)
    
    p <- map.boxplot(y=d$y, Group=d$Sample.Group, main="", facet.var=NULL, alpha=.1, add.pval=TRUE, plot.legend.only=FALSE, ylab="log10(B/P)", strip.text.size=5, y.size=10, 
                            x.size=9, show.x=TRUE, group.vars.df=d[,c("Resident.Continent","Birth.Continent","Ethnicity")])

    
    save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1)
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