require(reshape2) 
require(plyr)
require(ggplot2)
require(cowplot)

plot.b.p.ratio.x.bmi <- function(map0, otu, bug1, bug2, outputfn)
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

    save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}

# num.clip.months = number of months to clip and average at the start and end
plot.b.p.ratio.L <- function(map0, otu0, bug1, bug2, outputfn, num.clip.months=NULL)
{        
    otu0 <- otu0[rownames(map0),]
    
    bp <- get.log.taxa.ratio(otu0, bug1, bug2)
        
    p <- NULL
    d <- data.frame(month = map0$Sample.Order, bp = bp, subject = map0$Subject.ID)
    if(!is.null(num.clip.months)) # take mean of arbitrary num of start and last months
    {
        d$subject <- as.character(d$subject)

        # find first x months, and average their BP
        mins <- NULL
        for(m in 0:(num.clip.months-1))
            mins <- rbind(mins, aggregate(d$month, list(d$subject), FUN=function(xx) min(xx)+m))
        colnames(mins) <- c("subject","month")            
        mins.d <- merge(mins, d, by=c("subject","month"))
        mean.bp1 <- cbind(aggregate(mins.d$bp, list(mins.d$subject), FUN=mean), month="start")

        # find last x months, and average their BP
        maxs <- NULL
        for(m in 0:(num.clip.months-1))
            maxs <- rbind(maxs, aggregate(d$month, list(d$subject), FUN=function(xx) max(xx)-m))
        colnames(maxs) <- c("subject","month")
        maxs.d <- merge(maxs, d, by=c("subject","month"))
        mean.bp2 <- cbind(aggregate(maxs.d$bp, list(maxs.d$subject), FUN=mean), month="end")

        mean.d <- rbind(mean.bp1, mean.bp2)
        colnames(mean.d) <- c("subject", "bp", "month")
        d <- rbind(d, mean.d[,c("month", "bp", "subject")]) # add this to the full dataframe so we can add mean as a layer

        this.data <- d[d$month %in% c("start","end"),]
        direction <- as.factor(ifelse(this.data[this.data$month == "end", "bp"] - this.data[this.data$month == "start", "bp"] > 0, "increased", "decreased"))
        this.data <- cbind(this.data, direction=rep(direction, times=2))

        # add stats           
        d.bp <- this.data[this.data$month == "end", "bp"] - this.data[this.data$month == "start", "bp"] 
        print(shapiro.test(d.bp))
        pval <- signif(t.test(d.bp)$p.value,2) # one sample t.test with diffs
      
        this.data$month <- factor(this.data$month, levels=c("start","end"))
        this.data$subject <- as.factor(this.data$subject)
        cols <- alpha(c("#63C9D5","#AA6884"),.7)
        shapes <- c(21,16)
        names(cols) <- c("increased","decreased")
        names(shapes) <- c("start","end")
        
        p <- ggplot(data=this.data, aes(x=month, y=bp, group=subject, color = direction, shape=month)) + xlab("") +
            scale_color_manual(name="Bacteroides", values=cols) + scale_shape_manual(values=shapes, guide = 'none') +
            geom_line() + geom_point(size=4, fill="white") + 
            theme(legend.title = element_text(size=10, face="italic"), legend.text = element_text(size=8))
    } 
    else # if no clipping of months, plot all months
    {    
        cols <- alpha(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(d$subject))),.3)    
        p <- ggplot(data=d[d$month %in% 1:6,], aes(x=month, y=bp, color = subject)) + xlab("Month in the US") + scale_color_manual(values=cols) +
        geom_point(size=4) + theme(legend.position='none') + geom_line() + scale_x_continuous(breaks=c(1:6), labels=c(1:6))
        
        # use Robin's spline code to calculate p-value
        pval <- trendyspliner(data=d[d$month %in% 1:6,], xvar="month", yvar="bp", cases="subject", perm=99, quiet=T)$pval
        
    }
    p <- p + ggtitle("Bacteroides-Prevotella Ratio") + 
        ylab("log10(B/P)") +
        theme(legend.key = element_blank()) # removes borders from legend    
    p <- ggdraw(p) + draw_figure_label(label=paste0("p=",pval,"\n"), size=8, position="bottom.right")                

    save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
    
    invisible(p)
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

plot.b.p.ratio.all <- function(map0, otu, bug1, bug2, outputfn, g1, g2, g3)
{
    # don't include any samples that have completely zero of any because it screws up the log10 transform
    otu <- otu[otu[,bug1]!=0,]
    otu <- otu[otu[,bug2]!=0,]
    
    validsamples <- intersect(rownames(map0),rownames(otu))
    otu <- otu[validsamples,]
    map0 <- map0[validsamples,]
    
    bp <- get.log.taxa.ratio(otu, bug1, bug2)

    d <- data.frame(map0, y = bp)
    ylim <- range(d$y)
    d.1st <- d[d$Years.in.US > 0 & d$Years.in.US < 50,]
    d.2nd.control <- d[d$Years.in.US %in% c(50,60),]
    d.2nd.control$Sample.Group <- factor(d.2nd.control$Sample.Group)
    d.thai <-  d[d$Years.in.US == 0,]
    d.thai$Sample.Group <- factor(d.thai$Sample.Group)
    
    groupnames <- as.character(unique(d$Sample.Group))

    cols <- get.group.colors(groups=groupnames, alpha.val=1) 
    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 

    main <- ggplot(data=d.1st, aes(x=Years.in.US,y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, size=Sample.Group, alpha=Sample.Group), stroke=1, fill=NA) +
        geom_smooth(data=d.1st, aes(x=Years.in.US,y=y), color="black", size=.5) +
        scale_color_manual(name="Groups", values=cols) + 
        scale_alpha_manual(name="Groups", values=alphas) +
        scale_shape_manual(name="Groups", values=shapes) +
        scale_size_manual(name="Groups", values=sizes) + 
        ylim(ylim)  +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        ylab("") + xlab("Years in the US") + theme(legend.position='none') + 
        # remove all y axis
        theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.x = element_text(size=10), axis.text.x = element_text(size=8))

    left <- ggplot(d.thai, aes(Sample.Group, y, color=Sample.Group)) + geom_quasirandom(dodge.width=.75) + geom_boxplot(alpha=0, colour="black") +  
            scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("log10(Bacteroides-Prevotella Ratio)") + xlab("") +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),axis.text.x = element_text(size=8)) + scale_x_discrete(labels=c("HmongThai" = "Hmong\nThai", "KarenThai" = "Karen\nThai"))
    
    right <- ggplot(d.2nd.control, aes(Sample.Group, y, color=Sample.Group)) + 
            geom_quasirandom(dodge.width=.75, aes(color=Sample.Group, shape=Sample.Group, size=Sample.Group, alpha=Sample.Group)) +
            geom_boxplot(alpha=0, colour="black") +  
            scale_color_manual(name = "Groups", values = cols) + # set color for points from quasirandom
            scale_shape_manual(name="Groups", values=shapes) +
            scale_size_manual(name="Groups", values=sizes) + 
            scale_alpha_manual(name="Groups", values=alphas) +
            guides(fill=FALSE) + theme(legend.position='none') + ylim(ylim) +
            ylab("") + xlab("")  +
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
            # remove all y axis
            theme(axis.line.y=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.text.x = element_text(size=8))

    plots <- plot_grid(left, main, right, nrow=1, rel_widths=c(1.3,2,1), align="h") 
    title <- ggdraw() + draw_label("Bacteroides-Prevotella Ratio", fontface='bold')
    final <- plot_grid(title, plots, ncol=1, rel_heights=c(0.1, 1))
    save_plot(outputfn, final, useDingbats=FALSE, base_aspect_ratio = 1.5)
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