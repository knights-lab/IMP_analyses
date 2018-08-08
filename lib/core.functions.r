# plots 1 or more boxplots by BMI within each sample group against some arbitrary y-value
# e.g. used for plotting alpha diversity across groups by bmi class

require(beeswarm)
require(faraway)
require(RColorBrewer)
require(ggsignif)
require(ggbeeswarm)
require(ggplot2)
require(cowplot)

multiplot.boxplot.by.group.x.bmi <- function(map00, y.list, ylabs, mains, outputfn, parametric=TRUE, add.pvals=FALSE)
{
    p <- list()
    for(i in 1:length(ylabs))
    {
      p[[i]] <- plot.boxplot.by.group.x.bmi(map00=map00, y=y.list[[i]], ylab=ylabs[i], main=mains[i], parametric=parametric, add.pvals=add.pvals)
      # save legend (can prob do this just once)
      leg <- get_legend(p[[i]])
      # remove legends from all plots
      p[[i]] <- p[[i]] + theme(legend.position='none')
    }
    save_plot(gsub("pdf","legend\\.pdf",outputfn), leg, base_aspect_ratio=1)
    multiplot <- plot_grid(plotlist=p, nrow=1, ncol = length(mains))
    save_plot(outputfn, multiplot, ncol = length(mains), nrow = 1, base_aspect_ratio = 1)
}

# useful for annotations for geom_signif
get.signif.symbol <- function(pval, use.stars=FALSE)
{
    if(use.stars)
    {
        if(pval < .001)
            ret <- "***"
        else if (pval < .01)
            ret <- "**"
        else if (pval < .05)
            ret <- "*"
        else 
            ret <- "NS"
    }   
    else
    {
        ret <- paste0("p=", signif(pval, 2))
    }
    return(ret)
}

plot.boxplot.by.group.x.bmi <- function(map00, y, ylab, main, parametric=TRUE, add.pvals=FALSE) # vector y = variable on y axis 
{
    map0 <- data.frame(map00, y=y, stringsAsFactors=F)
    
    d <- map0[,c("y", "BMI.Class")]
    d[map0$Sample.Group %in% c("KarenThai", "HmongThai"), "Group"] <- "Thai"
    d[map0$Sample.Group %in% c("Karen1st", "Hmong1st"), "Group"] <- "1stGen"
    d[map0$Sample.Group == "Hmong2nd", "Group"] <- "2ndGen"
    d[map0$Sample.Group == "Control", "Group"] <- "Control"
    
    d$Group <- factor(d$Group, levels=c("Thai","1stGen","2ndGen","Control"))
    d$Group <- factor(d$Group) # remove any levels that aren't present
    d$BMI.Class <- factor(d$BMI.Class) # remove any levels that aren't present
        
    bmi.cols <- get.bmi.colors(alpha=.8)
    bmi.cols["Lean"] <- alpha(bmi.cols["Lean"], .5)
    
    d$pairwise.group <- as.factor(paste0(d$Group,".",d$BMI.Class))
    
    p <- ggplot(d, aes(Group, y, color=BMI.Class)) + geom_quasirandom(dodge.width=.75) + geom_boxplot(aes(fill = BMI.Class), alpha=0, colour="black") +  
      scale_color_manual(name = "BMI Class", values = bmi.cols) + # set color for points from quasirandom
        ggtitle(main) + 
        guides(fill=FALSE) + # turns off extra legend we get for boxplots created by BMI.Class
        ylab(ylab) + xlab("") + theme(axis.text.x = element_text(size=10), legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=10)) + 
        guides(colour = guide_legend(override.aes = list(size=3))) # make points in legend larger

    if(add.pvals)
        p <- add.pvals.to.plot(p, model=aov(d$y ~ d$Group), group.var="d$Group", alpha=.10, textsize=2, stepinc=.05)

    return(p)
}

plot.line.by.group.x.bmi <- function(map00, y, ylab, main, parametric=TRUE) # vector y = variable on y axis 
{
    map0 <- data.frame(map00, y=y, stringsAsFactors=F)
    
    d <- map0[,c("y", "BMI.Class")]
    d[map0$Sample.Group %in% c("KarenThai", "HmongThai"), "Group"] <- "Thai"
    d[map0$Sample.Group %in% c("Karen1st", "Hmong1st"), "Group"] <- "1st-Gen"
    d[map0$Sample.Group == "Hmong2nd", "Group"] <- "2nd-Gen"
    d[map0$Sample.Group == "Control", "Group"] <- "Control"
    
    d$Group <- factor(d$Group, levels=c("Thai","1st-Gen","2nd-Gen","Control"))
    d$Group <- factor(d$Group) # remove any levels that aren't present
    d$BMI.Class <- factor(d$BMI.Class) # remove any levels that aren't present
        
    bmi.cols <- get.bmi.colors(alpha=.7)
        
    p <- ggplot(d, aes(Group, y, color=BMI.Class, group=BMI.Class)) + geom_point(size=4) +
        #+ geom_quasirandom(dodge.width=.75) + geom_boxplot(aes(fill = BMI.Class), alpha=0, colour="black") +  
        scale_color_manual(name = "BMI Class", values = bmi.cols) + # set color for points from quasirandom
        ggtitle(main) + 
        guides(fill=FALSE) + # turns off extra legend we get for boxplots created by BMI.Class
        ylab(ylab) + xlab("") + theme(axis.text.x = element_text(size=10), legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=10)) + # make legend text smaller
        guides(colour = guide_legend(override.aes = list(size=3))) # make points in legend larger
    

    return(list(p=p, stats.label=""))
}


# plots longitudinal samples by some Y response variable, allowing for limiting to only first and last points
# num.clip.months = number of months to clip and average at the start and end
# y = response variable to plot, ordered by samples in map0
plot.response.L <- function(map0, y, outputfn=NULL, ylab, ggtitle, num.clip.months=NULL, show.stats=TRUE)
{        
    p <- NULL
    d <- data.frame(month = map0$Sample.Order, y = y, subject = map0$Subject.ID)
    if(!is.null(num.clip.months)) # take mean of arbitrary num of start and last months
    {
        d$subject <- as.character(d$subject)

        # find first x months, and average their response
        mins <- NULL
        for(m in 0:(num.clip.months-1))
            mins <- rbind(mins, aggregate(d$month, list(d$subject), FUN=function(xx) min(xx)+m))
        colnames(mins) <- c("subject","month")            
        mins.d <- merge(mins, d, by=c("subject","month"))
        mean.bp1 <- cbind(aggregate(mins.d$y, list(mins.d$subject), FUN=mean), month="start")

        # find last x months, and average their response
        maxs <- NULL
        for(m in 0:(num.clip.months-1))
            maxs <- rbind(maxs, aggregate(d$month, list(d$subject), FUN=function(xx) max(xx)-m))
        colnames(maxs) <- c("subject","month")
        maxs.d <- merge(maxs, d, by=c("subject","month"))
        mean.bp2 <- cbind(aggregate(maxs.d$y, list(maxs.d$subject), FUN=mean), month="end")

        mean.d <- rbind(mean.bp1, mean.bp2)
        colnames(mean.d) <- c("subject", "y", "month")
        d <- rbind(d, mean.d[,c("month", "y", "subject")]) # add this to the full dataframe so we can add mean as a layer

        this.data <- d[d$month %in% c("start","end"),]
        direction <- as.factor(ifelse(this.data[this.data$month == "end", "y"] - this.data[this.data$month == "start", "y"] > 0, "increased", "decreased"))
        this.data <- cbind(this.data, direction=rep(direction, times=2))

        # add stats           
        d.y <- this.data[this.data$month == "end", "y"] - this.data[this.data$month == "start", "y"] 
        print(shapiro.test(d.y))
        pval <- signif(t.test(d.y)$p.value,2) # one sample t.test with diffs
      
        this.data$month <- factor(this.data$month, levels=c("start","end"))
        this.data$subject <- as.factor(this.data$subject)
        fill.color= "black" #"#3A3096" #"#40004B"
        low.fill.color="#FF6600" #"#FFC900"
        cols <- c( alpha(fill.color,.5), alpha(low.fill.color,.8))
        shapes <- c(21,16)
        names(cols) <- c("increased","decreased")
        names(shapes) <- c("start","end")
        
        p <- ggplot(data=this.data, aes(x=month, y=y, group=subject, color = direction, shape=month)) + xlab("") +
            scale_color_manual(name="", values=cols) + scale_shape_manual(values=shapes, guide = 'none') +
            geom_line(size=1) + geom_point(size=4, fill="white", stroke=1, show.legend=FALSE) + theme(legend.position="bottom") + ggtitle(ggtitle)

    } 
    else # if no clipping of months, plot all months
    {    
        cols <- alpha(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(d$subject))),.3)    
        p <- ggplot(data=d, aes(x=month, y=y, color = subject)) + xlab("Month in the US") + scale_color_manual(values=cols) +
        geom_point(size=4) + theme(legend.position='none') + geom_line() + scale_x_continuous(breaks=c(min(d$month):max(d$month)), labels=c(min(d$month):max(d$month)))
        
        # use Robin's spline code to calculate p-value
        pval <- trendyspliner(data=d, xvar="month", yvar="y", cases="subject", perm=99, quiet=T)$pval
        
    }
    p <- p + ylab(ylab) 
        
    if(show.stats)
        # for showing P-values as title #p <- p + ggtitle(paste0("P=",pval)) + theme(plot.title = element_text(face="plain"))
        p <- ggdraw(p) + draw_figure_label(label=paste0("P=",pval,"\n"), size=8, position="bottom.right")

    if(is.null(outputfn))
        return(list(p=p, pval=pval))
    else
        save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
    

}
