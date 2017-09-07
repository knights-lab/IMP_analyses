# plots 1 or more boxplots by BMI within each sample group against some arbitrary y-value
# e.g. used for plotting alpha diversity across groups by bmi class

require(beeswarm)
require(faraway)
require(RColorBrewer)
require(ggsignif)
require(ggbeeswarm)


multiplot.boxplot.by.group.x.bmi <- function(map00, y.list, ylabs, mains, outputfn, parametric=TRUE)
{
    sink(paste0(outputfn, ".txt")) # set text file for writing stats results
    cat("Non-parametric two-sample test, adjusted for number of comparisons (P-value and Adjusted P-value)\n\n")
    p <- NULL
    for(i in 1:length(ylabs))
    {
      cat(mains[i], "\n")
      p[[i]] <- plot.boxplot.by.group.x.bmi(map00=map00, y=y.list[[i]], ylab=ylabs[i], main=mains[i], parametric=parametric)
      # save legend (can prob do this just once)
      legend <- get_legend(p[[i]])
      # remove legends from all plots
      p[[i]] <- p[[i]] + theme(legend.position='none')
      cat("\n\n")
    }
    sink()
    
    multiplot <- plot_grid(plot_grid(plotlist=p, ncol=(length(p)), nrow=1), plot_grid(legend), rel_widths=c(1,.11))
    save_plot(outputfn, multiplot, ncol = 2, nrow = 1, base_aspect_ratio = 1.3)
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

plot.boxplot.by.group.x.bmi <- function(map00, y, ylab, main, add.stats=TRUE, parametric=TRUE) # vector y = variable on y axis 
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
        
    bmi.cols <- alpha(c("#31625b","#fcd167","#c45565"),.5)
    names(bmi.cols) <- c("Lean", "Overweight", "Obese")
    
    p <- ggplot(d, aes(Group, y, color=BMI.Class)) + geom_boxplot(aes(fill = BMI.Class), alpha=0, colour="black") + geom_quasirandom(dodge.width=.75) +  
      scale_color_manual(name = "BMI Class", values = bmi.cols) + # set color for points from quasirandom
        ggtitle(main) + 
        guides(fill=FALSE) + # turns off extra legend we get for boxplots created by BMI.Class
        ylab(ylab) + xlab("") + theme(axis.text.x = element_text(size=10), legend.title=element_text(size=10, face="bold"), legend.text=element_text(size=10)) + # make legend text smaller
        guides(colour = guide_legend(override.aes = list(size=3))) # make points in legend larger
    
    if(add.stats){
        # within sample groups, test for differences in BMI class
        groups <- levels(d$Group)
        pvals <- NULL
        for(i in 1:length(groups))
        {
            this.pvals <- test.groups(d[d$Group==groups[i], "y"], d[d$Group==groups[i], "BMI.Class"], parametric)
            names(this.pvals) <- paste(groups[i],names(this.pvals),sep=":")
            pvals <- c(pvals, this.pvals)
        }
        p <- ggdraw(p) + draw_figure_label(label=paste(names(pvals), signif(pvals, 2), sep="=", collapse="; "), size=8, position="bottom.right")
    }
    return(p)
}
