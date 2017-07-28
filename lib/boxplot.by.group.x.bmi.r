# plots 1 or more boxplots by BMI within each sample group against some arbitrary y-value
# e.g. used for plotting alpha diversity across groups by bmi class

require(beeswarm)
require(faraway)
require(RColorBrewer)

# x = +2-level factor, y = continuous value
# simply runs wilcoxon test, OR manually loops through all possible combinations and corrects p-values for 3+ level group factors
test.samples <- function(y, x)
{
    comparisons <- combn(levels(x), 2)
    
    pvals <- NULL
    for(i in 1:ncol(comparisons))
    {

        this.y <- y[x %in% comparisons[,i]]
        this.x <- x[x %in% comparisons[,i]]
        m <- wilcox.test(this.y ~ this.x)
        pvals <- c(pvals, m$p.value)
        
        print(paste0(comparisons[1,i], " VS ", comparisons[2,i], ": ", m$p.value))
    }
    return(p.adjust(pvals, method="fdr"))

}


multiplot.boxplot.by.group.x.bmi <- function(map00, y.list, ylabs, mains, outputfn)
{
    p <- NULL
    for(i in 1:length(ylabs))
    {
        p[[i]] <- plot.boxplot.by.group.x.bmi(map00, y.list[[i]], ylabs[i], mains[i])
    }
    multiplot <- plot_grid(plotlist=p, ncol=length(p), nrow=1)
    save_plot(outputfn, multiplot, ncol = length(p), nrow = 1, base_aspect_ratio = 1.3)
}

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

plot.boxplot.by.group.x.bmi <- function(map00, y, ylab, main) # vector y = variable on y axis 
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
    
    print(d$BMI.Class)
    
    # within sample groups, test for differences in BMI class
    groups <- levels(d$Group)

    pvals <- NULL
    for(i in 1:length(groups))
    {
        pvals <- c(pvals, test.samples(d[d$Group==groups[i], "y"], d[d$Group==groups[i], "BMI.Class"]))
    }    
    pval.labels <- sapply(pvals, get.signif.symbol)
    print(pval.labels)
    
    ymax <- max(d$y)*1.05 # add ymax to make space for pvalues
    p <- ggplot(d, aes(Group, y)) + geom_boxplot(aes(fill = BMI.Class)) +  geom_point(aes(y=y, fill = BMI.Class), position=position_dodge(width=.75), color=alpha("black",.3), shape=1) +
        scale_fill_manual(name = "Sample Groups", values = c("#ffffcc","#a1dab4","#41b6c4")) +      
        ggtitle(main) + coord_cartesian(ylim = c(min(d$y), ymax)) +
        ylab(ylab) + xlab("BMI Class") + theme(axis.text.x = element_text(size=10))

    # add pval annotations
    for(i in 0:(length(pval.labels)-1))        
        p <- p + geom_signif(annotation = pval.labels[i+1], xmin = 0.8+i, xmax = 1.2+i, y_position = ymax, tip_length=0)

    return(p)
}
