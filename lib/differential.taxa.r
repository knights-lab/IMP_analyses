
# finds taxa that are differentiated by x.var (can be factor or continuous) then
# plots only those that are significant
# control.vars = names of variables to control for
plot.diff.taxa <- function(map0, taxa0, x.var, control.vars=NULL, sig.level=.05, outputfn.prepend, do.sqrt=TRUE)
{
    valid.samples <- intersect(rownames(map0), rownames(taxa0))
    map0 <- map0[valid.samples,]
    taxa0 <- taxa0[valid.samples,]

    # filter taxa
    prevalences <- apply(taxa0, 2, function(bug.col) mean(bug.col > 0))
    taxa0 <- taxa0[, prevalences >= .10]

    # additional filtering to allow for parametric tests
    if(do.sqrt) taxa0 <- asin(sqrt(taxa0))
    ret <- collapse.by.correlation(taxa0, .95)
    taxa0 <- taxa0[, ret$reps]

    ret <- test.features.parametric(taxa0, map0[, x.var], controls=map0[,control.vars, drop=F], sig.level=sig.level)      
    sig.taxa <- ret$features
   
    if(length(sig.taxa)==0)
        print("No significant taxa found")        
    else
    {
        for(i in 1:length(sig.taxa))
        {
            taxa_abbrev <- shorten.taxonomy(sig.taxa[i])
        
            ggdata <- data.frame(y=taxa0[,sig.taxa[i]], x=map0[,x.var])
            
            if(is.factor(map0[,x.var]))
            {
                # don't add stats, we'll add our own since we're only plotting sig taxa anyway
                p <- plot.boxplot.by.group.x.bmi(map0=map0, y=ggdata$y, ylab="", main=taxa_abbrev, add.stats=FALSE)
            }
            else
            {
                p <- ggplot(data=ggdata, aes(x=x, y=y)) + geom_point(size=2, colour=alpha("black",.5)) +
                        ggtitle(label=taxa_abbrev) + theme(legend.position="none") + ylab("Relative Abundance") + xlab(gsub("\\.", " ", x.var))

                # stat_smooth smooths on all data points
                p <- p + stat_smooth(aes(group = 1), method = "lm", se = T, size=2)
            }

            # add stats to bottom right
            p <- ggdraw(p) + draw_figure_label(label=paste0("FDR-adjusted (number of taxa) q = ", signif(ret$pvals[sig.taxa[i]], 2)), size=8, position="bottom.right")

            save_plot(paste(outputfn.prepend, taxa_abbrev, "pdf", sep="."), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
        
        }
    }
    
}