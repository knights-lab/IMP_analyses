multiplot.alphadiv.L <- function(map0, alpha_df, metrics, outputfn, num.clip.months=NULL)
{
    plist <- NULL
    for(i in 1:length(metrics))
    {
        valid_names <- intersect(rownames(map0), rownames(alpha_df))
        p <- plot.response.L(map0[valid_names,], y=alpha_df[valid_names, metrics[i]], outputfn=NULL, ylab="", ggtitle=metrics[i], num.clip.months=num.clip.months)
        plist[[i]] <- p
    }
    multiplot <- plot_grid(plotlist=plist, ncol=(length(plist)), nrow=1)
    save_plot(outputfn, multiplot, ncol = length(plist), nrow = 1, base_aspect_ratio = 1)
}

# x.var = Sample.Day or Sample.Order
plot.alphadiv.L.old <- function(map0, alpha, main, x.var="Sample.Day.Since.Arrival", xlab="Days since arrival",  num.clip.months=NULL)
{
    d <- data.frame(x = map0[,x.var], y = alpha, id = map0[,"Subject.ID"])
    cols <- alpha(colorRampPalette(brewer.pal(9, "Set1"))(length(unique(d$id))),.3)    
    
    p <- ggplot(data=d, aes(x, y, color = id)) + geom_line() + xlab(xlab) + ylab("") + ggtitle(main) + theme(legend.position='none') +
        scale_color_manual(values=cols) + geom_point(size=4)
    
    # draw loess line
    p <- p + geom_smooth(aes(group = 1), size=2, se = FALSE, color=alpha("black",.7))
    # add robin's permutation based test
    pvaltext <- paste("P = " , trendyspliner(data=d, xvar="x", yvar="y", cases="id", perm=999, quiet=T)$pval)
        
    p <- ggdraw(p) + draw_label(pvaltext, size=8, x=.85, y=.18)

    return(p)
}
