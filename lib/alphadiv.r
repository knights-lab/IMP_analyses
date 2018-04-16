multiplot.alphadiv.L <- function(map0, alpha_df, metrics, outputfn,x.var="Sample.Day.Since.Arrival", xlab="Days since arrival")
{
    p <- NULL
    for(i in 1:length(metrics))
    {
        p[[i]] <- plot.alphadiv.L(map0, alpha_df[rownames(map0), metrics[i]], metrics[i], x.var=x.var, xlab=xlab)
    }
    multiplot <- plot_grid(plotlist=p, ncol=(length(p)), nrow=1)
    save_plot(outputfn, multiplot, ncol = length(p), nrow = 1, base_aspect_ratio = 1.3)

}

# x.var = Sample.Day or Sample.Order
plot.alphadiv.L <- function(map0, alpha, main, x.var="Sample.Day.Since.Arrival", xlab="Days since arrival")
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