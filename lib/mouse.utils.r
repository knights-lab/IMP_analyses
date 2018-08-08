


# y = response , Group is grouping, group.vars = variables that determine groupings, facet.var = additional factor to facet by
mouse.boxplot <- function(y, Group, main, facet.var=NULL, alpha=.1, add.pval=FALSE, plot.legend.only=FALSE, ylab="", strip.text.size=5, title.size=14, y.size=11, x.size=11, group.vars.df, hide.box=TRUE)
{
    p <- map.boxplot.base(y=y, Group=Group, main=main, facet.var=facet.var, alpha=alpha, add.pval=add.pval, plot.legend.only=plot.legend.only, ylab=ylab, 
                        strip.text.size=strip.text.size, title.size=title.size, y.size=y.size, x.size=x.size,
                        hide.box=hide.box, fills=GROUP.FILLS, cols=GROUP.COLORS, alphas=GROUP.ALPHAS, shapes=GROUP.SHAPES, legend.labels = GROUP.NAMES.SHORT[levels(Group)], xlabels=GROUP.NAMES.SHORT[levels(Group)],
                        group.vars.df=group.vars.df)                    
    return(p)
}

# 
# mouse.boxplot <- function(y, Group, main, facet.var=NULL, alpha=.1, add.pval=FALSE, plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=10, x.size=9)
# {
#     ggdata <- data.frame(y=y, Group=Group)
#     if(!is.null(facet.var)) ggdata$facet.var <- facet.var
# 
# 
#     p <- ggplot(ggdata, aes(x=Group, y=y, group=Group)) + 
#         geom_boxplot(aes(fill=Group, alpha=Group), colour="black", width=.9) + # unfilled boxplot
#         geom_point(size=2, aes(shape=Group), stroke=1) + #, color="black") +
#         scale_color_manual(name="Groups", values=GROUP.COLORS) + 
#         scale_fill_manual(name="Groups", values=GROUP.FILLS) + 
#         scale_alpha_manual(name="Groups", values=GROUP.ALPHAS) +
#         scale_shape_manual(name="Groups", values=GROUP.SHAPES) +
#         ggtitle(main) + ylab(ylab) + xlab("") + 
#         scale_x_discrete(labels=GROUP.NAMES.SHORT[levels(ggdata$Group)]) +
#         theme(axis.text.y = element_text(size=y.size), axis.text.x = element_text(size=x.size),
#         #axis.text.x = element_blank(), axis.ticks.x=element_blank(), 
#         axis.title.x=element_blank(), axis.title.y=element_text(size=11),
#         strip.background = element_blank(),
#         strip.text.x = element_text(size = strip.text.size), 
#         plot.title = element_text(size=12))
#     
#     if(plot.legend.only)
#       p <- get_legend(p)
#     else
#     {
#         p <- p + theme(legend.position="none")    
#         if(!is.null(facet.var))
#             p <- p + facet_grid(. ~ facet.var)
#         if(add.pval) # only compare non-cohoused groups
#             p <- add.pvals.to.plot(p=p, model=aov(ggdata$y ~ ggdata$Group), group.var="ggdata$Group", alpha=alpha)
#     }
#     return(p)
# }

