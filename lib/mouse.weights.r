calc.weights <- function(ggdata)
{

    # for any map that didn't survive until Week 8, remove them
    ggdata$Mouse.ID <- as.character(ggdata$Mouse.ID)    
    valid.mouse.ids <- unique(ggdata[ggdata$Week==8, "Mouse.ID"])
    ggdata <- ggdata[ggdata$Mouse.ID %in% valid.mouse.ids,]

    mouse.ids <- unique(ggdata$Mouse.ID)
    for(i in 1:length(mouse.ids))
    {
        baseline.weight <- ggdata[ggdata$Week == 0 & ggdata$Mouse.ID == mouse.ids[i], "Body.Weight"]
        # % change
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Grams.Weight.Gain"] <- ggdata[ggdata$Mouse.ID == mouse.ids[i],"Body.Weight"] - baseline.weight
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Percent.Weight.Gain"] <- 100*(ggdata[ggdata$Mouse.ID == mouse.ids[i],"Grams.Weight.Gain"])/baseline.weight
    }

    return(ggdata)
    
}

plot.weights <- function(ggdata, group)
{
    ggdata$Group <- ggdata[,group]
    weight.all <- calc.weights(ggdata)

    weight.all.mean <- aggregate(Percent.Weight.Gain ~ Group + Week, weight.all, mean)
    weight.all.se <- aggregate(Percent.Weight.Gain ~ Group + Week, weight.all, FUN=function(xx) sd(xx)/sqrt(length(xx)))
    weight.all.df <- cbind(weight.all.mean[,1:2], Mean=weight.all.mean$Percent.Weight.Gain, SE=weight.all.se$Percent.Weight.Gain)

    p.weight.all <- ggplot(weight.all.df, aes(x=Week, y=Mean, group=Group, color=Group)) + 
            scale_color_manual(name = "", values = GROUP.COLORS) + 
            ylab("Percent Weight Gain") + xlab("Week") + scale_x_continuous(name="Week", breaks=seq(0,8,2)) + 
            geom_line(alpha=.7, size=2) + 
            geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2)
    p.legend <- get_legend(p.weight.all)
    p.weight.all <- p.weight.all + theme(legend.position='none')
    save_plot("mouse-weights-all.pdf", p.weight.all, base_aspect_ratio = 1.3)

    gg_weightwk8 <- ggdata[ggdata$Mouse.ID != "M21" & weight.all$Week %in% c(8,10),]
    # outlier M21 gained 100% of its weight - I think because it was a super tiny mouse and only 4-weeks old 

    p.weight.wk8 <- mouse.boxplot(y=gg_weightwk8$Percent.Weight.Gain, Group=gg_weightwk8$Group, main="% Weight Gain, Week 8", facet.var=NULL, 
                add.pval=TRUE, ylab="")

    save_plot("mouse-weights-endpoint.pdf", p.weight.wk8, base_aspect_ratio = 1.3)

}
    