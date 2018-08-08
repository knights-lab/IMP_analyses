plot.glucose <- function(ggdata, group, add.pval=TRUE, outputfn)
{
    ggdata$Group <- ggdata[,group]

#    ggdata$Week <- factor(as.character(ggdata$Week), labels=c("Baseline", "Post-Diet", "Post-Cohousing"), levels=c(0,8,10))
    ggdata$Week <- factor(as.character(ggdata$Week), labels=c("Baseline", "Post-Diet"), levels=c(0,8))

    p <- ggplot(ggdata, aes(x=Week, y=Fasting.Glucose, color=Group)) + geom_boxplot(aes(fill = Group), color="black") +
        facet_grid(. ~ Week, scales = "free", switch="x") +
        scale_fill_manual(name="", values = GROUP.COLORS) +
        ggtitle("Fasting Blood Glucose") + 
        ylab("(mg/dL)") + xlab("") +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank()) + theme(legend.position="none")
    save_plot("mouse-glucose-boxplot.pdf", p, base_aspect_ratio = 1.3)

    # plot by fold change
    base.blood <- ggdata[ggdata$Week=="Baseline",c("Mouse.ID","Fasting.Glucose", "Donor")]
    postdiet.blood <- ggdata[ggdata$Week=="Post-Diet",c("Mouse.ID","Fasting.Glucose","Week","Group","Cohoused")]
#    postcohouse.blood <- ggdata[ggdata$Week=="Post-Cohousing",c("Mouse.ID","Fasting.Glucose","Week","Group","Cohoused")]

    dfdiet <- merge(postdiet.blood, base.blood, by="Mouse.ID", suffixes=c(".post",".base"))    
 #   dfcohouse <- merge(postcohouse.blood, base.blood, by="Mouse.ID", suffixes=c(".post",".base"))
  #  fc.df <- rbind(dfdiet,dfcohouse)
    fc.df <- dfdiet
    fc.df$fold.change <- (fc.df$Fasting.Glucose.post - fc.df$Fasting.Glucose.base)/fc.df$Fasting.Glucose.base

    fc.df <- merge(fc.df, unique(ggdata[,c("Diet.Type","Donor.Type","Mouse.ID")]), by="Mouse.ID")

    p <- mouse.boxplot(y=fc.df$fold.change, Group=fc.df$Group, main="Fasting Blood Glucose Fold Change\n", add.pval=add.pval, ylab="(mg/dL)", group.vars.df=fc.df[,c("Diet.Type","Donor.Type")], hide.box=F)

    save_plot(outputfn, p, base_aspect_ratio = 1)
}