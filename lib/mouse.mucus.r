prep.colon.data <- function()
{
    setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/Mouse Experiments/Hunter Lab")

    new.map <- read.table("Mouse Colon Tracking.txt", sep="\t", header=T, stringsAsFactors=F)
    meas <-  read.table("Colon Mucus thickness Measurements.txt", sep="\t", header=T, stringsAsFactors=F)

    # found error - this should be fixed now
    # meas$Blinded.Mouse.Letter[meas$Blinded.Mouse.Letter=="B "] <- "B"

    # check that everything has exactly 15 measurements
    table(meas$Blinded.Mouse.Letter, meas$Fecal.Pellet.Number,meas$Fecal.Pellet.Section)

    final.df <- merge(meas, new.map[,c("Mouse.ID", "Letter.Assigned", "Reembedded")], by.x="Blinded.Mouse.Letter", by.y="Letter.Assigned")

    write.table(final.df, "mucus.txt", sep="\t", quote=F, col.names=T, row.names=F)

    # final dataset contains n=32 mice because n=3 died and never made it and we can't image n=3
}

normalize.mucus <- function(m)
{
    # normalize by researcher
    # x - min / max - min
    
    researchers <- unique(m$Initials)
    for(i in 1:length(researchers))
    {
        ix <- which(m$Initials == researchers[i])
        xx <- m[ix,"Measurement"]
        norm.xx <- (xx - min(xx)) / (max(xx) - min(xx))
        m[ix,"Normalized.Measurement"] <- norm.xx
    }
    
    # measurement for each image
    m.researcher <- aggregate(Normalized.Measurement ~ Fecal.Pellet.Number + Fecal.Pellet.Section + Mouse.ID, data=m, mean)

    # measurement for each pellet
    m.pellet <- aggregate(Normalized.Measurement ~ Fecal.Pellet.Number + Mouse.ID, data=m, mean)
 
    # measurement for each mouse
    m.mouse <- aggregate(Normalized.Measurement ~ Mouse.ID, data=m, mean)

    # let's return avg measurements for every image available    
    return(m.researcher)   
}

plot.mucus <- function(ggdata, add.pval=TRUE, highlight.reembedded=FALSE)
{
    if(highlight.reembedded)
    {
        main="Mucus Measurements\n"; add.pval=add.pval;
        ylab="um"; 
        #group.vars.df=ggdata[,c("Donor.Type","Diet.Type","Reembedded")]; 
        group.vars.df=NULL
        hide.box=FALSE

        facet.var=NULL; alpha=.1; plot.legend.only=F; 
        fills=GROUP.FILLS; cols=GROUP.COLORS; alphas=GROUP.ALPHAS; shapes=GROUP.SHAPES; legend.labels = GROUP.NAMES.SHORT[levels(ggdata$Group.End)];
        xlabels=GROUP.NAMES.SHORT[levels(ggdata$Group.End)];
        group.vars.df=group.vars.df; facet.var=NULL     
        strip.text.size=5; title.size=12; y.size=10; x.size=9

        cols.embed=c(Yes="red",No="black")

        d <- data.frame(y=ggdata$Normalized.Measurement, Group=ggdata$Group.End, Group2=ggdata$Reembedded)
        if(!is.null(group.vars.df)) d <- cbind(d, group.vars.df)
        if(!is.null(facet.var)) d$facet.var <- facet.var
        p <- ggplot(d, aes(x=Group, y=y, group=Group)) + 
           scale_fill_manual(name="Groups", values=fills, labels=legend.labels) + 
           scale_color_manual(name = "Group2", values = cols.embed, labels=legend.labels) +
            scale_alpha_manual(name="Groups", values=alphas, labels=legend.labels) +
            scale_shape_manual(name="Groups", values=shapes, labels=legend.labels) +
            ggtitle(main) + ylab(ylab) + xlab("") + 
            scale_x_discrete(labels=xlabels) + #
            theme(axis.text.y = element_text(size=y.size), axis.text.x = element_text(size=x.size),
            axis.title.x=element_blank(), axis.title.y=element_text(size=11),
            strip.background = element_blank(),
            strip.text.x = element_text(size = strip.text.size), 
            plot.title = element_text(size=title.size)) 

        p <- p + geom_boxplot(colour="black", width=.9, aes(fill=Group, alpha=Group)) + geom_quasirandom(dodge.width=.75, aes(shape=Group, color=Group2), size=2, stroke=1)

        p <- p + theme(legend.position="none")    

        if(add.pval){
            if(is.null(group.vars.df))
            {
                # plot with reembedding as an interaction to check if this matters
                print(summary(aov(d$y ~ d$Group * d$Group2)))
                p <- add.pvals.to.plot(p=p, model=aov(d$y ~ d$Group + d$Group2), group.var="d$Group", alpha=alpha)
            }
            else
            {
                f <- as.formula(paste0("y ~ ", paste(colnames(group.vars.df), collapse="+")))
                p <- add.pvals.to.plot(p=p, model=aov(f, data=d), group.var=colnames(group.vars.df), alpha=alpha, method="two-way")
            }
        }
    }
    else
    {
        p <- mouse.boxplot(y=ggdata$Normalized.Measurement, Group=ggdata$Group.End, 
            main="Mucus Measurements\n", add.pval=add.pval,
            ylab="um", group.vars.df=ggdata[,c("Donor.Type","Diet.Type")], hide.box=FALSE)

    }
    return(p)
}
