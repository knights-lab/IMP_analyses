calc.vc <- function(vc, type, return.summary=TRUE)
{
    villus <- vc[vc$Type==type,]
    villus <- aggregate(Measurement ~ Mouse.ID + Group.End, villus, FUN=mean)
    v.mean <- aggregate(Measurement ~ Group.End, villus, FUN=mean)
    v.sd <- aggregate(Measurement ~ Group.End, villus, FUN=sd)
    v <- merge(v.mean, v.sd, by=1)
    colnames(v) <- c("Group", "Mean", "SD")
    if(return.summary) return(v)
    else return(villus)
}   

plot.vc <- function(villus.crypt.df, add.pval=TRUE, outputfn)
{    
    ratio <- villus.crypt.df[villus.crypt.df$Type=="villus (um)", ]
    ratio$Measurement <- villus.crypt.df[villus.crypt.df$Type=="villus (um)", "Measurement"] / villus.crypt.df[villus.crypt.df$Type=="crypt (um)", "Measurement"]
    ratio$Type <- "VC Ratio"

    vc <- data.frame(calc.vc(ratio, "VC Ratio", return.summary=FALSE), type="VC", stringsAsFactors=F)    
    vc <- merge(vc, unique(villus.crypt.df[,c("Mouse.ID","Donor.Type","Diet.Type")]),by="Mouse.ID")
    p <- mouse.boxplot(y=vc$Measurement, Group=vc$Group.End, main="Villus-Crypt Ratio", facet.var=NULL, add.pval=add.pval, ylab="", y.size=10, group.vars.df=vc[,c("Donor.Type","Diet.Type")])
    save_plot(outputfn, p, base_aspect_ratio = 1)

}

plot.v.or.c <- function(villus.crypt.df, type="v", add.pval=TRUE, outputfn)
{    
    if(type=="v")
    {
        Type <- "villus (um)"
        main <- "Villus Height"
    }
    else if (type=="c")
    {
        Type <- "crypt (um)"
        main <- "Crypt Depth"
    }
    ratio <- villus.crypt.df[villus.crypt.df$Type==Type, ]

    vc <- data.frame(calc.vc(ratio, Type, return.summary=FALSE), stringsAsFactors=F)    
    vc <- merge(vc, unique(villus.crypt.df[,c("Mouse.ID","Donor.Type","Diet.Type")]),by="Mouse.ID")
    
    p <- mouse.boxplot(y=vc$Measurement, Group=vc$Group.End, main=main, facet.var=NULL, add.pval=add.pval, ylab="(um)", y.size=10, group.vars.df=vc[,c("Donor.Type","Diet.Type")])
    save_plot(outputfn, p, base_aspect_ratio = 1)

}