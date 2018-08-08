
# d = columns as cell types, which will be plotted in separate figures
# colnames.short = formatted names for cell types
# group.vars = one or more variables for group vars
plot.immune <- function(d, group.var, named.cells, outputfn, facet.var=NULL, add.pval=FALSE, grouping.vars)
{
    p <- list()
    for(i in 1:length(named.cells))
    {
        p[[i]] <- mouse.boxplot(y=d[,named.cells[i]], Group=d[,group.var], main=names(named.cells[i]), facet.var=facet.var, add.pval=add.pval, alpha=.05, group.vars.df=d[,grouping.vars], title.size=20)
    }
    # save legend from a plot that has all groups
    p.legend <- mouse.boxplot(y=d[,named.cells[1]], Group=d[,group.var], main=names(named.cells[1]), facet.var=facet.var, add.pval=add.pval, plot.legend.only=TRUE, group.vars.df=d[,grouping.vars])

    if(is.null(outputfn))
    {
        return(p)
    }
    else
    {
        save_plot(paste0(dirname(outputfn), "/", "legend.pdf"), p.legend)
        save_plot(outputfn, plot_grid(plotlist=p, ncol=3), base_height=6, base_aspect_ratio = 1.5)
    }
}

# hand selected
make.immune.plots.manuscript <- function(ggdata, immune, immune_apc, group.var, grouping.vars, add.pval=FALSE, outdir, selected.cells, cell.pop="IEL")
{
    ggdata$Group <- ggdata[,group.var]
    apc <- c(grep("^Percent", colnames(immune_apc), value=T), grep("^Count", colnames(immune_apc), value=T))
    immune.df <- merge(merge(immune, ggdata, by="Mouse.ID"), immune_apc[,c(apc,"Sample")], by="Sample", all.x=T)
       
    cell.percents <- grep("^Percent", colnames(immune.df), value=T)
    names(cell.percents) <- gsub("Percent\\.[Live\\.]*", "", cell.percents)

    cell.counts <- grep("^Count", colnames(immune.df), value=T)
    names(cell.counts) <-  gsub("Count\\.[Live\\.]*", "", cell.counts)

    cd45.cell.percents <- cell.percents
    cd45.cell.percents <- cd45.cell.percents[names(cd45.cell.percents)!="CD45+"]  # remove CD45 itself from the list  
    cd45.cell.fnames <- names(cd45.cell.percents) # formatted names for CD45 expressing cells (all of them)
    names(cd45.cell.fnames) <- cd45.cell.fnames

    CELLPOP <- immune.df[immune.df$Cell.Population==cell.pop,]
    CELLPOP45 <- CELLPOP
    CELLPOP45[cd45.cell.fnames] <- CELLPOP45[,cd45.cell.percents]/CELLPOP45[,cell.percents["CD45+"]]
    p <- plot.immune(CELLPOP45, named.cells=cd45.cell.fnames[selected.cells], group.var=group.var, outputfn=NULL, add.pval=add.pval, grouping.vars=grouping.vars)
    CELLPOP[,cell.percents["CD45+"]] <- CELLPOP[,cell.percents["CD45+"]]
    p.cd45 <- plot.immune(d=CELLPOP, group.var=group.var, named.cells=cell.percents["CD45+"], outputfn=NULL, add.pval=add.pval, grouping.vars=grouping.vars)
    
    return(c(p.cd45,p))
}

make.immune.plots <- function(ggdata, immune, immune_apc, group.var, grouping.vars, add.pval=FALSE, outdir)
{
    ggdata$Group <- ggdata[,group.var]
    # add APCs
    apc <- c(grep("^Percent", colnames(immune_apc), value=T), grep("^Count", colnames(immune_apc), value=T))
    immune.df <- merge(merge(immune, ggdata, by="Mouse.ID"), immune_apc[,c(apc,"Sample")], by="Sample", all.x=T)
    
    
    cell.percents <- grep("^Percent", colnames(immune.df), value=T)
    names(cell.percents) <- gsub("Percent\\.[Live\\.]*", "", cell.percents)

    cell.counts <- grep("^Count", colnames(immune.df), value=T)
    names(cell.counts) <-  gsub("Count\\.[Live\\.]*", "", cell.counts)

    cd45.cell.percents <- cell.percents
    cd45.cell.percents <- cd45.cell.percents[names(cd45.cell.percents)!="CD45+"]  # remove CD45 itself from the list  
    cd45.cell.fnames <- names(cd45.cell.percents) # formatted names for CD45 expressing cells (all of them)
    names(cd45.cell.fnames) <- cd45.cell.fnames

    cd3e.cell.percents <- grep("Percent", colnames(immune), value=T) # APCs do not express CD3e (use immune table only)
    cd3e.cell.percents <- cd3e.cell.percents[cd3e.cell.percents!="Percent.CD3e+"] # remove cd3e itself
    cd3e.cell.fnames <- gsub("Percent\\.[Live\\.]*", "", cd3e.cell.percents) # formatted names for CD3e expressing cells
    names(cd3e.cell.fnames) <- cd3e.cell.fnames

    IEL <- immune.df[immune.df$Cell.Population=="IEL",]
    LPL <- immune.df[immune.df$Cell.Population=="LPL",]
    
    if(outdir != "./") dir.create(file.path(outdir), showWarnings = FALSE)

    IEL45 <- IEL
    IEL45[cd45.cell.fnames] <- IEL45[,cd45.cell.percents]/IEL45[,cell.percents["CD45+"]]

    plot.immune(IEL45, named.cells=cd45.cell.fnames, group.var=group.var, outputfn=paste0(outdir,"IEL-per-CD45.pdf"),add.pval=add.pval, grouping.vars=grouping.vars)

    LPL45 <- LPL
    LPL45[cd45.cell.fnames] <- LPL45[,cd45.cell.percents]/LPL45[,cell.percents["CD45+"]]
    plot.immune(LPL45, named.cells=cd45.cell.fnames, group.var=group.var, outputfn=paste0(outdir,"LPL-per-CD45.pdf"),add.pval=add.pval, grouping.vars=grouping.vars )

    IEL3e <- IEL
    IEL3e[cd3e.cell.fnames] <- IEL3e[,cd3e.cell.percents]/IEL3e[,cell.percents["CD3e+"]]
    plot.immune(IEL3e, named.cells=cd3e.cell.fnames, group.var=group.var, outputfn=paste0(outdir,"IEL-per-CD3e.pdf"),add.pval=add.pval, grouping.vars=grouping.vars)

    LPL3e <- LPL
    LPL3e[cd3e.cell.fnames] <- LPL3e[,cd3e.cell.percents]/LPL3e[,cell.percents["CD3e+"]]
    plot.immune(LPL3e, named.cells=cd3e.cell.fnames, group.var=group.var, outputfn=paste0(outdir,"LPL-per-CD3e.pdf"),add.pval=add.pval, grouping.vars=grouping.vars)

    # as counts per square inch
    plot.immune(d=IEL, group.var=group.var, named.cells=cell.counts, outputfn=paste0(outdir,"IEL-per-inch.pdf"), add.pval=add.pval, grouping.vars=grouping.vars)
    plot.immune(d=LPL, group.var=group.var, named.cells=cell.counts, outputfn=paste0(outdir,"LPL-per-inch.pdf"), add.pval=add.pval, grouping.vars=grouping.vars)

    # as percentages of total cells 
    plot.immune(d=IEL, group.var=group.var, named.cells=cell.percents, outputfn=paste0(outdir,"IEL-percent-totalcells.pdf"), add.pval=add.pval, grouping.vars=grouping.vars)
    plot.immune(d=LPL, group.var=group.var, named.cells=cell.percents, outputfn=paste0(outdir,"LPL-percent-totalcells.pdf"), add.pval=add.pval, grouping.vars=grouping.vars)

}