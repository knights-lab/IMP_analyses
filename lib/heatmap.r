# adapted from mwas.heatmap used for marty blaser's diabetes mice experiment
source(paste0(LIBDIR,"heatmap.3.r"))
source(paste0(LIBDIR,"collapse-features.r"))
source(paste0(LIBDIR,"utils.r"))

require(vegan)
require(RColorBrewer)
require(gplots) 

# samples are columns
# hclusts samples within each group
cluster.columns<-function(otu0, distfun, new.treatments)
{
	sample.labels <- NULL
	# separate samples by unique values of new.treatments that we're interested in clustering by
	unique.new.treatments <- levels(new.treatments)

	for(i in 1:length(unique.new.treatments)){
		otu0.by.clustergroup <- otu0[new.treatments==unique.new.treatments[i], ]
 		data.dist.group <- distfun(otu0.by.clustergroup)
 		col.clus.group <- hclust(data.dist.group)

		d <- as.dendrogram(col.clus.group)
		# since the index order is per group only, we'll have to add the current length each time 
		# to get correct indices along the combined dendrogram
		sample.labels <- c(sample.labels, labels(d))
	}
	sample.labels
	
}

# OTUs will be displayed as rows, make sure that otu0.t is the otu table with OTUs as rows
cluster.rows<-function(otu0.t, distfun)
{
	# do the row-clustering (note that we need to explicitly reorder the dendrogram here)
	Rowv <- rowMeans(otu0.t)
	hcr <- hclust(distfun(otu0.t))
	ddr <- as.dendrogram(hcr)
	ddr <- reorder(ddr, Rowv)
	rowInd <- order.dendrogram(ddr)
#		row.labels <- labels(ddr)[rowInd] #THIS IS CRITICAL?
	otu.labels <- labels(ddr)

	list(ddr=ddr, row.labels=otu.labels)
}

# map0: rows must have already been reorder by the hclust
# color.list: list of named vectors that have been named according to color.var values
# num.breaks for gradients
create.color.bars<-function(map0, color.var, color.list, num.breaks=20)
{
	col.colors <- NULL
	for(i in 1:length(color.var))
	{
		reordered <- map0[,color.var[i]]
		if(is.numeric(reordered)){
			color.palette <- colorRampPalette(color.list[[color.var[i]]])
			colors <- color.palette(num.breaks)[as.numeric(cut(reordered,breaks = num.breaks))]
		} else {
				colors <- unname(color.list[[color.var[i]]][as.character(reordered)])
		}
		col.colors <- cbind(col.colors, colors)		
	}
	colnames(col.colors) <- color.var
	col.colors
}

# combines cluster variable values to create new groupings to cluster by
create.new.treatments<-function(map0, cluster.var, color.list)
{
	if(length(cluster.var) > 1){
		new.treatments <- as.factor(apply(map0[,cluster.var], 1, paste0, collapse="."))
		names(new.treatments) <- rownames(map0)	
	} else {
        new.treatments <- map0[,cluster.var]
		names(new.treatments) <- rownames(map0)
		# cluster.var must be a factor in order to cluster properly, apply original levels
		new.treatments <- as.factor(new.treatments)
		# reorder the factor levels to reflect user-defined cluster var order
		new.treatments <- factor(new.treatments, levels = names(color.list[[cluster.var]]))
	}
	new.treatments
}

# gets list of colors based on cluster.var
# override this function to generate other colors and variables
get.color.list<-function(map0)
{
    color.list <- list()
    alphas <- get.group.alphas()
    cols <- get.group.colors()
    
    color.list[["Sample.Group"]] <- unlist(lapply(names(cols), FUN=function(xx) alpha(cols[xx], alphas[xx])))
    names(color.list[["Sample.Group"]]) <- names(cols)
    color.list[["Sample.Group"]] <- color.list[["Sample.Group"]][names(color.list[["Sample.Group"]]) %in% unique(map0$Sample.Group)]

    color.list[["Age"]] <- setNames(c('white','#CC9900'), as.character(signif(range(map0$Age),2)))
    
    color.list[["Years.in.US"]] <- setNames(c('white','#005073'), as.character(signif(range(map0$Years.in.US),2)))
    return(color.list)
}


# This heatmap converts OTUs to presence/absence, and plots OTUs from left to right (highest to lowest prevalence) and samples from top to bottom (non-western to western)
#   sorts OTUs by prevalence in baseline groups (min prevalence cut off also only applied at baseline)
#   sorts samples by richness
make.heatmap.binary <- function(otu0, map0, group.var="Sample.Group", min.prevalence=0.25, baseline.groups, show.colnames=FALSE, outputfn)
{
    valid.samples <- intersect(rownames(map0), rownames(otu0))
    map0 <- map0[valid.samples,]
    otu0 <- otu0[valid.samples,]

    map0[,group.var] <- factor(map0[,group.var]) # remove levels not present after filtering

    # sort columns (otus) by prevalence of baseline group only
    this.taxa <- otu0[map0[,group.var] %in% baseline.groups,]
    prevalences <- apply(this.taxa, 2, function(bug.col) mean(bug.col > 0))
    # filter prevalence within the baseline sample group
    prevalences <- prevalences[prevalences >= min.prevalence]
    # only add unique values to final ordered taxa list
    sorted.taxa <- names(sort(prevalences, decreasing=T))
    otu0 <- otu0[,sorted.taxa]
    # convert table to presence absence
    otu0[otu0 > 0] <- 1

    # sort rows by richness (number of otus in sample)
    final.otu <- NULL
    for(this.group in levels(map0[,group.var]))
    {
        this.otu <- otu0[map0[,group.var]==this.group,,drop=F]

        if(nrow(this.otu) > 1) # only reorder if we have more than 1 sample!
            this.otu <- this.otu[order(rowSums(this.otu)),]
            
        final.otu <- rbind(final.otu, this.otu)
    }
    otu0 <- final.otu

    ok.to.draw=TRUE
    if(sum(is.infinite(otu0)) != 0 | sum(is.na(otu0)) != 0)
    {   
        print("Found infinites or NAs in your OTU table")
        ok.to.draw <- FALSE
    }
    if(ncol(otu0) < 2 | nrow(otu0) < 2)
    {
        print("Only one differential feature found. Cannot generate heatmap.")
        ok.to.draw <- FALSE
    }

    if(ok.to.draw)
    {
        color.list <- get.color.list(map0)
    
        # don't cluster OTUs - should already be ordered by prevalence
        col.labels <- colnames(otu0)
        row.labels <- rownames(otu0)
            
        # draw line between groups in order to see separation better
        last.groups <- aggregate(1:length(row.labels), list(map0[row.labels,group.var]), max)
        last.group.ix <- as.numeric(as.character(last.groups[,2]))
        last.group.ix <- last.group.ix[-which.max(last.group.ix)]
    
        row.colors <- t(create.color.bars(map0[row.labels,,drop=F], color.var=names(color.list), color.list))
        row.colors <- row.colors[rev(rownames(row.colors)),] # reverse order so group is closest to heatmap
        otu.hm <- as.matrix(otu0)
        draw.prevalence.heatmap(otu.hm=otu.hm, color.list=color.list, row.colors=row.colors, 
                                row.sep.ix=last.group.ix, show.colnames=show.colnames, outputfile=outputfn)
    }
}

draw.prevalence.heatmap <- function(otu.hm, color.list, row.colors, row.sep.ix, show.colnames=FALSE, outputfile)
{
	hm.colors <- colorRampPalette(c("#f2ece3","#365B6D"))(2)    # white and dark greygreen

    # note that these settings are different from trad heatmap, and make a huge difference
    lwid = c(1.2,.2,9,.5) #col widths
	lmat = rbind(c(0,0,4,0), c(3,1,2,0), c(0,0,5,0))	
    lhei = c(1, 8, 1) #row heights
    margins=c(.05,.05)
    legend1.x=0
    legend1.y=.95
	
    width=30
    height=12
	
	if(show.colnames) labCol <- colnames(otu.hm)
	else labCol <- NA

	pdf(file=outputfile, width=width, height=height)	    
	
	heatmap.3(otu.hm, Colv = NA, Rowv = NA, dendrogram='none',
            col = hm.colors, 
            RowSideColors = row.colors,
            rowAxisColors=NULL,
            trace = "none", density.info = "none", 
            key=FALSE, keysize=5, 
            lmat=lmat, lwid=lwid, lhei=lhei, margins=margins, 
            labCol=labCol, labRow=NA, 
            cexCol=.7, scale="none",
            rowsep=row.sep.ix, sepwidth=c(.1,.1), sepcolor="white") # indices for where to draw row separators and width of line

	# add the legend for the colored bars
    legend(x=legend1.x, y=legend1.y,legend=unname(unlist(lapply(color.list, names))), 
    fill=unlist(color.list), border=FALSE, bty="n", cex=.9, ncol=2, text.width=strwidth("1,000,000"))
	
	dev.off()
}


# keeps prevalent OTUs, collapse correlated OTUs, transforms rel abundances if parametric testing, keeps only significant OTUs
prep.otu.trad.heatmap<-function(otu0, groups, min.prevalence=0.25, sig.level=.25, parametric=TRUE, normalize=TRUE)
{
    # if not normalizing (e.g. clr), do not convert to relative abundance, filter low prevalence, or sqrt transform
    if(normalize){ 
        otu0 <- sweep(otu0, 1, rowSums(otu0), '/')	
        prevalences <- apply(otu0, 2, function(bug.col) mean(bug.col > 0))
        otu0 <- otu0[, prevalences >= min.prevalence]
        otu0 <- asin(sqrt(otu0))
    }
    ret <- collapse.by.correlation(otu0, .95)
    otu0 <- otu0[, ret$reps]

    # select features that are statistically different between groups (parametric)
    if(parametric) ret <- test.features.parametric(otu0, groups, sig.level=sig.level)
    else ret <- test.features.nonparametric(otu0, groups, sig.level=sig.level)

	otu0 <- otu0[,ret$features,drop=F]	
    return(list(otu=otu0, ret=ret))
}


# This makes a traditional heatmap (samples as columns, non-western to western and OTUs as rows)
# min.prevalence applied globally 
make.heatmap.traditional <- function(otu0, map0, outputfn, group.var="Sample.Group", filter.mode="standard", save.pvals=FALSE, min.prev=0.1, sig.level=.10, show.rownames=FALSE)
{
    valid.samples <- intersect(rownames(map0), rownames(otu0))
    map0 <- map0[valid.samples,]
    otu0 <- otu0[valid.samples,]

    map0[,group.var] <- factor(map0[,group.var]) # remove levels not present
    
    # useful to see distribution of samples
    print(table(map0[,group.var]))
    
    color.list <- get.color.list(map0)
    print(color.list)
    
	new.treatments <- create.new.treatments(map0, cluster.var=group.var, color.list)

    # transforms and filters otu using standard methods
    if(filter.mode!="none")
    {
        if(filter.mode=="standard") 
            prep.out <- prep.otu.trad.heatmap(otu0, groups=new.treatments, min.prevalence=min.prev, sig.level=sig.level)
        else if (filter.mode=="clr")
            prep.out <- prep.otu.trad.heatmap(otu0, groups=new.treatments, min.prevalence=min.prev, sig.level=sig.level, normalize=FALSE)
        otu0 <-prep.out$otu
        ret <- prep.out$ret
        
        if(save.pvals)
        {
            d <- data.frame(adj.pval=ret$adj.pvals, pval=ret$pvals, pathway.id=gsub(":.*","",names(ret$adj.pvals)), row.names=names(ret$adj.pvals))
            d <- d[order(d$adj.pval),] 
            pval.outfn <- gsub("\\.pdf","\\.pvals\\.txt",outputfn)
            cat("Pathway\t", file=pval.outfn)
            write.table(d, file=pval.outfn, sep="\t", quote=F, append=T)
            write.table(gsub(":.*","",ret$features), file=gsub("pvals","toppvals",pval.outfn), quote=F, row.names=F, col.names=F)
        }
    }

    ok.to.draw=TRUE
    if(sum(is.infinite(otu0)) != 0 | sum(is.na(otu0)) != 0)
    {   
        print("Found infinites or NAs in your OTU table")
        ok.to.draw <- FALSE
    }
    if(ncol(otu0) < 2 | nrow(otu0) < 2)
    {
        print("Only one differential feature found. Cannot generate heatmap.")
        ok.to.draw <- FALSE
    }

    if(ok.to.draw)
    {
        distfun <- function(x) as.dist((1-cor(t(x)))/2)

        # standardize each feature (taxa) so it stands out more (important so that color scheme magnifies properly)
        # rescale for better color gradient - makes the mean = 0 and sd = 1
        otu0 <- apply(otu0, 2, function(x) (x-mean(x))/sd(x))	
            #otu0 <- apply(otu0, 2, function(x) (x-mean(x))/(max(abs(min(x)), abs(max(x)))))	# divide by abs of min or max (depending on which is bigger)

        #cluster the columns (samples) first
        col.labels <- cluster.columns(otu=otu0, distfun=distfun, new.treatments=new.treatments)

        #then cluster the rows (OTUs)
        ret <- cluster.rows(otu=as.matrix(t(otu0))[,col.labels], distfun=distfun)
        row.dendrogram <- ret$ddr 

        # draw line between groups in order to see separation better
        last.groups <- aggregate(1:length(col.labels), list(map0[col.labels,group.var]), max)
        last.group.ix <- as.numeric(as.character(last.groups[,2]))
        last.group.ix <- last.group.ix[-which.max(last.group.ix)]
    
        col.colors <- create.color.bars(map0[col.labels,,drop=F], color.var=names(color.list), color.list)
        # flip this
        otu.hm <- as.matrix(t(otu0))

       draw.trad.heatmap(otu.hm=otu.hm[,col.labels], rowv = row.dendrogram, color.list=color.list, col.colors=col.colors, 
                        col.sep.ix=last.group.ix, outputfile=outputfn, show.rownames=show.rownames)
    }
}

# otu.hm must be ordered by row and col and rescaled as desired already
# set row.labels to NA if hide row names
draw.trad.heatmap <- function(otu.hm, rowv, color.list, col.colors, col.sep.ix, outputfile, show.rownames=FALSE)
{
	hm.colors <- colorRampPalette(c("#00441B", "white","#40004B"))(27)

	lwid = c(1.9,9,4) #col widths
	lhei = c(1.1,.3+.2*length(color.list),7.1,1) #row heights
	lmat = rbind(c(3,0,0), c(0,1,0), c(0,2,0), c(5,4,0))# 1=color side bars?, 2=heatmap?, 3=col dendogram, 4=row dendogram, 5 = key
	margins=c(4,6)
    legend1.x=0
    legend1.y=.95

    width=30
    height=16
    	
    if(show.rownames) labRow <- rownames(otu.hm)
	else labRow <- NA

		
	pdf(file=outputfile, width=width, height=height)	    
	
	heatmap.3(otu.hm, Colv = NA, Rowv = rowv, dendrogram='none',
            col = hm.colors, labRow=labRow,
            ColSideColors = col.colors,
            rowAxisColors=NULL,
            trace = "none", density.info = "none", 
            key=TRUE, keysize=0.75, 
            lmat=lmat, lwid=lwid, lhei=lhei, margins=margins, #labCol=NA,
            cexCol=.7, scale="none",
            colsep=col.sep.ix, sepwidth=c(.1,.1), sepcolor="white") # indices for where to draw row separators and width of line

	# add the legend for the colored bars
    legend(x=legend1.x, y=legend1.y,legend=unname(unlist(lapply(color.list, names))), 
    fill=unlist(color.list), border=FALSE, bty="n", cex=.9, ncol=2, text.width=strwidth("1,000,000"))
	
	dev.off()
}







