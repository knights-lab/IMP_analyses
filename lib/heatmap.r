# adapted from mwas.heatmap used for marty blaser's diabetes mice experiment
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/heatmap.3.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/collapse-features.r")
require(vegan)
require(RColorBrewer)
require(gplots) 



# hclusts the columns of an otu0 table according to groupings of samples
cluster.samples<-function(otu0, distfun, new.treatments)
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

cluster.otus<-function(otu0, distfun)
{
    # hclusts the rows, so make sure to transpose
    otu0 <- as.matrix(t(otu0))
	# do the row-clustering (note that we need to explicitly reorder the dendrogram here)
	Rowv <- rowMeans(otu0)
	hcr <- hclust(distfun(otu0))
	ddr <- as.dendrogram(hcr)
	ddr <- reorder(ddr, Rowv)
	rowInd <- order.dendrogram(ddr)
#		row.labels <- labels(ddr)[rowInd] #THIS IS CRITICAL?
	otu.labels <- labels(ddr)

	list(ddr=ddr, otu.labels=otu.labels)
}

# map0: rows must have already been reorder by the hclust
# color.list: list of named vectors that have been named according to color.var values
create.color.bars<-function(map0, color.var, color.list)
{
	col.colors <- NULL
	for(i in 1:length(color.var))
	{
		reordered <- map0[,color.var[i]]
		if(is.numeric(reordered)){
			color.palette <- colorRampPalette(color.list[[color.var[i]]])
			num.breaks <- 20 #? for gradients
			colors <- color.palette(num.breaks)[as.numeric(cut(reordered,breaks = num.breaks))]
		} else {
				colors <- unname(color.list[[color.var[i]]][as.character(reordered)])
		}
		col.colors <- cbind(col.colors, colors)		
	}
	colnames(col.colors) <- color.var
	t(col.colors)

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

make.heatmap <- function(taxa0, map0, min.prevalence=0.25, presence.absence=FALSE, baseline.groups=NULL, outputfn, order.by.num.otus=FALSE, override.rescale=NULL)
{
    valid.samples <- intersect(rownames(map0), rownames(taxa0))
    map0 <- map0[valid.samples,]
    taxa0 <- taxa0[valid.samples,]

    taxa0 <- sweep(taxa0, 1, rowSums(taxa0), '/')
    
     prevalences <- apply(taxa0, 2, function(bug.col) mean(bug.col > 0))
     taxa0 <- taxa0[,prevalences >= min.prevalence]
    
    # sort otus by prevalence within each group
    sorted.taxa <- NULL
    map0$Sample.Group <- factor(map0$Sample.Group) # remove levels not present
    groups <- levels(map0$Sample.Group)
    if(is.null(baseline.groups)){
        for(i in 1:length(groups))
        {
            this.taxa <- taxa0[map0$Sample.Group==groups[i],]
            prevalences <- apply(this.taxa, 2, function(bug.col) mean(bug.col > 0))
            # only add unique values to final ordered taxa list
            sorted.taxa <- union(sorted.taxa, names(sort(prevalences, decreasing=T)))
        }
    } else {
            this.taxa <- taxa0[map0$Sample.Group %in% baseline.groups,]
            prevalences <- apply(this.taxa, 2, function(bug.col) mean(bug.col > 0))
            # only add unique values to final ordered taxa list
            sorted.taxa <- union(sorted.taxa, names(sort(prevalences, decreasing=T)))
    }
    
    taxa0 <- taxa0[,sorted.taxa]
   
 
    # convert taxa table to 0 or 1
    rescale=TRUE
    if(presence.absence)
    {
        taxa0[taxa0 > 0] <- 1
        rescale <- FALSE
    }
    if(!is.null(override.rescale)) rescale <- override.rescale
    
    cluster.var <- c("Sample.Group")
    color.list <- list()
    # note that otu0 test differentiation will relevel factors based on this order (important for 2-group comparisons when looking at coefficients, 
    # and important for ordering the color groupings)
    # names must match columns exactly!

    color.list[["Sample.Group"]] <- get.group.colors(alpha.val=.8)
    
    color.list[["Years.in.US"]] <- setNames(c('#f3f6fb','#7c736c'), c("new arrival", "long-term"))

    color.list[["Sample.Group"]] <- color.list[["Sample.Group"]][names(color.list[["Sample.Group"]]) %in% unique(map0$Sample.Group)]
    
    core.mb.heatmap(otu0=taxa0, map0=map0, cluster.var=cluster.var, color.list=color.list, heatmap.title="", 
            outputfile=outputfn, rescale=rescale, order.by.num.otus=order.by.num.otus)
}


# this heatmap plots OTUs from left to right (highest to lowest prevalence) and
# samples from top to bottom (non-western to western)
# purpose: is to show depletion of microbes as we move down the rows
core.mb.heatmap <- function(otu0, map0, cluster.var, color.list, heatmap.title="", outputfile, rescale=TRUE, order.by.num.otus=FALSE)
{
    color.var <- names(color.list)
	new.treatments <- create.new.treatments(map0, cluster.var, color.list)

    # don't cluster OTUs either - should already be ordered by prevalence
    otu.labels <- colnames(otu0)
# 	ret <- cluster.otus(otu0=otu0[sample.labels,], distfun=distfun)
# 	otu.dendrogram <- ret$ddr 
# 	otu.labels <- ret$otu.labels
    
    # let's try clustering the rows for any samples that don't have a valid Years.in.US
    distfun <- function(x) as.dist((1-cor(t(x)))/2)
    partial.sample.labels <- cluster.samples(otu0=otu0[,otu.labels], distfun=distfun, new.treatments=map0$Sample.Group)
    partial.map <- map0[partial.sample.labels, c("Sample.Group","Years.in.US")] 
    
    # override sample labels for any 1st gen
    g <- levels(map0$Sample.Group)
    sample.labels <- NULL
    group.end.ix <- NULL # useful for marking separation between groups on the heatmap
    for(i in 1:length(g))
    {
        this.map <- partial.map[partial.map$Sample.Group==g[i],]        
    	
        if(order.by.num.otus)
        {
            this.otu <- otu0[rownames(this.map),]
            sample.labels <- c(sample.labels, rownames(this.otu[order(rowSums(this.otu), decreasing=TRUE),]))
        }
        else
        {
            if(g[i] %in% c("HmongThai","KarenThai","Hmong2nd"))
                sample.labels <- c(sample.labels, rownames(partial.map)[partial.map$Sample.Group==g[i]])
            else
                sample.labels <- c(sample.labels, rownames(this.map[order(this.map$Years.in.US),]))
        }
        group.end.ix <- c(group.end.ix, length(sample.labels))
    }  
	group.end.ix <- group.end.ix[-length(group.end.ix)] # don't need the last index
	
	# create each of the colored bars, reorder to previous column clustering
	col.colors <- create.color.bars(map0[sample.labels,,drop=F], color.var, color.list)
	

    lwid = c(1.2,.2,9,.5) #col widths
	lmat = rbind(c(0,0,4,0), c(3,1,2,0), c(0,0,5,0))
	        # 1 = color bars; 2 = heatmap
#	hm.colors <- colorRampPalette(c("white","plum"))(27)	# for presence absence only
	hm.colors <- colorRampPalette(c("#f2ece3","#764e80"))(2)	  # presence absence only
#
#	hm.colors <- colorRampPalette(c("blue","white","red"))(27)
#	hm.colors <- colorRampPalette(c("#2166ac","white","#b2182b"))(20)
	main <- heatmap.title
	
		lhei = c(1, 8, 1) #row heights
		margins=c(.05,.05)
		legend1.x=0
		legend1.y=.95
width= 30 #IMP
height=12 #IMP
	
	# standardize each feature (taxa) so it stands out more (important so that color scheme magnifies properly) (within samples)

	if(rescale)
	{
#	    otu0 <- apply(otu0, 2, function(x) (x-mean(x))/sd(x))	# makes the mean = 0 and sd = 1
	    otu0 <- apply(otu0, 2, function(x) (x-mean(x))/(max(abs(min(x)), abs(max(x)))))	# divide by abs of min or max (depending on which is bigger)
		otu0 <- sweep(otu0, 1, rowSums(otu0), '/')
        # try quantile norm for samples (instead) (appears that 1 sample is super high intensity -- messing up the others)
         
        
	    #otu0 <- apply(otu0, 2, function(x) (x-min(x))/(max(x)-min(x))) # 0 to 1
	    hm.colors <- colorRampPalette(c("blue","white","red"))(27)
	}
	pdf(file=outputfile, width=width, height=height)	
#dev.new(width=width, height=height) # for interactive plots only

    otu.for.heatmap <- as.matrix(otu0) 

    cnames <- rep("",length(otu.labels))

    cnames[colnames(otu.for.heatmap) == bacteroides] <- "B"
    cnames[colnames(otu.for.heatmap) == prevotella] <- "P"

	heatmap.3(otu.for.heatmap[sample.labels, otu.labels], Colv = NA, Rowv = NA, dendrogram='none',
	col = hm.colors, 
  RowSideColors = col.colors, # for now only
	rowAxisColors=NULL,
  trace = "none", density.info = "none", 
	main = main, key=FALSE, keysize=5, 
	lmat=lmat, lwid=lwid, lhei=lhei, margins=margins,
	labCol=cnames, labRow=NA, # remove all row and col labels
	cexCol=.7, scale="none",
	rowsep=group.end.ix, sepwidth=c(.1,.1), sepcolor="white") # indices for where to draw row separators and width of line

# super useful for clicking and identifying a location (for legend
# 		coords <- locator(1)
# 		print(coords)

	# add the legend for the colored bars
	legend(x=legend1.x, y=legend1.y,legend=unname(unlist(lapply(color.list, names))), 
	fill=unlist(color.list),border=FALSE, bty="n", cex=.9, ncol=2, text.width=strwidth("1,000,000"))#
																						#originally: ncol=3, text.width=c(0, .04, .06 )
																						#text.width= strwidth("1,000")
																						# PMP: ncol=2, text.width= strwidth("1,000,000")

	dev.off()
}







