# adapted from mwas.heatmap used for marty blaser's diabetes mice experiment
source(paste0(LIBDIR,"heatmap.3.r"))
source(paste0(LIBDIR,"collapse-features.r"))
source(paste0(LIBDIR,"utils.r"))
library(RColorBrewer)
require(vegan)
require(RColorBrewer)
require(gplots) 

# hclusts each row by treatment group (samples)
cluster.rows.by.group<-function(otu0, distfun, new.treatments)
{
	row.labels <- NULL
	# separate samples by unique values of new.treatments that we're interested in clustering by
	unique.new.treatments <- levels(new.treatments)

	for(i in 1:length(unique.new.treatments)){
		otu0.by.clustergroup <- otu0[new.treatments==unique.new.treatments[i], ]
 		data.dist.group <- distfun(otu0.by.clustergroup)
 		clus.group <- hclust(data.dist.group)

		d <- as.dendrogram(clus.group)
		# since the index order is per group only, we'll have to add the current length each time 
		# to get correct indices along the combined dendrogram
		row.labels <- c(row.labels, labels(d))
	}
	row.labels
	
}

# use for clustering OTUs
cluster.rows<-function(otu0.t, distfun)
{
	# do the row-clustering (note that we need to explicitly reorder the dendrogram here)
	Rowv <- rowMeans(otu0.t)
	hcr <- hclust(distfun(otu0.t))
	ddr <- as.dendrogram(hcr)
	ddr <- reorder(ddr, Rowv) # reorder dendrogram by OTU mean counts
	rowInd <- order.dendrogram(ddr)
    # row.labels <- labels(ddr)[rowInd] #not sure if this is needed
	otu.labels <- labels(ddr)

	list(ddr=ddr, row.labels=otu.labels)
}

# creates a data.frame of colors (samples as rows, columns as variables to color by)
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
	# this line makes the prevalence heatmap break: rownames(col.colors) <- rownames(map0)
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
    #alphas <- get.group.alphas()
    cols <- get.group.colors()
    cols["Hmong1st"] <- alpha(cols["Hmong1st"], .5)
    cols["Control"] <- alpha(cols["Control"], .7)
    
    color.list[["Sample.Group"]] <- cols
    color.list[["Sample.Group"]] <- color.list[["Sample.Group"]][names(color.list[["Sample.Group"]]) %in% unique(map0$Sample.Group)]

    color.list[["Age"]] <- setNames(c('white','#CC9900'), as.character(signif(range(map0$Age),2)))
    
    color.list[["Years.in.US"]] <- setNames(c('white','#005073'), as.character(signif(range(map0$Years.in.US),2)))


    return(color.list)
}

# for use with paired test on Longitudinal samples only
# get.color.list<-function(map0)
# {
#     color.list <- list()
#     
#     color.list[["Sample.Month"]] <- c("1"="black","6"="red")
# 
#     return(color.list)
# }


# This heatmap converts OTUs to presence/absence, and plots OTUs from left to right (highest to lowest prevalence) and samples from top to bottom (non-western to western)
#   sorts OTUs by prevalence in baseline groups (min prevalence cut off also only applied at baseline)
#   sorts samples by richness

make.heatmap.binary <- function(otu0, map0, group.var="Sample.Group", min.prevalence=0.25, baseline.groups, show.colnames=FALSE, sig.level=NULL, outputfn, is.otu=FALSE, taxamapfn=NULL)
{
    valid.samples <- intersect(rownames(map0), rownames(otu0))
    map0 <- map0[valid.samples,]
    otu0 <- otu0[valid.samples,]

    map0[,group.var] <- factor(map0[,group.var]) # remove levels not present after filtering
        
    prep.ret <- prep.otu.binary.heatmap(otu0=otu0, map0=map0, group.var=group.var, baseline.groups=baseline.groups, min.prevalence=min.prevalence, sig.level=sig.level)
    otu0 <- prep.ret$otu
    output.df <- prep.ret$output.df
        
    # do some special wrangling here so we can see what taxa the OTUs map to
    if(is.otu)
    {
        tax <- read.table(file=taxamapfn, header=F, sep="\t", colClasses="character")
        colnames(tax) <- c("otuid","taxa")
        output.df$feature <- as.character(output.df$feature)
        output.df <- merge(output.df, tax, by.x="feature", by.y="otuid")
        # reorder to match original row
        rownames(output.df) <- output.df$feature
        output.df <- output.df[prep.ret$output.df$feature,]
        output.df$short.taxa <- shorten.taxonomy(output.df$taxa, min.char=3)
    }
    write.table(output.df, file=gsub("pdf","txt",outputfn), row.names=F, quote=F, sep="\t")

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

    if(sum(is.infinite(otu0)) != 0 | sum(is.na(otu0)) != 0 | ncol(otu0) < 2 | nrow(otu0) < 2)
        print("Your OTU table is too small or contains NAs. Not generating heatmap.")
    else
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

# loops through taxa, and does a fisher's test using a prevalence otu, then corrects by fdr
test.features.fishers <- function(otu0, groups, sig.level)
{
    pval <- NULL
    for(i in 1:ncol(otu0))
    {
        count.table <- table(otu0[,i], groups)
        
        # if no zero counts at all, add it manually (note, no need to add 1s manually due to prevalence filter)
        if(nrow(count.table) != 2) count.table <- rbind(count.table, "0"=c(0,0))
            
        pval[i] <- fisher.test(count.table)$p.value
    }
    qval <- p.adjust(pval, method="fdr")


    sig.features <- qval[qval < sig.level]
    names(sig.features) <- colnames(otu0)[qval < sig.level]
    return(sig.features)
}

prep.otu.binary.heatmap <- function(otu0, map0, group.var, baseline.groups, min.prevalence=0.25, sig.level=NULL)
{    
    # convert table to presence absence
    otu0[otu0 > 0] <- 1

    # filter min prevalence based on the baseline group only
    if(!is.null(min.prevalence)){ 
        # filter prevalence within the baseline sample group
        this.taxa <- otu0[map0[,group.var] %in% baseline.groups,]
        prevalences <- apply(this.taxa, 2, function(bug.col) mean(bug.col > 0))
        prevalences <- prevalences[prevalences >= min.prevalence]
        otu0 <- otu0[,names(prevalences)]
    }
    output.df <- data.frame(prevalences, row.names=names(prevalences))
    
    if(!is.null(sig.level)){
        # select features that are statistically different between groups (fisher's exact)
        sig.features <- test.features.fishers(otu0, map0[,group.var], sig.level=sig.level)
        otu0 <- otu0[,names(sig.features),drop=F]	

        # regenerate new prevalences with non-sig OTUs removed
        this.taxa <- otu0[map0[,group.var] %in% baseline.groups,]
        prevalences <- apply(this.taxa, 2, function(bug.col) mean(bug.col > 0))
        
        this.taxa.other <- otu0[!(map0[,group.var] %in% baseline.groups),]
        prevalences.other <- apply(this.taxa.other, 2, function(bug.col) mean(bug.col > 0))
    }
    # now sort based on prevalences of OTUs that remain, within baseline group only
    sorted.taxa <- names(sort(prevalences, decreasing=T))
    otu0 <- otu0[,sorted.taxa]
    output.df <- data.frame(feature=sorted.taxa, qval=sig.features[sorted.taxa], prevalence.baseline=prevalences[sorted.taxa], prevalence.nonbaseline=prevalences.other[sorted.taxa], stringsAsFactors=F)

    # sort rows by richness (number of otus in sample)
    final.otu <- NULL
    for(this.group in levels(map0[,group.var]))
    {
        this.otu <- otu0[map0[,group.var]==this.group,,drop=F]

        if(nrow(this.otu) > 1) # only reorder if we have more than 1 sample!
            this.otu <- this.otu[order(rowSums(this.otu)),]
            
        final.otu <- rbind(final.otu, this.otu)
    }
    
    return(list(otu=final.otu, output.df=output.df))
}

# keeps prevalent OTUs, collapse correlated OTUs, transforms rel abundances if parametric testing, keeps only significant OTUs
<<<<<<< HEAD
prep.otu.trad.heatmap<-function(otu0, groups, min.prevalence=0.25, sig.level=.25, test="parametric", normalize=TRUE, paired=FALSE)
=======
prep.otu.trad.heatmap<-function(otu0, groups, min.prevalence=0.25, sig.level=.25, test="parametric", normalize=TRUE)
>>>>>>> d870284bddbbb62205128d8e65887970bb65e795
{
    # if not normalizing (e.g. clr), do not convert to relative abundance, filter low prevalence, or sqrt transform
    if(normalize){ 
        otu0 <- sweep(otu0, 1, rowSums(otu0), '/')	
        prevalences <- apply(otu0, 2, function(bug.col) mean(bug.col > 0))
        otu0 <- otu0[, prevalences >= min.prevalence]
    }
    ret <- collapse.by.correlation(otu0, .95)
    otu0 <- otu0[, ret$reps]

    # select features that are statistically different between groups (parametric)
    if(test=="parametric") {
<<<<<<< HEAD
        otu0 <- asin(sqrt(otu0))
        ret <- test.features.parametric(otu0, groups, sig.level=sig.level, paired=paired)
        otu0 <- otu0[,ret$features,drop=F]	
    } else if(test=="nonparametric") {
        ret <- test.features.nonparametric(otu0, groups, sig.level=sig.level, paired=paired)
=======
        ret <- test.features.parametric(otu0, groups, sig.level=sig.level)
        otu0 <- otu0[,ret$features,drop=F]	
    } else if(test=="nonparametric") {
        ret <- test.features.nonparametric(otu0, groups, sig.level=sig.level)
>>>>>>> d870284bddbbb62205128d8e65887970bb65e795
	    otu0 <- otu0[,ret$features,drop=F]	
	}

    return(list(otu=otu0, ret=ret))
}


make.bubble.heatmap <- function(otu0, map0, outputfn, group.var="Sample.Group", filter.mode="standard", save.pvals=FALSE, min.prev=0.1, sig.level=.10, 
                            show.rownames=FALSE, color.rownames=FALSE, show.samplenames=FALSE, rescale="standard", as.bubble=FALSE, cov.df=NULL, baseheight=4, basewidth=15)
{
    make.heatmap.traditional.base(otu0=otu0, map0=map0, outputfn=outputfn, group.var=group.var, filter.mode=filter.mode, save.pvals=save.pvals, min.prev=min.prev, sig.level=sig.level, 
                            show.rownames=show.rownames, color.rownames=color.rownames, show.samplenames=show.samplenames, rescale=rescale, as.bubble=as.bubble, cov.df=cov.df, baseheight=baseheight, basewidth=basewidth)

}

make.heatmap.traditional <- function(otu0, map0, outputfn, group.var="Sample.Group", filter.mode="standard", save.pvals=FALSE, min.prev=0.1, sig.level=.10, 
                            show.rownames=FALSE, color.rownames=FALSE, show.samplenames=FALSE, rescale="standard", paired=FALSE)
{
    make.heatmap.traditional.base(otu0=otu0, map0=map0, outputfn=outputfn, group.var=group.var, filter.mode=filter.mode, save.pvals=save.pvals, min.prev=min.prev, sig.level=sig.level, 
                            show.rownames=show.rownames, color.rownames=color.rownames, show.samplenames=show.samplenames, rescale=rescale, as.bubble=FALSE, cov.df=NULL, paired=FALSE)
}

# This makes a traditional heatmap (samples as columns, non-western to western and OTUs as rows)
# min.prevalence applied globally 
<<<<<<< HEAD
make.heatmap.traditional.base <- function(otu0, map0, outputfn, group.var="Sample.Group", filter.mode="standard", save.pvals=FALSE, min.prev=0.1, sig.level=.10, 
                            show.rownames=FALSE, color.rownames=FALSE, show.samplenames=FALSE, rescale="standard", as.bubble=FALSE, cov.df=NULL, baseheight=4, basewidth=15, paired=FALSE)
=======
make.heatmap.traditional <- function(otu0, map0, outputfn, group.var="Sample.Group", filter.mode="standard", save.pvals=FALSE, min.prev=0.1, sig.level=.10, show.rownames=FALSE)
>>>>>>> d870284bddbbb62205128d8e65887970bb65e795
{
    valid.samples <- intersect(rownames(map0), rownames(otu0))
    map0 <- map0[valid.samples,]
    otu0 <- otu0[valid.samples,]

    map0[,group.var] <- factor(map0[,group.var]) # remove levels not present

    # useful to see distribution of samples
    print(table(map0[,group.var]))

    color.list <- get.color.list(map0)
<<<<<<< HEAD

    new.treatments <- create.new.treatments(map0, cluster.var=group.var, color.list)

print(new.treatments)
=======
>>>>>>> d870284bddbbb62205128d8e65887970bb65e795

    new.treatments <- create.new.treatments(map0, cluster.var=group.var, color.list)
    
    # transforms and filters otu using standard methods
    if(filter.mode!="none")
    {
        if(filter.mode=="standard") 
            prep.out <- prep.otu.trad.heatmap(otu0, groups=new.treatments, min.prevalence=min.prev, sig.level=sig.level, paired=paired)
        else if (filter.mode=="clr")
            prep.out <- prep.otu.trad.heatmap(otu0, groups=new.treatments, min.prevalence=min.prev, sig.level=sig.level, normalize=FALSE, paired=paired)
        else if (filter.mode=="nonparametric")
            prep.out <- prep.otu.trad.heatmap(otu0, groups=new.treatments, min.prevalence=min.prev, sig.level=sig.level, test="nonparametric", normalize=TRUE, paired=paired)

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

<<<<<<< HEAD
    if(sum(is.infinite(as.matrix(otu0))) != 0 | sum(is.na(otu0)) != 0 | ncol(otu0) < 2 | nrow(otu0) < 2)
=======
    if(sum(is.infinite(otu0)) != 0 | sum(is.na(otu0)) != 0 | ncol(otu0) < 2 | nrow(otu0) < 2)
>>>>>>> d870284bddbbb62205128d8e65887970bb65e795
        print("Your OTU table is too small or contains NAs. Not generating heatmap.")
    else
    {

        distfun <- function(x) as.dist((1-cor(t(x)))/2)

        # standardize each feature (taxa) so it stands out more (important so that color scheme magnifies properly)
        # rescale for better color gradient - makes the mean = 0 and sd = 1
        if(rescale=="standard")
            otu0 <- apply(otu0, 2, function(x) (x-min(x))/(max(x) - min(x)))	# this rescaling seems to work nicest
        else if (rescale=="zero-mean")
            otu0 <- apply(otu0, 2, function(x) (x-mean(x))/sd(x)) # zero-mean, unit variance (best displayed as 3 colors)	
            #otu0 <- apply(otu0, 2, function(x) (x-mean(x))/(max(x) - min(x)))	# mean normalization

        # cluster the samples within each group
        sample.labels <- cluster.rows.by.group(otu=otu0, distfun=distfun, new.treatments=new.treatments)

        # pass in transposed OTU so we can cluster the OTUs 
        ret <- cluster.rows(otu=as.matrix(t(otu0[sample.labels,])), distfun=distfun)
        otu.dendrogram <- ret$ddr 
        otu.labels <- ret$row.labels

        last.group.ix <- NULL
        # draw line between groups in order to see separation better
<<<<<<< HEAD
        last.groups <- aggregate(1:length(sample.labels), list(map0[sample.labels,group.var]), max)
        last.group.ix <- as.numeric(as.character(last.groups[,2]))
        last.group.ix <- last.group.ix[-which.max(last.group.ix)]
        col.colors <- create.color.bars(map0[sample.labels,,drop=F], color.var=names(color.list), color.list)

        # flip: in a traditional heatmap, let's have samples as columns and features as rows
        otu.hm <- as.matrix(t(otu0))


        if(as.bubble)
        {        
            dcols <- get.group.colors()
            dcols["Hmong1st"] <- alpha(dcols["Hmong1st"], .7)
            dcols["Control"] <- alpha(dcols["Control"], .7)            
            # make sure to reorder the otu table based on clustering above
            draw.bubble.heatmap(otu.hm=otu.hm[otu.labels, sample.labels], cols = dcols[as.character(map0[sample.labels,group.var])], 
                        outputfile=outputfn, cov.df=cov.df, show.samplenames=show.samplenames, baseheight=baseheight, basewidth=basewidth)
        
        }else {
            # dendrogram will automatically reorder the rows
            # make sure to pass in the new ordered sample labels
           draw.trad.heatmap(otu.hm=otu.hm[,sample.labels], rowv = otu.dendrogram, color.list=color.list, col.colors=col.colors, 
                        col.sep.ix=last.group.ix, outputfile=outputfn, show.rownames=show.rownames, color.rownames=color.rownames, show.samplenames=show.samplenames,rescale=rescale)
        }
    }
}

# This makes a traditional heatmap (samples as columns, non-western to western and OTUs as rows)
# plots everything as-is, no testing for diff features, no filtering, nuthin'!
make.heatmap.traditional.nogroups <- function(otu0, map0, outputfn, show.rownames=FALSE, sample.clustering=TRUE, rescale="standard")
{
    valid.samples <- intersect(rownames(map0), rownames(otu0))
    map0 <- map0[valid.samples,]
    otu0 <- otu0[valid.samples,]

    color.list <- get.color.list(map0)
    
    if(sum(is.infinite(otu0)) != 0 | sum(is.na(otu0)) != 0 | ncol(otu0) < 2 | nrow(otu0) < 2)
        print("Your OTU table is too small or contains NAs. Not generating heatmap.")
    else
    {
        distfun <- function(x) as.dist((1-cor(t(x)))/2)

        # standardize each feature (taxa) so it stands out more (important so that color scheme magnifies properly)
        # rescale for better color gradient - makes the mean = 0 and sd = 1        
        otu0 <- apply(otu0, 2, function(x) (x-mean(x))/sd(x))	
        
        #cluster the samples
        sample.labels <- rownames(otu0)

        if(sample.clustering)
        {
            ret <- cluster.rows(otu=otu0, distfun=distfun)
            sample.labels <- ret$row.labels
        }
        # then cluster the OTUs
        ret <- cluster.rows(otu=as.matrix(t(otu0[sample.labels,])), distfun=distfun)
        otu.dendrogram <- ret$ddr 
        otu.labels <- ret$row.labels
        
        last.group.ix <- NULL

        col.colors <- create.color.bars(map0[sample.labels,,drop=F], color.var=names(color.list), color.list)
=======
        last.groups <- aggregate(1:length(col.labels), list(map0[col.labels,group.var]), max)
        last.group.ix <- as.numeric(as.character(last.groups[,2]))
        last.group.ix <- last.group.ix[-which.max(last.group.ix)]
        col.colors <- create.color.bars(map0[col.labels,,drop=F], color.var=names(color.list), color.list)
>>>>>>> d870284bddbbb62205128d8e65887970bb65e795
        # flip this
        otu.hm <- as.matrix(t(otu0))

       draw.trad.heatmap(otu.hm=otu.hm[,sample.labels], rowv = otu.dendrogram, color.list=color.list, col.colors=col.colors, 
                        col.sep.ix=last.group.ix, outputfile=outputfn, show.rownames=show.rownames, rescale=rescale)
    }
}

# otu.hm must be ordered by row and col and rescaled as desired already
# rowname.colors = vector of colors mapped to rownames
draw.trad.heatmap <- function(otu.hm, rowv, color.list, col.colors, col.sep.ix, outputfile, show.rownames=FALSE, color.rownames=FALSE, show.samplenames=TRUE, rescale="standard")
{
    if(rescale %in% c("none","standard"))
        hm.colors <- colorRampPalette(c("white","#40004B"))(27) # if not rescaled, then only use two colors
	else        
	    hm.colors <- colorRampPalette(c("#00441B", "white","#40004B"))(27)

	
	lwid = c(1.9,9,4) #col widths
	lhei = c(1.1,.3+.2*length(color.list),7.1,1) #row heights
	lmat = rbind(c(3,0,0), c(0,1,0), c(0,2,0), c(5,4,0))# 1=color side bars?, 2=heatmap?, 3=col dendogram, 4=row dendogram, 5 = key
	margins=c(4,6)
    legend1.x=0
    legend1.y=.95

    width=30
    height=14 #16
    	
    if(show.rownames) labRow <- rownames(otu.hm)
	else labRow <- NA

    if(show.samplenames) labCol <- colnames(otu.hm)
	else labCol <- NA
		
    rowname.colors <- NULL
    if(color.rownames)
    {
        level = "L1"
        ont <- read.table("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/project_043_shotgun/humann2/metacyc/metacyc.ontology.txt", sep="\t", comment="", head=T, quote="", colClasses="character")
        otu.ont <- merge(as.data.frame(gsub(".*: ", "", rownames(otu.hm))), ont, by.x=1, by.y="Pathways", all.x=T)
        otu.ont[is.na(otu.ont[,level]), level] <- "No Ontology"
        text.cols <- c("#e41a1c","#f781bf","#4daf4a","#984ea3","#ff7f00","#377eb8") 
        sorted.levels <- sort(unique(otu.ont[,level]))
        sorted.levels <- sorted.levels[sorted.levels != "No Ontology"]
        names(text.cols) <- sorted.levels
        text.cols <- c(text.cols, "No Ontology"="#483D41")

        rowname.colors <- text.cols[as.character(otu.ont[,level])] 
    } 
        
	pdf(file=outputfile, width=width, height=height)	    
	heatmap.3(otu.hm, Colv = NA, Rowv = rowv, dendrogram='none',
            col = hm.colors, labRow=labRow, labCol=labCol,
            ColSideColors = col.colors,
            trace = "none", density.info = "none", 
            key=TRUE, keysize=0.75, 
            lmat=lmat, lwid=lwid, lhei=lhei, margins=margins, 
            cexCol=.7, scale="none", rowAxisColors=rowname.colors,
            colsep=col.sep.ix, sepwidth=c(.1,.1), sepcolor="white") # indices for where to draw row separators and width of line

	# add the legend for the colored bars
    legend(x=legend1.x, y=legend1.y,legend=unname(unlist(lapply(color.list, names))), 
    fill=unlist(color.list), border=FALSE, bty="n", cex=.9, ncol=2, text.width=strwidth("1,000,000"))
	
	if(color.rownames)
	{
		legend(x=.705, y=.08,legend=gsub("\\(.*","",names(text.cols)),fill=text.cols, border=FALSE, bty="n", y.intersp = 0.8, cex=0.7, text.col="dimgray")
    }
	dev.off()
}

<<<<<<< HEAD
# cov.df = same as OTU in samples and features, but instead shows coverage
# cols = colors with names as sample id
draw.bubble.heatmap <- function(otu.hm, cols=NULL, outputfile, cov.df, show.samplenames, baseheight=4, basewidth=15)
{
    # manually stick B and P together as a group
    otu.levels <- c(grep("^B", rownames(otu.hm), value=T), grep("^P", rownames(otu.hm), value=T))

    # manually order longitudinal samples by name
    sample.levels <- colnames(otu.hm)
    sample.levels[grep("^L", sample.levels, value=F)] <- sort(grep("^L", sample.levels, value=T))

    otu.long <- melt(otu.hm)
    colnames(otu.long) <- c("otu.id", "sample.id", "count")
    otu.long$otu.id <- as.character(otu.long$otu.id)
    otu.long$sample.id <- as.character(otu.long$sample.id)
=======
# This makes a traditional heatmap (samples as columns, non-western to western and OTUs as rows)
# plots everything as-is, no testing for diff features, no filtering, nuthin'!
make.heatmap.traditional.nogroups <- function(otu0, map0, outputfn, show.rownames=FALSE, min.prevalence=.1)
{
    valid.samples <- intersect(rownames(map0), rownames(otu0))
    map0 <- map0[valid.samples,]
    otu0 <- otu0[valid.samples,]

    otu0 <- sweep(otu0, 1, rowSums(otu0), '/')	
    prevalences <- apply(otu0, 2, function(bug.col) mean(bug.col > 0))
    otu0 <- otu0[, prevalences >= min.prevalence]

    color.list <- get.color.list(map0)
    
    if(sum(is.infinite(otu0)) != 0 | sum(is.na(otu0)) != 0 | ncol(otu0) < 2 | nrow(otu0) < 2)
        print("Your OTU table is too small or contains NAs. Not generating heatmap.")
    else
    {
        distfun <- function(x) as.dist((1-cor(t(x)))/2)

        # standardize each feature (taxa) so it stands out more (important so that color scheme magnifies properly)
        # rescale for better color gradient - makes the mean = 0 and sd = 1
        #apply(otu0, 2, function(x) print(sd(x)))
        
        otu0 <- apply(otu0, 2, function(x) (x-mean(x))/sd(x))	
        
        #cluster the samples
        ret <- cluster.rows(otu=otu0, distfun=distfun)
        col.labels <- ret$row.labels

        #then cluster the OTUs
        ret <- cluster.rows(otu=as.matrix(t(otu0))[,col.labels], distfun=distfun)
        row.dendrogram <- ret$ddr 

        last.group.ix <- NULL

        col.colors <- create.color.bars(map0[col.labels,,drop=F], color.var=names(color.list), color.list)
        # flip this
        otu.hm <- as.matrix(t(otu0))

       draw.trad.heatmap(otu.hm=otu.hm[,col.labels], rowv = row.dendrogram, color.list=color.list, col.colors=col.colors, 
                        col.sep.ix=last.group.ix, outputfile=outputfn, show.rownames=show.rownames)
    }
}

>>>>>>> d870284bddbbb62205128d8e65887970bb65e795

    # otu is otus as rows and samples as cols - flip covdf to match
    cov.long <- melt(t(cov.df[colnames(otu.hm), rownames(otu.hm)]))
    colnames(cov.long) <- c("otu.id", "sample.id", "coverage")
    cov.long$otu.id <- as.character(cov.long$otu.id)
    cov.long$sample.id <- as.character(cov.long$sample.id)

    ggdata <- merge(otu.long, cov.long, by=c("sample.id","otu.id"))

    xcols <- data.frame(sample.id=sample.levels, group = names(cols), cols=cols, stringsAsFactors=F)
    xcols$group <- factor(xcols$group, levels=SAMPLE.GROUP.NAMES[SAMPLE.GROUP.NAMES %in% unique(xcols$group)])
     freqdf <- as.data.frame(table(xcols$group))
     freq <- freqdf$Freq
     names(freq) <- freqdf$Var1
     sample.order <- unlist(lapply(freq, function(xx) 1:xx))
    xcols$sample.order <- sample.order

    ggdata <- merge(ggdata, xcols, by="sample.id") # add colors as groups, and a per group sample order

    # reset the levels as the established order of otu.hm so that it appears correctly in plot
    ggdata$otu.id <- factor(ggdata$otu.id, levels=otu.levels)
    ggdata$sample.id <- factor(ggdata$sample.id, levels=sample.levels)

#write.table(ggdata, file="ggdata.txt", sep="\t", quote=F)

    # anything that is covered at less than 1%, consider this to not be there
    ggdata$coverage[ggdata$coverage < 0.01] <- NA # setting to 0 still plots the point, need to set as NA

    plist <- list()
    groups <- levels(ggdata$group)
    for(group in groups)
    {
        this.ggdata <- ggdata[ggdata$group %in% group,]
        group.color <- this.ggdata$cols[1]
        p <- ggplot(this.ggdata, aes(x=sample.id, y=otu.id)) + 
            geom_point(aes(size = coverage, color=count)) +
            scale_color_gradient(low="#F0F5EF", high="#40004B") +  
            # scale_size(range = c(0,8)) + 
            scale_size_area(max_size = 8) + 
            scale_x_discrete(position = "bottom") + 
            facet_grid(. ~ group, scales="free", drop=TRUE) + 
            theme(  axis.title=element_blank(), axis.ticks.x = element_blank(),
                    strip.text = element_text(colour = 'white'),
                    strip.background=element_rect(fill=group.color),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_blank()) 

        if(show.samplenames)
            p <- p + theme(axis.text.x = element_text(angle=90, size=5))
        else 
            p <- p + theme(axis.text.x = element_blank())

        if(group != groups[length(groups)]) # if not last group, remove all borders, axes
        {
            p <- p + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line=element_blank()) 
        }
        else # if last group, add y-axis only
        {
            # custom y-axis with colors
            cols <- c("#539798", "#A73E5D")
            names(cols) <- c("P","B")
            ylabs <- substring(otu.levels,1,1) 
            ycols <- cols[ylabs]
            ylabs[ylabs=="B"] <- paste0("B",sum(ylabs=="B"):1)
            ylabs[ylabs=="P"] <- paste0("P",sum(ylabs=="P"):1)

            p <- p + scale_y_discrete(position = "right",labels=ylabs) + 
                    theme(axis.text.y = element_text(size=10,face="bold",colour=ycols), axis.line.x=element_blank()) 
            # save legend
            leg <- get_legend(p)
        }
        p <- p + theme(legend.position='none')
        plist[[group]] <- p
    }
    
    # save strains list with coverage values, reverse the labels to match the heatmap (goes bottom to up)
    write.table(cbind(Legend.ID=rev(ylabs), Genome.ID=rev(otu.levels), Coverage=colMeans(cov.df[,rev(otu.levels)])), file=gsub("\\.pdf", "\\.strains\\.txt", outputfile), row.names=F, col.names=T, quote=F, sep="\t")
#print(freqdf)
#print(names(plist))
    # save legend
    save_plot(gsub("\\.pdf", "\\.legend\\.pdf", outputfile), leg, base_aspect_ratio=1)
    pgrid.cs <- plot_grid(plotlist=plist[c("HmongThai","Hmong1st","Control")], nrow=1, rel_widths=freq[c("HmongThai","Hmong1st","Control")], align="h")
    pgrid.kt <- plist[["Karen1st"]]
    
    save_plot(gsub("\\.pdf", "\\.CS\\.pdf", outputfile), pgrid.cs, base_height=baseheight, base_width=basewidth)
    save_plot(gsub("\\.pdf", "\\.KT\\.pdf", outputfile), pgrid.kt, base_height=baseheight, base_width=basewidth/4)
}





