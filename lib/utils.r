get.pch <- function(map)
{    
    lookup <- c(19,17) # point type
    names(lookup) <- sort(unique(map$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    pch <- lookup[as.character(map$Ethnicity)] 
    return(pch)
}


"shorten.taxonomy" <- function(ids,delim=';'){
	ids <- gsub('[kpcofgs]__','',ids)
	newids <- ids
	ids <- strsplit(ids,delim)
	for(i in seq_along(ids)){
		n <- length(ids[[i]])
		j <- n
		while(ids[[i]][j] == 'NA' || ids[[i]][j] == '') j <- j - 1 # NA instead of "Other" here
		newids[i] <- ids[[i]][j]
	}
	return(newids)
}

# x can be continuous or factor (vector)
# allows for control variables to be included (passed in as a data.frame)
# fits a linear model and assumes data is parametric (or has been transformed to be)
test.features.parametric <- function(otu, x, controls=NULL, sig.level, paired=FALSE)
{
    base.df <- x
    f <- "feature ~ x"
    if(!is.null(controls))
    {
        f <- paste0(f, " + ", paste(colnames(controls), collapse=" + "))
        base.df <- data.frame(x, controls)
    }
    ff <- as.formula(f)

    pvals <- apply(otu, 2, function(feature) {
            c <- summary(lm(ff, data.frame(feature, base.df)))$coefficients; return(c[2,4]);})
	
	adj.pvals <- p.adjust(pvals, "fdr")
	
	diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

	list(features=diff.features, pvals=adj.pvals)
}

# allows only for group comparisons, therefore x must be a factor
test.features.nonparametric <- function(otu, x, sig.level)
{
    stopifnot(is.factor(x))
    
    pvals <- apply(otu, 2, function(feature) 
        (kruskal.test(feature~x, data.frame(feature=feature, x=x)))$p.value)
	
	adj.pvals <- p.adjust(pvals, "fdr")
	
	diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

	list(features=diff.features, pvals=adj.pvals)
}

# allows only for two-group comparisons, therefore x must be a factor
# samples1 and samples2 are paired samples ordered by subject
test.features.paired <- function(otu, samples1, samples2, sig.level, parametric=FALSE)
{
    if(parametric) pvals <- apply(otu, 2, function(feature) t.test(x=feature[samples1], y=feature[samples2], paired=T)$p.value)
    else pvals <- apply(otu, 2, function(feature) wilcox.test(x=feature[samples1], y=feature[samples2], paired=T, exact=F)$p.value)
	
	adj.pvals <- p.adjust(pvals, "fdr")
	
	diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

	list(features=diff.features, pvals=adj.pvals)
}


# x = +2-level factor, y = continuous value
# simply runs wilcoxon test, OR manually loops through all possible combinations and corrects p-values for 3+ level group factors
# better than kruskal so we can see direction
test.groups <- function(y, x, parametric=TRUE)
{
    comparisons <- combn(levels(x), 2)
    
    pvals <- NULL
    pnames <- NULL
    for(i in 1:ncol(comparisons))
    {
        this.y <- y[x %in% comparisons[,i]]
        this.x <- x[x %in% comparisons[,i]]

        if(parametric) m <- t.test(this.y ~ this.x)
        else m <- wilcox.test(this.y ~ this.x)
        pvals <- c(pvals, m$p.value)
        pnames <- c(pnames, paste0(comparisons[1,i], "-", comparisons[2,i]))
    }
    names(pvals) <- pnames
    return(pvals)
}

# depending on what's available (map, otutable, and distance matrix), generates appropriate DM and PCs
# using samples that are currently in map0
# add.samples.dm: allows for additional samples to be included in dm that you don't care to keep in the mapping
prep.dm <- function(map0, otu0=NULL, dm=NULL, method="euclidean", axis1=1, axis2=2, add.samples.dm=NULL)
{
    if(!is.null(dm))
    {    
        valid_samples <- intersect(rownames(dm), rownames(map0))
        map0 <- map0[valid_samples,]
        if(!is.null(add.samples.dm)) 
            valid_samples <- c(valid_samples, add.samples.dm)
        dm <- dm[valid_samples, valid_samples]
        ddm <- as.dist(dm)
    }
    else
    {
        valid_samples <- intersect(rownames(otu0), rownames(map0))
        map0 <- map0[valid_samples,]
        if(!is.null(add.samples.dm)) 
            valid_samples <- c(valid_samples, add.samples.dm)
        otu0 <- otu0[valid_samples,]        
        ddm <- vegdist(otu0, method=method)
        dm <- as.matrix(ddm) # for use in stats later
    }
    pc <- cmdscale(ddm,max(axis1,axis2))
    return(list(map0=map0, pc=pc, dm=dm))
}


