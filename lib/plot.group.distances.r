# plots bar graphs showing within and between differences, and permutation based tests for pvalues

plot.group.distances() <- function(map, dm, between.groups.fn)
{
    # peek at the data
    # table(map[,c("Ethnicity","BMI.Class","Residence.Class")])
    
    ######## within-group variation ########
    
    # generate all combinations for within group distances
    combdm <- expand.grid(Ethnicity=unique(map$Ethnicity), BMI.Class = unique(map$BMI.Class), Residence.Class=unique(map$Residence.Class))
    combdm <- combdm[order(combdm$Ethnicity, combdm$Residence.Class, combdm$BMI.Class),]
    rownames(combdm) <- 1:nrow(combdm)
    # grab samples for each group
    groups <- apply(combdm, 1, function(xx) 
                    {rownames(map[map$Ethnicity==xx[1] & map$BMI.Class==xx[2] & map$Residence.Class==xx[3],])}) 
    
    # calc within-group distances for all groups
    groups.within.dm <- NULL
    groups.bmi <- NULL
    for(i in 1:length(groups))
    {
        if(length(groups[[i]])>2) # remove any groups with small samples
        {        
            within.dist <- get.within.dist(dm, groups[[i]])
            groups.within.dm <- rbind(groups.within.dm, data.frame(Distance=within.dist,combdm[i,], Group=rep(paste(unlist(combdm[i,]),collapse="\n"), length(within.dist))))
            groups.bmi <- c(groups.bmi,as.character(combdm[i,"BMI.Class"]))
        }
    }
    groups.within.dm$Distance <- as.numeric(groups.within.dm$Distance)
    # plot it
    lookup <- c("green","orange","red")
    names(lookup) <- as.character(levels(map$BMI.Class))
    cols <- lookup[groups.bmi]
    par(cex.axis=.4)
    beeswarm(groups.within.dm$Distance ~ groups.within.dm$Group, col=cols,xlab="Groups",ylab="Unweighted Unifrac Distance")
    boxplot(groups.within.dm$Distance ~ groups.within.dm$Group,names=NA,ylab="",add=T)

    ######## between-group variation ########

    # read in table with specific between-group comparisons
    bgroups <- read.table(between.groups.fn,header=T,sep="\t", quote="", as.is=T)
    # set factor levels for variables appropriately so can be compared to mapping metadata
    bgroups$Residence.Class <- factor(bgroups$Residence.Class, levels=levels(map$Residence.Class))
    bgroups$BMI.Class <- factor(bgroups$BMI.Class, levels=levels(map$BMI.Class))
    bgroups$Ethnicity <- factor(bgroups$Ethnicity, levels=levels(map$Ethnicity))
        
    # calc between-group distances for each group comparison
    groups.between.dm <- NULL
    for(i in 1:length(unique(bgroups$Comparison)))
    {
        # make names for the groups for labeling later
        xx <- bgroups[bgroups$Comparison==i,]
        col.to.compare <- colnames(xx)[which(xx[1,]!=xx[2,])]
        cols <- unlist(lapply(xx[1,1:3],as.character)) # required to prevent factors to convert automatically to levels
        cols[col.to.compare] <- paste(xx[,col.to.compare],collapse="-v-")
        label <- paste(cols,collapse="\n")
        bgroups[bgroups$Comparison==i,"GroupLabel"] <- label

        group1 <- rownames(map[map$Ethnicity==xx[1,1] & map$BMI.Class==xx[1,2] & map$Residence.Class==xx[1,3],])
        group2 <- rownames(map[map$Ethnicity==xx[2,1] & map$BMI.Class==xx[2,2] & map$Residence.Class==xx[2,3],])

        # always put the larger group as the ref        
        if(length(group1) >= length(group2))
        {
            query <- group1
            ref <- group2
        }
        else
        {
            query <- group2
            ref <- group1
        }
        
        between.dist <- get.between.dist(dm,samples=query,ref_samples=ref)   
        # store this in a data.frame for easier plotting later
        groups.between.dm <- rbind(groups.between.dm, data.frame(sampleid=names(between.dist), 
                                Distance=between.dist,GroupLabel=rep(label, length(between.dist))))
    }    
    rownames(groups.between.dm) <- make.names(1:nrow(groups.between.dm))
    groups.between.dm$Distance <- as.numeric(groups.between.dm$Distance)
    # plot it
    par(cex.axis=.5)
    beeswarm(groups.between.dm$Distance ~ groups.between.dm$GroupLabel, col="gray",xlab=NA,ylab="Unweighted Unifrac Distance")
    boxplot(groups.between.dm$Distance ~ groups.between.dm$GroupLabel,names=NA,ylab="",add=T)

    # permutation test to calc differences between groups
    grouplabels <- unique(groups.between.dm$GroupLabel)
    
    # test groups 1,2 and 2,3 and 4,5 and 7,8
    tests <- rbind(c(1,2), c(2,3), c(4,5), c(7,8), c(5,6), c(6,9), c(7,9))
    pvals <- vector("numeric",length=nrow(tests))
    for(i in 1:nrow(tests))
    {
        group1 <- groups.between.dm[groups.between.dm$GroupLabel==grouplabels[tests[i,1]], "Distance"]
        names(group1) <- rownames(groups.between.dm[groups.between.dm$GroupLabel==grouplabels[tests[i,1]],])
        group2 <- groups.between.dm[groups.between.dm$GroupLabel==grouplabels[tests[i,2]], "Distance"]
        names(group2) <- rownames(groups.between.dm[groups.between.dm$GroupLabel==grouplabels[tests[i,2]],])
        pvals[i] <- permute.t.test(group1,group2)
    }
}


# all elements in groups must be named by sampleid
permute.t.test <- function(group1, group2, n=9999)
{
    obs.diff <- mean(group1) - mean(group2)
    
    all.data <- c(group1, group2)  # combine into one vector
    
    diffs <- vector("numeric", n)  # vector to store means
    for(i in 1:n) 
    {
        random.data <- sample(all.data)
        names(random.data) <- names(all.data)
        
        diffs[i] <- mean(random.data[names(group1)]) - mean(random.data[names(group2)])
    }
    diffs <- c(diffs, obs.diff)
        
 #    hist(diffs)
#     abline(v = obs.diff, lty = 2)
#     abline(v = -obs.diff, lty = 2)
    
    signif((sum(diffs <= -abs(obs.diff)) + sum(diffs >= abs(obs.diff)))/n)    
}

get.within.dist<-function(dm, samples)
{
    within <- numeric(length(samples))
    for(i in 1:length(samples))
    {
        within[i] <- mean(as.numeric(dm[samples[i], samples[-i]]))
    }
    names(within) <- samples
    return(within)
}

get.between.dist<-function(dm, samples, ref_samples)
{
    between <- numeric(length(samples))
    for(i in 1:length(samples))
    {
        between[i] <- mean(as.numeric(dm[samples[i], ref_samples]))
    }
    names(between) <- samples
    return(between)
}


		
		
		
		
		
		
		
		