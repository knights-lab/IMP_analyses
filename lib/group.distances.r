# plots bar graphs showing within and between differences, and permutation based tests for pvalues

# 1. within-group distances: Thai-Hmong, Thai-Karen, 1stHmong<10yrs, 1stHmong>10yrs 1stKaren, 2ndHmong
    # this shows what the spread is like in each group (very similar to each other or not?)
plot.within.group.distances <- function(map0, dm, fn, ylab)
{
    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    dm <- dm[valid_samples, valid_samples]

    ######## within-group variation ########
    # let's make a new grouping so that we separate <10 yr Hmong1st from the others
    map0[,"Sample.Group.2"] <- as.character(map0$Sample.Group)
    map0[map0$Sample.Group == "Hmong1st" & map0$Years.in.US <= 10,"Sample.Group.2"] <- "Hmong1st\n<= 10 yrs"
    map0[map0$Sample.Group == "Hmong1st" & map0$Years.in.US > 10,"Sample.Group.2"] <- "Hmong1st\n> 10 yrs"
    map0$Sample.Group.2 <- factor(map0$Sample.Group.2, levels=c("KarenThai","HmongThai","Karen1st","Hmong1st\n<= 10 yrs","Hmong1st\n> 10 yrs","Hmong2nd","Control")) 
    
    groups.within.dm <- NULL
    groups <- unique(map0$Sample.Group.2)
    for(i in 1:length(groups))
    {
        within.dist <- get.within.dist(dm, rownames(map0)[map0$Sample.Group.2==groups[i]])
        groups.within.dm <- rbind(groups.within.dm, 
                            data.frame(Distance=within.dist, Group=rep(groups[i], length(within.dist)), row.names=names(within.dist)))
    }
    groups.within.dm$Distance <- as.numeric(groups.within.dm$Distance)
    print(test.samples(y=groups.within.dm$Distance, x=groups.within.dm$Group))
    
    pdata <- data.frame(y = groups.within.dm$Distance, group = groups.within.dm$Group, row.names=rownames(groups.within.dm))
    p <- ggplot(pdata, aes(x = group, y = y)) + 
        geom_boxplot() + geom_jitter(width=.3, shape=21, colour=alpha("black",.5)) +
        ylab(ylab) + xlab("") + theme(axis.text.x = element_text(size=8)) +
        ggtitle("Within-Group Variability")
    save_plot(plot=p, fn, useDingbats=FALSE, base_aspect_ratio = 1.3 )

    invisible(list(plot=p, data=pdata))
}

plot.between.group.distances <- function(map0, dm, fn, ylab)
{
    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    dm <- dm[valid_samples, valid_samples]

    # between-group distances: 
    #   Thai-Hmong vs Thai-Karen
    #   1stHmong vs 1stKaren (0-10 yrs only)
    #       see if greater diffs between ethnicities after being in US for 10 yrs
    #   1stHmong<10 vs 2ndHmong
    #   Thai-Hmong vs 1stHmong<10
    #       see if greater diffs between thai-hmong and new immigrants, or 2ndgen and new immigrants
    #  
    between <- vector(mode="list")
    between[["HmongThai"]] <- rownames(map0)[map0$Sample.Group=="HmongThai"]
    between[["KarenThai"]] <- rownames(map0)[map0$Sample.Group=="KarenThai"]
    between[["Karen1st"]] <- rownames(map0)[map0$Sample.Group=="Karen1st"]
    between[["Hmong2nd"]] <- rownames(map0)[map0$Sample.Group=="Hmong2nd"]
    between[["Hmong1st<=10"]] <- rownames(map0)[map0$Sample.Group == "Hmong1st" & map0$Years.in.US <= 10]
    between[["Hmong1st>10"]] <- rownames(map0)[map0$Sample.Group == "Hmong1st" & map0$Years.in.US > 10]
    between[["Hmong1stObese"]] <- rownames(map0)[map0$Sample.Group == "Hmong1st" & map0$BMI.Class == "Obese"]
    between[["Hmong1stLean"]] <- rownames(map0)[map0$Sample.Group == "Hmong1st" & map0$BMI.Class == "Normal"]
    between[["Karen1stObese"]] <- rownames(map0)[map0$Sample.Group == "Karen1st" & map0$BMI.Class == "Obese"]
    between[["Karen1stLean"]] <- rownames(map0)[map0$Sample.Group == "Karen1st" & map0$BMI.Class == "Normal"]

        
    comparisons <- data.frame(group1=c("HmongThai","Karen1st","Hmong1st<=10", "Karen1stObese", "Hmong1stObese", "Hmong2nd"),
                            group2=c("KarenThai","KarenThai", "Karen1st", "Karen1stLean", "Hmong1stLean", "HmongThai"),
                            stringsAsFactors=F)
        
    groups.between.dm <- NULL
    for(i in 1:nrow(comparisons))
    {
        # make names for the groups for labeling later
        group1 <- between[[comparisons$group1[i]]]
        group2 <- between[[comparisons$group2[i]]]

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
        groups.between.dm <- rbind(groups.between.dm, data.frame(row.names=names(between.dist),
                                sampleid=names(between.dist), Distance=between.dist,
                                GroupLabel=rep(paste(comparisons$group1[i],comparisons$group2[i],sep="\n"), length(between.dist))))
    }    
    rownames(groups.between.dm) <- make.names(1:nrow(groups.between.dm))
    groups.between.dm$Distance <- as.numeric(groups.between.dm$Distance)
    
    # plot it
    pdata <- data.frame(y = groups.between.dm$Distance, group = groups.between.dm$GroupLabel, row.names=rownames(groups.between.dm)) 
    p <- ggplot(pdata, aes(x = group, y = y)) + 
        geom_boxplot() + geom_jitter(width=.3, shape=21, colour=alpha("black",.5)) +
        ylab(ylab) + xlab("") + theme(axis.text.x = element_text(size=8)) +
        ggtitle("Between-Group Variability")
    save_plot(plot=p, fn, useDingbats=FALSE, base_aspect_ratio = 1.3 )


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
    print(pvals)
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


		
		
		
		
		
		
		
		