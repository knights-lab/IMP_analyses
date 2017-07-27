
require(faraway)
# y = continuous var
# x.group - vector - how to group along the x axis
# x.subgroup - vector - how to group within groups along x axis (e.g. group by BMI class along decades in US)
boxplot.by.subgroup <- function(map0, x.group, x.subgroup=NULL, subgroup.cols, y, plot.title) #takes in a dataframe of mapping file + y = variable of interest ordered with map
{
    if(is.null(x.subgroup))
        d <- data.frame(Sub.Group=x.group, y=y, Group=x.group)
    else
        d <- data.frame(Sub.Group=x.subgroup, y=y, Group=x.group)
 
    p <- ggplot(d, aes(Group, y)) + geom_boxplot(aes(fill = Sub.Group)) +  geom_point(aes(y=y, fill = Sub.Group), position=position_dodge(width=.75), color="black", shape=21) +
        ggtitle(plot.title) + 
        theme(legend.position="none") +  
        ylab("log10 + 1") + xlab("") + theme(axis.text.x = element_text(size=9))

    if(!is.null(x.subgroup))
        p <- p + scale_fill_manual(name = "Groups", values = subgroup.cols)
        
    return(p)

}


# nutrients are already rownamed by samples
# map is a dataframe that has already been filtered for samples desired
plot.nutrients <- function(map0, ethnicity, nutrients0, nutrient.vars, independent.var="Years.in.US", control.var = "", bins=seq(0,45,5))
{
    valid_samples <- intersect(rownames(nutrients0), rownames(map0))
    map0 <- map0[valid_samples,]
    map0 <- map0[rownames(map0[map0$Ethnicity %in% ethnicity,]),]
    nutrients0 <- nutrients0[rownames(map0),nutrient.vars]
 
    # define first gen group here, to allow for different mapping files to be passed in
    # Sample ID named vector of Subject IDs
    firstgen_cs_samples <- rownames(map0)[map0$Sample.Group %in% c("Hmong1st","Karen1st") & (is.na(map0$Sample.Order) | map0$Sample.Order==1)]
    
    lookup <- alpha(c("green","orange","red"), .5)
    names(lookup) <- levels(map0$BMI.Class) 
    cols <- lookup[as.character(map0$BMI.Class)] 
    names(cols) <- rownames(map0)
    borders <- cols
    
    colbrewer <- alpha(brewer.pal(9,"RdPu"),.7)

    map0$Years.in.US[is.na(map0$Years.in.US)] <- -2
    map0$Years.in.US[map0$Years.in.US==0] <- 45
    cols[map0$Years.in.US==45] <- alpha(colbrewer[7], .3)
    borders[map0$Years.in.US==45] <- alpha(colbrewer[7], .3)
    cols[map0$Years.in.US==-2] <- colbrewer[1]
    borders[map0$Years.in.US==-2] <- colbrewer[7]
    
    # add 1 to all values in nutrients so we can log transform
    nutrients0 <- log10(nutrients0 + 1)
    
    bin.labels <- paste0("<",bins[-1])
    map0$Decade <- cut(map0$Years.in.US, bins, labels=bin.labels)
    map0$Decade <- as.character(map0$Decade)
    map0[map0$Years.in.US==-2,"Decade"] <- "Thai"

    if(sum(map0$Years.in.US==45)==0) # if no US born
        map0$Decade <- factor(map0$Decade, levels=c("Thai", bin.labels)) 
    else
    {
        map0[map0$Years.in.US==45,"Decade"] <- "USborn"
        map0$Decade <- factor(map0$Decade, levels=c("Thai", bin.labels, "USborn")) 
    }
    
    # do stats on firstgen only
    map00 <- map0[firstgen_cs_samples,]

    fn <- paste("Nutrient",paste0(ethnicity,collapse=""),independent.var,sep="_")
    cat("\n",file=paste(fn,".txt",sep=""), append=F)

    p <- NULL
    for(i in 1:length(nutrient.vars))
    {
        d <- data.frame(nutrient=nutrients0[,nutrient.vars[i]], x=map0[,independent.var])
                
         if(control.var == "BMI.Class")
        {       
            subgroup.cols = c("#ffffcc","#a1dab4","#41b6c4")
            p[[i]] <- boxplot.by.subgroup(map0, x.group=map0$Decade, x.subgroup=map0[,control.var], y=d$nutrient, plot.title=nutrient.vars[i], subgroup.cols = subgroup.cols)
        }
         else if(independent.var=="Sample.Group")
           p[[i]] <- boxplot.by.subgroup(map0, x.group=map0$Sample.Group, y=d$nutrient, plot.title=nutrient.vars[i])
           
    }

    multiplot <- plot_grid(plotlist=p, ncol=4, nrow=2)
    save_plot(paste(fn,".pdf",sep=""), multiplot,
          ncol = 4, # we're saving a grid plot of 2 columns
          nrow = 2, # and 2 rows
          # each individual subplot should have an aspect ratio of 1.3
          base_aspect_ratio = 1.3
          )
    
}

