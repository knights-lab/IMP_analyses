require(faraway)
# plot BMI and Waist-Height over years

plot.WHR.boxplot <- function(map0, bins=seq(0,10,2), fn)
{
    map0$Years.in.US[is.na(map0$Years.in.US)] <- -2
    map0$Years.in.US[map0$Years.in.US==0] <- 45

    bin.labels <- paste0("<",bins[-1])
    map0$Decade <- cut(map0$Years.in.US, bins, labels=bin.labels)
    map0$Decade <- as.character(map0$Decade)
    map0[map0$Years.in.US==-2,"Decade"] <- "Thai"

    if(sum(map0$Years.in.US==45)==0) # if no US born
        map0$Decade <- factor(map0$Decade, levels=c("Thai", bin.labels)) 
    else {
        map0[map0$Years.in.US==45,"Decade"] <- "USborn"
        map0$Decade <- factor(map0$Decade, levels=c("Thai", bin.labels, "USborn")) 
    }

    pdata<-map0
    p <- ggplot(pdata, aes(x = Decade, y = Waist.Height.Ratio)) + 
        geom_boxplot() + geom_jitter(width=.1) +
        ylab("") + xlab("Years in the US") +
        ggtitle("Waist-to-Height-Ratio")
    
    save_plot(plot=p, fn, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}

# bins the years then plots BMI class proportions over time
# freq = whether to plot y as absolute counts (TRUE) or as ratios (FALSE)
plot.BMI.barplot <- function(map0, bins=seq(0,10,2), freq=T, main="BMI Trends")
{
    gen1 <- map0[map0$Group == "1st",]
    bin.labels <- bins[-1]
    gen1$Decade <- cut(gen1$Years.in.US, bins, labels=bin.labels)
    gen1$Decade <- as.character(gen1$Decade)
    
    map0$Group.New <- as.character(map0$Group)
    map0[map0$Group == "1st","Group.New"] <- gen1$Decade
    map0$Group.New <- factor(map0$Group.New, levels=c("Pre", bin.labels, "2nd", "Control"))
    map0$Group.New <- factor(map0$Group.New) # remove levels not present
    
    # reverse the BMI order for stacking purposes
   map0$BMI.Class <- factor(map0$BMI.Class, levels = rev(levels(map0$BMI.Class)))
    #map0$BMI.Class <- factor(map0$BMI.Class)
    
    if(freq==TRUE) {position <- "stack" } else { position <- "fill" }

    colors <- get.bmi.colors()

    available.decades <- as.character(sort(as.numeric(unique(gen1$Decade))))
    x.breaks <- c(1, 1.5, 1:length(available.decades) + 1.5)
    this.levels <- levels(factor(map0$Group))
    extra.levels <- length(this.levels)-2
    while(extra.levels != 0) # if extra groups beyond firstgen are here, add them but not at the .5 marker
    {        
        x.breaks <- c(x.breaks, ceiling(x.breaks[length(x.breaks)]))
        extra.levels <- extra.levels - 1
    }

    x.labels <- c(this.levels[1], 0, available.decades, this.levels[-(1:2)])

#  print(cbind(x.labels,x.breaks))  
print(table(map0[,c("Group.New","BMI.Class")]))

    p <- ggplot(map0, aes(x = as.numeric(Group.New), fill = BMI.Class)) + 
        geom_bar(position=position) + scale_fill_manual(name="", values=colors) +
        ylab("") + xlab("") + facet_grid(. ~ Group, scales="free", space="free") +
        ggtitle(main) + theme(strip.background = element_blank(),strip.text.x = element_blank()) +
        scale_x_continuous(breaks = x.breaks, labels=x.labels) 

    return(p)
 }