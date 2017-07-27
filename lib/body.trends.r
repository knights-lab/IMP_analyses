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
plot.BMI.barplot <- function(map0, bins=seq(0,10,2), fn, freq=T)
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

    if(freq)
    {
        pdata<-map0
        ylab<-"Number of people"
        p <- ggplot(pdata, aes(x = Decade, fill = BMI.Class)) + 
            geom_bar() 
    }
    else
    {
        pdata <- table(map0$Decade,map0$BMI.Class)    
        pdata <- pdata/rowSums(pdata)
        pdata <- melt(pdata)
        colnames(pdata) <- c("Decade","BMI.Class","Ratio")
        ylab<-"Proportion of people"
        p <- ggplot(pdata, aes(x = Decade, y = Ratio, fill = BMI.Class)) + 
            geom_bar(stat="summary", fun.y="identity") 
    }
    
    colors <- c("#ffffcc","#a1dab4","#41b6c4")
    #colors <- alpha(c("green","yellow","red"),.7) # old colors
    
    p <- p +
        scale_fill_manual(labels = c("Lean","Overweight","Obese"), values=colors) + # legend labels and colors
        ylab(ylab) + xlab("Years in the US") +
        ggtitle("BMI Trends Over Time in the US")
    
    save_plot(plot=p, fn, useDingbats=FALSE, base_aspect_ratio = 1.3 )
 }