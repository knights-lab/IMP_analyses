library(ggplot2)
library(RColorBrewer)
require(ggbeeswarm)
library(cowplot)

setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

mice_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/mice.txt"
chow_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/chow.txt"

mice <- read.table(file=mice_fn, sep="\t", header=T, as.is=T)
chow <- read.table(file=chow_fn, sep="\t", header=T, as.is=T)

# rename colnames
colnames(mice)[5] <- "Fasting.Glucose"
colnames(mice)[6] <- "Weight"

    mouse.ids <- unique(mice$Mouse.ID)

    baseline.weights <- mice[mice$Round==1,]
    ggdata <- mice[mice$Round %in% 2:5,]        

    all_temp <- NULL
    for(i in 1:length(mouse.ids))
    {
        temp <- ggdata[ggdata$Mouse.ID == mouse.ids[i],"Weight"] - baseline.weights[baseline.weights$Mouse.ID == mouse.ids[i], "Weight"]
        names(temp) <- paste( ggdata[ggdata$Mouse.ID == mouse.ids[i], "Mouse.ID"], ggdata[ggdata$Mouse.ID == mouse.ids[i], "Round"], sep=".")
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Weight.Gain"] <- temp
        
        all_temp <- c(all_temp, temp)
    }

    cols <- alpha(c("red", "blue"),.5)

    ggdata$Round <- as.factor(ggdata$Round)
    ggdata$Donor <- as.factor(ggdata$Donor)

    p <- ggplot(ggdata, aes(x=Round, y=Weight.Gain, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") + geom_quasirandom(dodge.width=.75, size=3) +  
      scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Change in Weight Relative to Baseline") + 
        guides(fill=FALSE) + 
        xlab("Timepoint") + ylab("Weight (g)")
        
    save_plot("mouse-weights.pdf", p, base_aspect_ratio = 1.3)
