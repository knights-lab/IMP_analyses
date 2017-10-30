library(ggplot2)
library(RColorBrewer)
require(ggbeeswarm)
library(cowplot)

setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

mice_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/mice.txt"
chow_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/chow.txt"
vc_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/villus-crypt.txt"
immune_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/immune-cell-counts.txt"


mice <- read.table(file=mice_fn, sep="\t", header=T, as.is=T)
chow <- read.table(file=chow_fn, sep="\t", header=T, as.is=T)
villus_crypt <- read.table(file=vc_fn, sep="\t", header=T, as.is=T)
immune <- read.table(file=immune_fn, sep="\t", header=T, as.is=T, check.names=F)

# rename colnames
colnames(mice)[which(colnames(mice) %in% c("Fasting.Glucose..mg.dL.", "Weight..g."))] <- c("Fasting.Glucose", "Weight")
colnames(chow)[which(colnames(chow) %in% c("Feed.Start.Weight..g.", "Feed.End.Weight..g.", "Feed.Consumed..g."))] <- c("Feed.Start", "Feed.End", "Feed.Diff")

### weights
    mouse.ids <- unique(mice$Mouse.ID)

    baseline <- mice[mice$Week==0,]
    ggdata <- mice

#    all_temp <- NULL
    for(i in 1:length(mouse.ids))
    {
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Weight.Gain"] <- ggdata[ggdata$Mouse.ID == mouse.ids[i],"Weight"] - baseline[baseline$Mouse.ID == mouse.ids[i], "Weight"]
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Glucose.Diff"] <- ggdata[ggdata$Mouse.ID == mouse.ids[i],"Fasting.Glucose"] - baseline[baseline$Mouse.ID == mouse.ids[i], "Fasting.Glucose"]
 
        # for debugging purposes only
#         temp <- ggdata[ggdata$Mouse.ID == mouse.ids[i],"Weight"] - baseline[baseline$Mouse.ID == mouse.ids[i], "Weight"]
#         names(temp) <- paste( ggdata[ggdata$Mouse.ID == mouse.ids[i], "Mouse.ID"], ggdata[ggdata$Mouse.ID == mouse.ids[i], "Week"], sep=".")
#         all_temp <- c(all_temp, temp)      
    }
    cols <- alpha(c("#869c69", "#cb0323", "#72553f"),.5)

    ggdata$Week.Factor <- as.factor(ggdata$Week)
    ggdata$Donor <- as.factor(ggdata$Donor)

    p <- ggplot(ggdata[ggdata$Week != 0,], aes(x=Week.Factor, y=Weight.Gain, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") + geom_quasirandom(dodge.width=.75, size=3) +  
      scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Weight Gain") + 
        guides(fill=FALSE) + ylab("Change in Body Weight (g)") + xlab("Week")
 
    save_plot("mouse-weights-boxplot.pdf", p, base_aspect_ratio = 1.3)

    #try plotting weights as individual line plots?
     p <- ggplot(data=ggdata[ggdata$Week != 0,], aes(x=Week, y=Weight.Gain, group=Mouse.ID, color=Donor)) + geom_line() + scale_color_manual(values=cols) + 
         scale_x_continuous(name="Week", breaks=seq(0,10,2)) + ylab("Change in Body Weight (g)") + ggtitle("Body Weight Gain") + geom_point(size=4)

    save_plot("mouse-weights-line.pdf", p, base_aspect_ratio = 1.3)


### glucose

    p <- ggplot(ggdata[ggdata$Glucose.Diff != 0 & !is.na(ggdata$Glucose.Diff),], aes(x=Week.Factor, y=Glucose.Diff, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") +
        geom_quasirandom(dodge.width=.75, size=3) +  
        scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Fasting Blood Glucose Change") + 
        guides(fill=FALSE) + ylab("Fasting Blood Glucose (mg/dL)") + xlab("")

    p <- ggdraw(p) + draw_figure_label(paste0("Week 5, P = ",signif(t.test(ggdata[!is.na(ggdata$Glucose) & ggdata$Week=="8", "Glucose.Diff"] ~ ggdata[!is.na(ggdata$Glucose) & ggdata$Week=="8", "Donor"])$p.value, 2)), 
                                        size=10, position="bottom.right")
        
    save_plot("mouse-glucose-boxplot.pdf", p, base_aspect_ratio = 1.3)

     p <- ggplot(ggdata[!is.na(ggdata$Fasting.Glucose),], aes(x=Week, y=Fasting.Glucose, group=Mouse.ID, color=Donor)) + geom_line() + scale_color_manual(values=cols) + 
         scale_x_continuous(name="Week", breaks=seq(0,10,2)) + ylab("Fasting Blood Glucose (mg/dL)") + ggtitle("Fasting Blood Glucose") + geom_point(size=4)

    save_plot("mouse-glucose-line.pdf", p, base_aspect_ratio = 1.3)


### chow
    ggdata2 <- chow
    ggdata2[ggdata2$Donor!="Thailand.US.D1", "Avg.Feed.Diff"] <- ggdata2[ggdata2$Donor!="Thailand.US.D1", "Feed.Diff"]/5 
    ggdata2[ggdata2$Donor=="Thailand.US.D1", "Avg.Feed.Diff"] <- ggdata2[ggdata2$Donor=="Thailand.US.D1", "Feed.Diff"]/4 
    ggdata2$Donor <- factor(ggdata2$Donor, levels=c("Thailand.D1", "US.D1", "Thailand.US.D1"))
    p <- ggplot(ggdata2, aes(x=Week, y=Avg.Feed.Diff, group=Donor, color=Donor)) + geom_line() + scale_color_manual(values=cols) +
        ggtitle("Food Intake") + ylab("Average Food Consumed per Mouse (g)") + geom_point(size=4)

    save_plot("mouse-chow-weights.pdf", p, base_aspect_ratio = 1.3)
    
    
### villus heights
    villus_crypt$Mouse.ID <- factor(villus_crypt$Mouse.ID, levels = c(paste0("M", c(3,4,5,6,7,10))))
    villus <- villus_crypt[villus_crypt$Type=="villus",]
    p <- ggplot(villus, aes(x=Mouse.ID, y=Measurement, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") +
        geom_quasirandom(dodge.width=.75, size=3) +  
        scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Villus Height, Duodenum") + 
        guides(fill=FALSE) + ylab("Villus Height (um)") + xlab("")


    vc.mean <- aggregate(villus$Measurement,list(villus$Mouse.ID), mean)
    colnames(vc.mean) <- c("Mouse.ID","mean")
    vc.mean <- cbind(vc.mean, Donor=unique(crypts[,c("Mouse.ID","Donor")])$Donor)
     p <- ggdraw(p) + draw_figure_label(paste0("P = ",signif(t.test(vc.mean$mean ~ vc.mean$Donor)$p.value, 2)), 
                                         size=10, position="bottom.right")
                
    save_plot("mouse-villi.pdf", p, base_aspect_ratio = 1.3)

### crypt depths

    crypts <- villus_crypt[villus_crypt$Type=="crypt",]
    p <- ggplot(crypts, aes(x=Mouse.ID, y=Measurement, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") +
        geom_quasirandom(dodge.width=.75, size=3) +  
        scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Crypt Depth, Duodenum") + 
        guides(fill=FALSE) + ylab("Crypt Depth (um)") + xlab("")

    vc.mean <- aggregate(crypts$Measurement,list(crypts$Mouse.ID), mean)
    colnames(vc.mean) <- c("Mouse.ID","mean")
    vc.mean <- cbind(vc.mean, Donor=unique(crypts[,c("Mouse.ID","Donor")])$Donor)
     p <- ggdraw(p) + draw_figure_label(paste0("P = ",signif(t.test(vc.mean$mean ~ vc.mean$Donor)$p.value, 2)), 
                                         size=10, position="bottom.right")
        
    save_plot("mouse-crypts.pdf", p, base_aspect_ratio = 1.3)

### villus-to-crypt ratios

    vc <- cbind(villus_crypt[villus_crypt$Type=="crypt",c("Mouse.ID","Sample.ID","Donor")], crypt=villus_crypt[villus_crypt$Type=="crypt","Measurement"], villus=villus_crypt[villus_crypt$Type=="villus","Measurement"])
    vc[,"vc.ratio"] <- vc$villus/vc$crypt
    p <- ggplot(vc, aes(x=Mouse.ID, y=vc.ratio, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") +
        geom_quasirandom(dodge.width=.75, size=3) +  
        scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Villus-to-Crypt Ratio") + 
        guides(fill=FALSE) + ylab("") + xlab("")

    vc.mean <- aggregate(vc$vc.ratio,list(vc$Mouse.ID), mean)
    colnames(vc.mean) <- c("Mouse.ID","mean.vc.ratio")
    vc.mean <- cbind(vc.mean, Donor=unique(vc[,c("Mouse.ID","Donor")])$Donor)
     p <- ggdraw(p) + draw_figure_label(paste0("P = ",signif(t.test(vc.mean$mean.vc.ratio ~ vc.mean$Donor)$p.value, 2)), 
                                         size=10, position="bottom.right")
    save_plot("villus-crypt-ratio.pdf", p, base_aspect_ratio = 1.3)

### R01 figure 5

    LPL <- colnames(immune)[9:13]
    immune.long <- melt(immune)
    data.LPL <- immune.long[immune.long$variable %in% LPL,]

    pvals<-NULL
    for(i in 1:length(LPL))
        pvals[i] <- paste0("P=", signif(t.test(immune[,LPL[i]] ~ immune$Donor)$p.value, 2))

    ##immune
    p <- ggplot(data.LPL, aes(x=variable, y=value, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") +
        geom_quasirandom(dodge.width=.75, size=3) +  
        scale_color_manual(name = "Donor", values = cols, labels=c("Thai Donor","US Donor")) + # set color for points from quasirandom
        ggtitle("Lymphocyte % in SI lamina propria compartment") + 
        guides(fill=FALSE) + ylab("") + xlab("") + theme(legend.title=element_blank(), 
        legend.text = element_text(size=8), 
        axis.text=element_text(size=7), plot.title = element_text(size=11)) +
        scale_x_discrete(labels = paste0(LPL, "\n", pvals))

    p.legend <- get_legend(p)
    p <- p + theme(legend.position="none")
    ## weight
    p2 <- ggplot(ggdata[ggdata$Week == 8,], aes(x=Week.Factor, y=Weight.Gain, color=Donor)) + geom_boxplot(aes(fill = Donor), alpha=0, colour="black") + geom_quasirandom(dodge.width=.75, size=3) +  
      scale_color_manual(name = "Donor", values = cols) + # set color for points from quasirandom
        ggtitle("Weight Gain (g)")  + theme(legend.position="none") + 
        guides(fill=FALSE) + ylab("") + xlab("") + theme(plot.title = element_text(size=10), axis.text=element_text(size=8),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

    p2 <- ggdraw(p2) + draw_figure_label(paste0("P = ",signif(t.test(ggdata[ggdata$Week=="8", "Weight.Gain"] ~ ggdata[ggdata$Week=="8", "Donor"])$p.value, 2)), size=10, position="bottom.right")
   
    ## food
    food.data <- ggdata2[ggdata2$Donor %in% c("Thailand.D1", "US.D1"),]    
    p3 <- ggplot(food.data, aes(x=Week, y=Avg.Feed.Diff, group=Donor, color=Donor)) + geom_line() + scale_color_manual(values=cols) +
        ggtitle("Average Food Intake (g)") + ylab("") + geom_point(size=4) + theme(legend.position="none", axis.text=element_text(size=8),
        axis.title=element_text(size=10), plot.title = element_text(size=10))

    pp <- plot_grid(p, plot_grid(p2, p3, p.legend, ncol=3, labels=c("B","C", ""), rel_widths=c(1,1.1,.5)), labels=c("A",""), ncol=1, rel_heights=c(1.1,1))
    
    #pp <- ggdraw(pp) + draw_label("Figure 5", x = 1, y = 0, hjust=1.3, vjust=-.5, size = 12)

    save_plot("Figure 5 - Mouse Exp.pdf", pp, base_aspect_ratio=1.2)
    

