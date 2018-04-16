library(lme4)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(ggbeeswarm)
library(cowplot)
library(scales)
library(ggsignif)
library(reshape)

source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/mouse.villus.crypt.r")

setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

### load & prep data
    vc_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/all-cpsr-measurements.txt"
    mice_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/All Mouse Data.txt"
    immune_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/immune-cell-counts.txt"
    body_comp_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse_exp/Body Composition.txt"

    mice <- read.table(file=mice_fn, sep="\t", header=T, as.is=T, quote="")
    villus_crypt <- read.table(file=vc_fn, sep="\t", header=T, as.is=T)
    immune <- read.table(file=immune_fn, sep="\t", header=T, as.is=T, check.names=F, quote="", comment="")
    body_comp <- read.table(file=body_comp_fn, sep="\t", header=T, as.is=T)

    mouse.id.levels <- paste0("M", sort(as.numeric(unique(gsub("M", "", mice$Mouse.ID)))))
    mice$Mouse.ID <- factor(mice$Mouse.ID, levels=mouse.id.levels)
    villus_crypt$Mouse.ID <- factor(villus_crypt$Mouse.ID, levels=mouse.id.levels)
    immune$Mouse.ID <- factor(immune$Mouse.ID, levels=mouse.id.levels)
    body_comp$Mouse.ID <- factor(body_comp$Mouse.ID, levels=mouse.id.levels)

    colnames(mice) <- gsub("\\.\\.g\\.", "", colnames(mice))
    colnames(body_comp) <- gsub("\\.\\.g\\.", "", colnames(body_comp))
    
    # get the start and end dates for all mice
    mice$Date <- as.Date(mice$Date, format="%m/%d/%y")
    endpoint <- aggregate(Date ~ Mouse.ID, mice, FUN=max)
    baseline <- aggregate(Date ~ Mouse.ID, mice, FUN=min)
    endpoint$Date <- as.Date(endpoint$Date)
    baseline$Date <- as.Date(baseline$Date)
    endpoint$Mouse.ID <- factor(endpoint$Mouse.ID, levels=mouse.id.levels)
    baseline$Mouse.ID <- factor(baseline$Mouse.ID, levels=mouse.id.levels)

    # transform dates into experiment week
    for(i in 1:nrow(baseline))
    {
      this.mouse <- baseline$Mouse.ID[i]
      mice[mice$Mouse.ID == this.mouse, "Week"] <- (as.numeric(mice[mice$Mouse.ID == this.mouse, "Date"] - baseline$Date[i]))/7
    }
  
    # rename cohoused groups by which donor.type they originally came from
    donor.df <- merge(unique(mice[mice$Donor.Type == "Cohoused","Donor",drop=F]), unique(mice[mice$Week==0,c("Donor","Donor.Type")]), by="Donor")
    donors <- donor.df[,2]
    names(donors) <- donor.df[,1]
    mice[mice$Donor.Type == "Cohoused","Donor.Type"] <- paste(donors[mice[mice$Donor.Type == "Cohoused","Donor"]], mice[mice$Donor.Type == "Cohoused","Donor.Type"], sep=".")
    
    mice$Group <- factor(paste(mice$Donor.Type, mice$Diet.Type, sep="."), 
                    levels=c("Thai.HighFiber", "Thai.LowFiber", "US.HighFiber", "US.LowFiber", "Thai.Cohoused.HighFiber", "US.Cohoused.HighFiber", "Thai.Cohoused.LowFiber", "US.Cohoused.LowFiber"))

    # global vars
    # Thai.HighFiber=darkorange, Thai.LowFiber=lightorange, US.HighFiber=darkblue, US.LowFiber=lightblue, Cohoused.HighFiber=darkgray, Cohoused.LowFiber=lightgray
    GROUP.COLORS <<- c("#F1651D", "#FEA360", "#0A97B7", "#55CFD6", "#92867D", "#DAD5D2", "#949494", "#D6D6D6")
    names(GROUP.COLORS) <- levels(mice$Group)


### VILLI & CRYPT
    # get endpoint measurements and groupings
    vc_end <- merge(villus_crypt, merge(endpoint, mice[!(mice$Cage.ID %in% c("6","7")),c("Mouse.ID","Group","Date")], by=c("Mouse.ID","Date"), all.x=T), by="Mouse.ID")
    plot.vc(vc_end)
    
### WEIGHT
    ggdata <- mice

    # for any mice that didn't survive until Week 8, remove them
    ggdata$Mouse.ID <- as.character(ggdata$Mouse.ID)    
    valid.mouse.ids <- unique(ggdata[ggdata$Week==8, "Mouse.ID"])
    ggdata <- ggdata[ggdata$Mouse.ID %in% valid.mouse.ids,]
  
    mouse.ids <- unique(ggdata$Mouse.ID)
    for(i in 1:length(mouse.ids))
    {
        baseline.weight <- ggdata[ggdata$Week == 0 & ggdata$Mouse.ID == mouse.ids[i], "Body.Weight"]
        # % change
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Grams.Weight.Gain"] <- ggdata[ggdata$Mouse.ID == mouse.ids[i],"Body.Weight"] - baseline.weight
        ggdata[ggdata$Mouse.ID == mouse.ids[i],"Percent.Weight.Gain"] <- 100*(ggdata[ggdata$Mouse.ID == mouse.ids[i],"Grams.Weight.Gain"])/baseline.weight
    }
    weight.all <- ggdata
    weight.all.mean <- aggregate(Percent.Weight.Gain ~ Group + Week, weight.all, mean)
    weight.all.se <- aggregate(Percent.Weight.Gain ~ Group + Week, weight.all, FUN=function(xx) sd(xx)/sqrt(length(xx)))
    weight.all.df <- cbind(weight.all.mean[,1:2], Mean=weight.all.mean$Percent.Weight.Gain, SE=weight.all.se$Percent.Weight.Gain)
    
    p.weight.all <- ggplot(weight.all.df, aes(x=Week, y=Mean, group=Group, color=Group)) + 
            scale_color_manual(name = "", values = GROUP.COLORS) + 
            ylab("Percent Weight Gain") + xlab("Week") + scale_x_continuous(name="Week", breaks=seq(0,8,2)) + 
            geom_line(alpha=.7, size=2) + 
            geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2)
    p.legend <- get_legend(p.weight.all)
    p.weight.all <- p.weight.all + theme(legend.position='none')
    save_plot("mouse-weights-all.pdf", p.weight.all, base_aspect_ratio = 1.3)

    gg_weightwk8 <- ggdata[ggdata$Mouse.ID != "M21" & weight.all$Week %in% c(8,10),]
    # outlier M21 gained 100% of its weight - I think because it was a super tiny mouse and only 4-weeks old 

    p.weight.wk8 <- ggplot(gg_weightwk8, aes(x=Group, y=Percent.Weight.Gain, group=Group, color=Group)) + geom_boxplot(aes(fill=Group), colour="black", width=.9) + 
            scale_fill_manual("Group", values = GROUP.COLORS) +
            ggtitle("% Weight Gain, Week 8") + ylab("") + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                    strip.background = element_blank(),strip.text.y = element_blank()) + theme(legend.position="none")                    
    p.weight.wk8 <- add.pvals.to.plot(p=p.weight.wk8, model=aov(Percent.Weight.Gain ~ Group + Donor, data=gg_weightwk8), group.var="Group", alpha=.15, textsize=3)
    save_plot("mouse-weights-endpoint.pdf", p.weight.wk8, base_aspect_ratio = 1.3)
    
### glucose - endpoints only

    ggdata <- mice[!is.na(mice$Fasting.Glucose),]
    ggdata$Week <- factor(ggdata$Week, labels=c("Baseline", "Post-Diet", "Post-Cohousing"))
    p <- ggplot(ggdata, aes(x=Week, y=Fasting.Glucose, color=Group)) + geom_boxplot(aes(fill = Group), color="black") +
        facet_grid(. ~ Week, scales = "free", switch="x") +
        scale_fill_manual(name="", values = GROUP.COLORS) +
        ggtitle("Fasting Blood Glucose") + 
        ylab("(mg/dL)") + xlab("") +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank()) + theme(legend.position="none")
    save_plot("mouse-glucose-boxplot.pdf", p, base_aspect_ratio = 1.3)

    # plot by fold change
    base.blood <- ggdata[ggdata$Week=="Baseline",c("Mouse.ID","Fasting.Glucose", "Donor")]
    postdiet.blood <- ggdata[ggdata$Week=="Post-Diet",c("Mouse.ID","Fasting.Glucose","Week","Group")]
    postcohouse.blood <- ggdata[ggdata$Week=="Post-Cohousing",c("Mouse.ID","Fasting.Glucose","Week","Group")]

    dfdiet <- merge(postdiet.blood, base.blood, by="Mouse.ID", suffixes=c(".post",".base"))    
    dfcohouse <- merge(postcohouse.blood, base.blood, by="Mouse.ID", suffixes=c(".post",".base"))
    fc.df <- rbind(dfdiet,dfcohouse)
    fc.df$fold.change <- (fc.df$Fasting.Glucose.post - fc.df$Fasting.Glucose.base)/fc.df$Fasting.Glucose.base

    p <- ggplot(fc.df, aes(x=Group, y=fold.change, color=Group)) + geom_boxplot(aes(fill = Group), color="black") +
        scale_fill_manual(name="", values = GROUP.COLORS) +
        ggtitle("Fasting Blood Glucose Fold Change") + 
        ylab("% change") + xlab("") +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank()) + theme(legend.position="none")

    p <- add.pvals.to.plot(p=p, model=aov(fold.change ~ Group + Donor, data=fc.df), group.var="Group", alpha=.15, textsize=3)

    save_plot("mouse-glucose-foldchange.pdf", p, base_aspect_ratio = 1.3)

### food consumption
  chow <- unique(mice[,c("Mouse.ID", "Cage.ID", "Food.Weight.Out", "Food.Weight.In", "Week")])
  weeks <- unique(mice$Week)
  gg_chow <- NULL
  # do this by Mouse.ID because of some issues with adding food to cages after sac'd mice are removed for fasting
  for(i in 2:(length(weeks)))
  {
    last.week <- weeks[i-1]
    this.week <- weeks[i]
    
    cc <- merge(chow[chow$Week==this.week, c("Mouse.ID", "Cage.ID","Food.Weight.Out")], chow[chow$Week==last.week, c("Mouse.ID", "Food.Weight.In")], by=c("Mouse.ID"))
    cc$FoodDiff <- cc$Food.Weight.In - cc$Food.Weight.Out
    cc$Week <- this.week
    
    gg_chow <- rbind(gg_chow, cc)
  }
                  
    # calculate how many mice there were in each cage (alive at the time of food being weighed)
    gg_chow$CageFood <- paste(gg_chow$Cage.ID, gg_chow$FoodDiff, sep=":")
    mice.counts <- table(gg_chow[,c("Week","CageFood")])
    mice.counts.long <- melt(mice.counts, id.vars=(CageFood))
    mice.counts.long <- mice.counts.long[mice.counts.long$value!=0,]
    
    gg_chowc <- merge(mice.counts.long, unique(gg_chow[,c("CageFood", "Cage.ID", "FoodDiff", "Week")]), by=c("CageFood", "Week"))
    
    colnames(gg_chowc)[which(colnames(gg_chowc) == "value")] <- "Num.Mice"

    # now make some manual modifications to account for the mice that died
    #M25  0	0	NA	NA	TFSCS026	Thai	6	TD.86489	LowFiber	NA	NA	11/6/17 --> 2 days before Week 2
    #M26	0	1	NA	NA	IMP.264	US	7	TD.86489	LowFiber	NA	NA	11/2/17 --> 4 days before Week 2
    #M28	0	2	NA	NA	IMP.264	US	7	TD.86489	LowFiber	NA	NA	11/22/17 --> on the end of Week 4

    gg_chowc[gg_chowc$Cage.ID=="6" & gg_chowc$Week==2, "Num.Mice"] <- gg_chowc[gg_chowc$Cage.ID=="6" & gg_chowc$Week==2, "Num.Mice"] + ((14-2)/14) # died 2 days prior
    gg_chowc[gg_chowc$Cage.ID=="7" & gg_chowc$Week==2, "Num.Mice"] <- gg_chowc[gg_chowc$Cage.ID=="6" & gg_chowc$Week==2, "Num.Mice"] + ((14-4)/14) # died 4 days prior
    gg_chowc[gg_chowc$Cage.ID=="7" & gg_chowc$Week==4, "Num.Mice"] <- gg_chowc[gg_chowc$Cage.ID=="7" & gg_chowc$Week==4, "Num.Mice"] + (14/14) # died 0 days prior

    gg_chowc$Food.Per.Mouse <- gg_chowc$FoodDiff/gg_chowc$Num.Mice

    avg_chow <- merge(gg_chowc, unique(mice[,c("Cage.ID","Group")]), by="Cage.ID")

    # manually add starting data points for cohousing lines
    cohouse.wk8 <- rbind(avg_chow[avg_chow$Cage.ID %in% c("1", "2") & avg_chow$Week==8,],
            avg_chow[avg_chow$Cage.ID %in% c("3", "4") & avg_chow$Week==8,])
    cohouse.wk8$Cage.ID <- c(rep("C1-2",2), rep("C3-4",4))
    cohouse.wk8$Group <- c(rep("Cohoused.LowFiber",2), rep("Cohoused.HighFiber",4))
    avg_chow <- rbind(avg_chow, cohouse.wk8)
    
    # cage 5 is an outlier, exclude for now
    # ******** let's not look at the 2nd donor pair until the other diet is available
    # avg_chow <- avg_chow[!(avg_chow$Cage.ID %in% c("5", "6", "7")),]

    p_chow <- ggplot(data=avg_chow, aes(x=Week, y=Food.Per.Mouse, group=Group, color=Group)) + scale_color_manual(name = "Group", values = GROUP.COLORS) + 
            geom_line(data = aggregate(Food.Per.Mouse ~ Group + Week, avg_chow, mean), aes(x=Week, y=Food.Per.Mouse, group=Group, color=Group), alpha=1, size=2) +
            ggtitle("Food Consumed") + ylab("Per Mouse (g)") + xlab("Week") + scale_x_continuous(name="Week", breaks=seq(0,10,2)) + theme(legend.position='none')

    save_plot("mouse-chow-all.pdf", p_chow, base_aspect_ratio = 1.3)
    
### percent weight gain by food consumed
    # A low FCR (food conversion ratio) means that they are more efficient users of food --> less food required to gain weight

    chow.df <- aggregate(Food.Per.Mouse ~ Group + Week, avg_chow, mean) # average food consumed per mouse within group

    # gram individual mouse weights so that we can include error bars in FCR calculation
    weight.food.df <- merge(weight.all[,c("Group","Week","Grams.Weight.Gain", "Donor")], chow.df, by=c("Group","Week"))
    weight.food.df$FCR <- weight.food.df$Food.Per.Mouse/weight.food.df$Grams.Weight.Gain
    
    weight.food.mean <- aggregate(FCR ~ Group + Week, weight.food.df, mean)
    weight.food.se <- aggregate(FCR ~ Group + Week, weight.food.df, FUN=function(xx) sd(xx)/sqrt(length(xx)))
    weight.food.summary.df <- cbind(weight.food.mean[,1:2], Mean.FCR=weight.food.mean$FCR, SE.FCR=weight.food.se$FCR)

    # *** outlier with Week 2 and US.HighFiber, let's drop it
    p.weight.food <- ggplot(weight.food.summary.df[!(weight.food.summary.df$Group == "US.HighFiber" & weight.food.summary.df$Week==2),], aes(x=Week, y=Mean.FCR, group=Group, color=Group)) + 
            scale_color_manual(name = "", values = GROUP.COLORS) + 
            ylab("Food Consumed / Weight Gain") + xlab("Week") + scale_x_continuous(name="Week", breaks=seq(0,8,2)) + 
            ggtitle("Food Conversion Rate") + geom_line(size=2) + geom_errorbar(aes(ymin=Mean.FCR-SE.FCR, ymax=Mean.FCR+SE.FCR), width=.2)
    p.legend <- get_legend(p.weight.food)
    p.weight.food <- p.weight.food + theme(legend.position='none')
    save_plot("mouse-weights-by-foodconsumed.pdf", p.weight.food, base_aspect_ratio = 1.3)
    
    # Week 8 only
    weight.food.wk8 <- weight.food.df[weight.food.df$Week %in% c(8,10),]
    p.weight.food.wk8 <- ggplot(weight.food.wk8, aes(y=FCR, fill=Group, x=Group)) + scale_fill_manual(name = "", values = GROUP.COLORS) + 
            ylab("Food Consumed / Weight Gain") + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
            ggtitle("Food Conversion Rate, endpoints") + geom_boxplot(color="black") + theme(legend.position='none')
    p.weight.food.wk8 <- add.pvals.to.plot(p=p.weight.food.wk8, model=aov(FCR ~ Group + Donor, data=weight.food.wk8), group.var = "Group", alpha=.15, textsize=3)
    
    save_plot("mouse-weights-by-foodconsumed-endpoints.pdf", p.weight.food.wk8, base_aspect_ratio = 1.3)

### body composition
    
    gg_bodycomp <- merge(body_comp[,c("Mouse.ID","Percent.Fat","Percent.Lean")], endpoint, by="Mouse.ID")
    gg_bodycomp <- merge(gg_bodycomp, mice[, c("Mouse.ID", "Date", "Group", "Donor", "Donor.Type", "Cage.ID", "Diet")], by=c("Mouse.ID", "Date"))
    
    p_lean <- ggplot(gg_bodycomp, aes(x=Group, y=Percent.Lean, fill=Group)) + geom_boxplot(color="black", width=.9) + scale_fill_manual(values=GROUP.COLORS) +
      ggtitle("% Lean Mass") + ylab("") + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
      theme(legend.position="none")
    p_lean <- add.pvals.to.plot(p=p_lean, model=aov(Percent.Lean ~ Group + Donor, data=gg_bodycomp), group.var="Group", alpha=.15, textsize=3)
    save_plot("mouse-percent-lean.pdf", p_lean, base_aspect_ratio = 1.3)
    
    p_fat <- ggplot(gg_bodycomp, aes(x=Group, y=Percent.Fat, fill=Group)) + geom_boxplot(color="black", width=.9) + scale_fill_manual(values=GROUP.COLORS) +
            ggtitle("% Body Fat") + ylab("") + xlab("") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
            theme(legend.position="none")
    # build model that controls for the specific donor, THEN do a post-hoc pairwise comparison test to see what group is significantly different
    p_fat <- add.pvals.to.plot(p=p_fat, model=aov(Percent.Fat ~ Group + Donor, data=gg_bodycomp), group.var="Group", alpha=.15, textsize=3)
    save_plot("mouse-percent-fat.pdf", p_fat, base_aspect_ratio = 1.3)


### immune cell types
    # note that counts reported are per inch (total counts will need to be recalculated)
    
    base.cell.type <- "Count.Live.CD45+" #"Count.CD3e+"
    
    count.cell.names <- grep("^Count", colnames(immune), value=T)
    
    IEL <- immune[immune$Cell.Population=="IEL",]
        
    # divide total count of cell types by total counts of CD45 (basic immune cell types) instead of using "per inch" estimates (biased)
    ratio.names <- paste0(gsub("Count\\.", "", count.cell.names))    
    IEL[,ratio.names] <- IEL[,count.cell.names]/IEL[,base.cell.type]

    end.groups <- merge(endpoint, mice[,c("Mouse.ID","Group","Date")], by=c("Mouse.ID","Date"), all.x=T)
    end.groups$Mouse.ID <- as.character(end.groups$Mouse.ID)
    gg_iel <- merge(IEL, end.groups, by="Mouse.ID")
    
    p <- list()
    cell.types <- ratio.names
    for(i in 1:length(cell.types))
    {
        gd <- gg_iel[,c(cell.types[i], "Mouse.ID", "Group")]
        colnames(gd)[1] <- "Cell.Type"
        p[[i]] <- ggplot(gd, aes(x=Group, y=Cell.Type, group=Group, color=Group)) + geom_boxplot(aes(fill=Group), colour="black", width=.9) + 
            scale_fill_manual(values = GROUP.COLORS) + 
            ggtitle(cell.types[i]) + ylab("") + xlab("") + 
            theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), strip.background = element_blank(),
                strip.text.x = element_blank(), legend.position="none", plot.title = element_text(size=10))

        p[[i]] <- add.pvals.to.plot(p=p[[i]], model=aov(gd$Cell.Type ~ gd$Group), group.var="gd$Group")

    }
    save_plot(paste0("immune-cells-per-", base.cell.type,".pdf"), plot_grid(plotlist=p, ncol=3), base_aspect_ratio = 1.8)
    
# if there are pairwise comparisons that are significant, add them to the plot
add.pvals.to.plot <- function(p, model, group.var, alpha=.10, textsize=1.8)
{
    tuke <- TukeyHSD(model)[[group.var]][,"p adj"] 
    sig.pvals <- signif(tuke[tuke < alpha],2)
    comparisons <- strsplit(names(sig.pvals),"-")

    # if there are significant comparisons, show them
    if(length(sig.pvals) > 0) 
        p <- p + geom_signif(comparisons = comparisons, annotation=as.character(sig.pvals), color="black", tip_length=.01, textsize=textsize, step_increase=.08)

    p
}

