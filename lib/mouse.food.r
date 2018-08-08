# calculates the average food eaten per mouse based on all of the different cages
calc.food <- function(ggdata)
{
    chow <- unique(ggdata[,c("Mouse.ID", "Cage.ID", "Food.Weight.Out", "Food.Weight.In", "Week")])
    weeks <- unique(ggdata$Week)
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
    map.counts <- table(gg_chow[,c("Week","CageFood")])
    map.counts.long <- melt(map.counts, id.vars=(CageFood))
    map.counts.long <- map.counts.long[map.counts.long$value!=0,]

    gg_chowc <- merge(map.counts.long, unique(gg_chow[,c("CageFood", "Cage.ID", "FoodDiff", "Week")]), by=c("CageFood", "Week"))

    colnames(gg_chowc)[which(colnames(gg_chowc) == "value")] <- "Num.Mice"

    # now make some manual modifications to account for the mice that died
    #M25  0	0	NA	NA	TFSCS026	Thai	6	TD.86489	LowFiber	NA	NA	11/6/17 --> 2 days before Week 2
    #M26	0	1	NA	NA	IMP.264	US	7	TD.86489	LowFiber	NA	NA	11/2/17 --> 4 days before Week 2
    #M28	0	2	NA	NA	IMP.264	US	7	TD.86489	LowFiber	NA	NA	11/22/17 --> on the end of Week 4

    gg_chowc[gg_chowc$Cage.ID=="6" & gg_chowc$Week==2, "Num.Mice"] <- gg_chowc[gg_chowc$Cage.ID=="6" & gg_chowc$Week==2, "Num.Mice"] + ((14-2)/14) # died 2 days prior
    gg_chowc[gg_chowc$Cage.ID=="7" & gg_chowc$Week==2, "Num.Mice"] <- gg_chowc[gg_chowc$Cage.ID=="6" & gg_chowc$Week==2, "Num.Mice"] + ((14-4)/14) # died 4 days prior
    gg_chowc[gg_chowc$Cage.ID=="7" & gg_chowc$Week==4, "Num.Mice"] <- gg_chowc[gg_chowc$Cage.ID=="7" & gg_chowc$Week==4, "Num.Mice"] + (14/14) # died 0 days prior

    gg_chowc$Food.Per.Mouse <- gg_chowc$FoodDiff/gg_chowc$Num.Mice
    
    avg_chow <- merge(gg_chowc, unique(ggdata[,c("Cage.ID","Group")]), by="Cage.ID")

    # manually add starting data points for cohousing lines
#     cohouse.wk8 <- rbind(avg_chow[avg_chow$Cage.ID %in% c("1", "2") & avg_chow$Week==8,],
#         avg_chow[avg_chow$Cage.ID %in% c("3", "4") & avg_chow$Week==8,])
#     cohouse.wk8$Cage.ID <- c(rep("C1-2",2), rep("C3-4",4))
#     cohouse.wk8$Group <- c(rep("Cohoused.LowFiber",2), rep("Cohoused.HighFiber",4))
#     avg_chow <- rbind(avg_chow, cohouse.wk8)

    return(avg_chow)
}

plot.food <- function(ggdata, group)
{
    ggdata$Group <- ggdata[,group]
    avg_chow <- calc.food(ggdata)
    
    p_chow <- ggplot(data=avg_chow, aes(x=Week, y=Food.Per.Mouse, group=Group, color=Group)) + scale_color_manual(name = "Group", values = GROUP.COLORS) + 
        geom_line(data = aggregate(Food.Per.Mouse ~ Group + Week, avg_chow, mean), aes(x=Week, y=Food.Per.Mouse, group=Group, color=Group), alpha=1, size=2) +
        geom_point(data = aggregate(Food.Per.Mouse ~ Group + Week, avg_chow, mean), aes(shape=Group), color="black", size=3) + scale_shape_manual(name="Group",values=GROUP.SHAPES) +
        ggtitle("Average Food Consumed") + ylab("Per Mouse (g)") + xlab("Week") + scale_x_continuous(name="Week", breaks=seq(0,10,2)) + theme(legend.position='none')

    save_plot("food-consumption-L.pdf", p_chow, base_aspect_ratio = 1)
}

plot.feed.efficiency.L <- function(ggdata, group)
{
    ggdata$Group <- ggdata[,group]
    avg_chow <- calc.food(ggdata)

    weight.all <- calc.weights(ggdata)

    chow.df <- aggregate(Food.Per.Mouse ~ Group + Week, avg_chow, mean) # average food consumed per mouse within group

    # grab individual mouse weights so that we can include error bars in FE calculation
    weight.food.df <- merge(weight.all[,c("Mouse.ID","Cohoused","Group","Week","Grams.Weight.Gain", "Donor")], chow.df, by=c("Group","Week"))
    weight.food.df$FE <- weight.food.df$Grams.Weight.Gain/weight.food.df$Food.Per.Mouse
    
    weight.food.mean <- aggregate(FE ~ Group + Week, weight.food.df, mean)
    weight.food.se <- aggregate(FE ~ Group + Week, weight.food.df, FUN=function(xx) sd(xx)/sqrt(length(xx)))
    weight.food.summary.df <- cbind(weight.food.mean[,1:2], Mean.FE=weight.food.mean$FE, SE.FE=weight.food.se$FE)

    p.weight.food <- ggplot(weight.food.summary.df, aes(x=Week, y=Mean.FE, group=Group, color=Group)) + 
            scale_color_manual(name = "Group", values = GROUP.COLORS) + 
            ylab("Weight Gain / Food Consumed") + xlab("Week") + scale_x_continuous(name="Week", breaks=seq(0,10,2)) + 
            ggtitle("Feed Efficiency\n") + geom_line(size=2) + geom_errorbar(aes(ymin=Mean.FE-SE.FE, ymax=Mean.FE+SE.FE), width=.2) +
            geom_point(aes(shape=Group), color="black", size=3) + scale_shape_manual(name="Group",values=GROUP.SHAPES)

    p.legend <- get_legend(p.weight.food)
    p.weight.food <- p.weight.food + theme(legend.position='none')
    save_plot("feed-efficiency-L.pdf", p.weight.food, base_aspect_ratio = 1)
}

# use group = group.end instead of start
plot.feed.efficiency.wk8 <- function(ggdata, group, add.pval=TRUE, outputfn)
{
    ggdata$Group <- ggdata[,group]
    avg_chow <- calc.food(ggdata)

    weight.all <- calc.weights(ggdata)

    chow.df <- aggregate(Food.Per.Mouse ~ Group + Week, avg_chow, mean) # average food consumed per mouse within group

    # grab individual mouse weights so that we can include error bars in FE calculation
    weight.food.df <- merge(weight.all, chow.df, by=c("Group","Week"))
    weight.food.df$FE <- weight.food.df$Grams.Weight.Gain/weight.food.df$Food.Per.Mouse
    
    # Week 8 only
    weight.food.wk8 <- weight.food.df[weight.food.df$Week %in% c(8,10),]
    p.weight.food.wk8 <- mouse.boxplot(y=weight.food.wk8$FE, Group=weight.food.wk8$Group, 
            main="Feed Efficiency\n", add.pval=add.pval,
            ylab="Weight Gain / Food Consumed", group.vars.df=weight.food.wk8[,c("Donor.Type","Diet.Type")], hide.box=FALSE)

    save_plot(outputfn, p.weight.food.wk8, base_aspect_ratio = 1)
}