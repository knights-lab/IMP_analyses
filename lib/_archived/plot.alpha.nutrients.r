# this just plots nutrients by BMI class -- consider getting rid of this
plot.alphadiv.nutrients <- function(map0, fn, nutrients0, nutrient.name="g Fiber per 1000 Calories") #takes in a dataframe of mapping file + alpha div metric
{
    d <- map0[,c("alphadiv", "BMI.Class")]
    d[map0$Sample.Group %in% c("KarenThai", "HmongThai"), "Group"] <- "Thai"
    d[map0$Sample.Group %in% c("Karen1st", "Hmong1st"), "Group"] <- "1st-Gen"
    d[map0$Sample.Group == "Hmong2nd", "Group"] <- "2nd-Gen"
    
    valid_samples <- intersect(rownames(d), rownames(nutrients0))
    d <- d[valid_samples,]
    nutrient <- log10(nutrients0[valid_samples,nutrient.name]+1)

    d<-cbind(d,nutrient)
    #plot(d$alphadiv,nutrient)
    
    d$Group <- factor(d$Group, levels=c("Thai","1st-Gen","2nd-Gen"))
    d$Group <- factor(d$Group) # remove any levels that aren't present
    
    # check for normality, then do the appropriate test
    groups <- levels(d$Group)
    for(i in 1:length(groups))
    {
        m <- test.samples(d[d$Group==groups[i], "alphadiv"], d[d$Group==groups[i], "BMI.Class"])
        print(paste0("Testing Diversity in Group: ", groups[i]))
        print(m)
    }

    p <- ggplot(d, aes(Group, nutrient)) + geom_boxplot(aes(fill = BMI.Class)) +  geom_point(aes(y=nutrient, fill = BMI.Class), position=position_dodge(width=.75), color=alpha("black",.3), shape=1) +
        scale_fill_manual(name = "Sample Groups", values = c("#ffffcc","#a1dab4","#41b6c4")) +      #, labels = c("0" = "T", "1" = "Bar"))
        ggtitle("Nutrient consumption by BMI Class") + 
        ylab(nutrient.name) + xlab("BMI Class")

    save_plot(fn, p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}