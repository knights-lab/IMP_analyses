#OLD
plot.alphadiv <- function(map00, alpha, metric) 
{
    map0 <- data.frame(map00, alphadiv=alpha[,metric], stringsAsFactors=F)
    
    d <- map0[,c("alphadiv", "BMI.Class")]
    d[map0$Sample.Group %in% c("KarenThai", "HmongThai"), "Group"] <- "Thai"
    d[map0$Sample.Group %in% c("Karen1st", "Hmong1st"), "Group"] <- "1st-Gen"
    d[map0$Sample.Group == "Hmong2nd", "Group"] <- "2nd-Gen"
    d[map0$Sample.Group == "Control", "Group"] <- "Control"
    
    d$Group <- factor(d$Group, levels=c("Thai","1st-Gen","2nd-Gen","Control"))
    d$Group <- factor(d$Group) # remove any levels that aren't present
    
    # check for normality, then do the appropriate test
    groups <- levels(d$Group)
    for(i in 1:length(groups))
    {
        m <- test.samples(d[d$Group==groups[i], "alphadiv"], d[d$Group==groups[i], "BMI.Class"])
        print(paste0("Testing Diversity in Group: ", groups[i]))
        print(m)
    }
    
    d$all.groups <- paste(d$Group,d$BMI.Class,sep=".")
    print(test.samples(d$alphadiv, as.factor(d$all.groups)))

    g <- expand.grid(levels(d$Group), levels(d$BMI.Class))
    

    p <- ggplot(d, aes(Group, alphadiv)) + geom_boxplot(aes(fill = BMI.Class)) +  geom_point(aes(y=alphadiv, fill = BMI.Class), position=position_dodge(width=.75), color=alpha("black",.3), shape=1) +
        scale_fill_manual(name = "Sample Groups", values = c("#ffffcc","#a1dab4","#41b6c4")) +      #, labels = c("0" = "T", "1" = "Bar"))
        ggtitle("Alpha diversity by sample group") + 
        ylab(metric) + xlab("BMI Class") + theme(axis.text.x = element_text(size=10))

    save_plot(paste0("alphadiv_", metric, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}