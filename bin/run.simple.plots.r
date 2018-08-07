
# ALPHA
    # plot alpha as simple line graphs (for presentations only)
    i=2
    m <- map[map$BMI.Class %in% c("Lean","Obese") & map$Ethnicity==c("Hmong","Caucasian"),]
    y <- alphadiv[rownames(m), alpha.metrics[i]]
    d <- cbind(m,y)
    dd <- aggregate(d$y, list(d$Sample.Group, d$BMI.Class), FUN=mean)
    colnames(dd) <- c("Sample.Group", "BMI.Class","y")
    dd[,1] <- as.character(dd[,1])
    dd[,2] <- as.character(dd[,2])
    ret <- plot.line.by.group.x.bmi(map00=dd[,1:2], y=dd$y, ylab=alpha.metrics[i], main="", parametric=TRUE)

### body trends
    # plot everyone, not just those that have been sequenced
    plot.BMI.barplot(map_all[c(karenthai_all,karen_firstgen_cs_all),], bins=0:10, fn="BMI_barplot_Karen.pdf", freq=F)   
    plot.BMI.barplot(map_all[c(hmongthai_all,hmong_firstgen_cs_all,hmong_secondgen_cs_all),], bins=seq(0,45,5), fn="BMI_barplot_Hmong.pdf", freq=F)
    plot.BMI.barplot(map_all[c(karenthai_all,karen_firstgen_cs_all),], bins=0:10, fn="BMI_barplot_Karen_abs.pdf", freq=T)   
    plot.BMI.barplot(map_all[c(hmongthai_all,hmong_firstgen_cs_all,hmong_secondgen_cs_all),], bins=seq(0,45,5), fn="BMI_barplot_Hmong_abs.pdf", freq=T)
    
    plot.BMI.barplot(map_all[cs_all,], bins=seq(0,45,5), fn="BMI_barplot_ALL.pdf", freq=F)   

### BP Ratio
    plot.b.p.barplot(map[c(karenthai,karen_firstgen_cs),], taxa, bins=0:10, fn="b.p.barplot.karen.pdf")
    plot.b.p.barplot(map[cs,], taxa, bins=seq(0,45,5), fn="b.p.barplot.pdf")
