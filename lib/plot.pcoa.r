
require(vegan)
require(faraway)
require(RColorBrewer)
plot.pcoa <- function(map, dm, plot.title)
{
    map0 <- map

    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    dm <- dm[valid_samples, valid_samples]

    ddm <- as.dist(dm)
    # plot clusters
    # get PCoA scores first
    pc <- cmdscale(ddm,2)
    
    firstgen <- rownames(map0)[map0$Sample.Group %in% c("Hmong1st","Karen1st")]
    firstgen.years <- cut(map0[firstgen,"Years.in.US"], breaks=8)
    names(firstgen.years) <- firstgen
    
    secondgen <- rownames(map0)[map0$Sample.Group=="Hmong2nd"]
    secondgen.years <- rep("US-born",length(secondgen))
    names(secondgen.years) <- secondgen
    
    thai <- rownames(map0)[map0$Sample.Group %in% c("KarenThai","HmongThai")]
    thai.years <- rep("Thailand",length(thai))
    names(thai.years) <- thai

    residence <- factor(c(firstgen.years, secondgen.years, thai.years), levels=c("Thailand", as.character(1:8), "US-born")) 
    levels(residence)[2:9] <- c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40")

    map0[,"residence"] <- residence[rownames(map0)]    

    # adonis on residence class as factors for ALL samples
    adonis.residence <- adonis(dm ~ map0$residence)
    print(adonis.residence)
    adonis.residence.pval <- adonis.residence$aov.tab[1,"Pr(>F)"]
    # adonis on continuous years.in.us for 1st generation ONLY
    adonis.years <- adonis(dm[firstgen,firstgen] ~ map0[firstgen,"Years.in.US"])
    print(adonis.years)
    adonis.years.pval <- adonis.years$aov.tab[1,"Pr(>F)"]
    
    cols <- alpha(brewer.pal(9,"RdPu"),.7)
    
    cols2 <- vector(mode="character", length=nrow(map0))
    cols2[map0$Years.in.US > 0 & map0$Years.in.US < 5] <- cols[4] #"#fa9fb5"
    cols2[map0$Years.in.US >= 5 & map0$Years.in.US < 10] <- cols[6] #"#99d8c9"
    cols2[map0$Years.in.US >= 10 & map0$Years.in.US < 20] <- cols[7] #"#fa9fb5"
    cols2[map0$Years.in.US >= 20 & map0$Years.in.US < 30] <- alpha(cols[8],.7) #"#f768a1"
    cols2[map0$Years.in.US >= 30] <- alpha(cols[9],.3) #"#ae017e"
    cols2[map0$Years.in.US == 0] <- cols[9] #"#49006a" #magenta
    cols2[is.na(map0$Years.in.US)] <- cols[1] #"#014636" #bluegreen
    names(cols2) <- rownames(map0)
    
    borders <- cols2
    borders[is.na(map0$Years.in.US)] <- cols[9]
    borders[map0$Years.in.US == 0] <- "black"
     
    lookup <- c(21,24) # point type
    names(lookup) <- sort(unique(map0$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    pch2 <- lookup[as.character(map0$Ethnicity)] 
    names(pch2) <- rownames(map0)

    pdf(file=paste0("pcoa - ", plot.title, ".pdf"),useDingbats=F)
    plot(pc[,1], pc[,2], bg=cols2, col=borders, xlab="PC1",ylab="PC2",pch=pch2, main=plot.title)

#     legend(pch=c(21,24,21,24,21), pt.bg=c(cols[1],cols[1],"gray","gray", cols[9]), col=c(cols[9],cols[9],"gray","gray", "black"), x="bottomleft", cex=.8, 
#             legend=c("Hmong-Thailand", "Karen-Thailand", "Hmong-US", "Karen-US","US-born Hmong"))
    legend(pch=c(rep(21,6),24,24,24,24,24), pt.bg=c(cols[c(1,4,6,7,8,9)],"white","white","white","white","gray","gray"), col=c(cols[c(9,4,6,7,8)],"black","white","white","white","white","gray","gray"), x="bottomleft", cex=.7, 
            legend=c("Thailand", "0-5 yrs", "5-10 yrs", "20-30 yrs","> 30 yrs", "US born","","","","","Hmong","Karen"), bty="n", ncol=2, bg="transparent")
    # add adonis pvalue to bottom-right
    legend(x="bottomright",legend=paste0("p=",adonis.years.pval), bty ="n", pch=NA, cex=.7) 
    dev.off()

}

