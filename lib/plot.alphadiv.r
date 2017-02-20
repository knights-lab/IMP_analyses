# d = matrix or data.frame of subsetted samples and metadata of interest with alpha div
# plots beeswarm of alpha div by BMI.class and colored by Years.in.US, US-born colored Gray
require(beeswarm)
require(faraway)
require(RColorBrewer)
plot.alphadiv <- function(map_alpha) #takes in a dataframe of mapping file + alpha div metric
{
    map0 <- map_alpha
    
    # check for normality 
    normality.pval <- shapiro.test(map0$alphadiv)$p.value

    if(normality.pval <= 0.05) {
        m <- kruskal.test(map0$alphadiv ~ map0$BMI.Class)
    } else {
        m <- pairwise.t.test(map0$alphadiv, map0$BMI.Class, p.adjust.method="bonf")
    }
    print(m)

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

    # try a beeswarm plot
    pdf(file="beeswarm_alphadiv_PD.pdf",useDingbats=F)
    beeswarm(map0$alphadiv ~ map0$BMI.Class, pwbg=cols2, pwcol=borders, pwpch=pch2,
        main="Biodiversity by BMI Class", ylab="PD Whole Tree", xlab="")
    bxplot(map0$alphadiv ~ map0$BMI.Class,add=T)

    legend(pch=c(21,24,21,24,21), pt.bg=c(cols[1],cols[1],"gray","gray", cols[9]), col=c(cols[9],cols[9],"gray","gray", "black"), x="topright", cex=.7, 
            legend=c("Hmong-Thailand", "Karen-Thailand", "Hmong-US", "Karen-US","US-born Hmong"))

#     legend(pch=c(19,17,19,19,17,19,19), col=c("black","black","blue","green", "green", "white","white"), x="topright", cex=.8, legend=c("Hmong-US", "Karen-US", "US-born Hmong", "Hmong-Thailand", "Karen-Thailand","",""))
#     # add gradient image for age colors
#     legend_image <- as.raster(rbPal(10)[10:1])
#     text(y=58, x = c(3.01,3.265,3.53), labels = c("0","Yrs in US","40"), cex=c(.6,.7,.6))
#     rasterImage(legend_image, 3.15, 60, 3.2, 65, angle=270)

    dev.off()

}