# d = matrix or data.frame of subsetted samples and metadata of interest with alpha div
# plots beeswarm of alpha div by BMI.class and colored by Years.in.US, US-born colored Gray
plot.alphadiv <- function(d)
{
    # experiment with coloring by fraction hitting GG99
    #fraction <- read.table("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/Sequencing/sequences_090716/4. EMBALMER_output/fraction_hitting_GG99.txt", row=1, head=F, quote="", sep=" ")
    #imp_fraction <- fraction[all_samples,1] # load in fraction of hits that hit GG 99 that was previously calculated
    #cols2 <- rbPal(10)[as.numeric(cut(d$BMI, breaks = 10))]
    #names(cols2) <- rownames(d)
    #cols2 <- rbPal(10)[as.numeric(cut(imp_fraction, breaks = 10))]

    # check for normality 
    normality.pval <- shapiro.test(d$alphadiv)$p.value

    if(normality.pval <= 0.05) {
        m <- kruskal.test(d$alphadiv ~ d$BMI.Class)
    } else {
        m <- pairwise.t.test(d$alphadiv, d$BMI.Class, p.adjust.method="bonf")
    }
    print(m)

    rbPal <- colorRampPalette(c('wheat', "red"))
    cols2 <- rbPal(10)[as.numeric(cut(d$Years.in.US, breaks = 10))]
    names(cols2) <- rownames(d)

    lookup <- c(19,17) # point type
    names(lookup) <- sort(unique(d$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    pch2 <- lookup[as.character(d$Ethnicity)] 
    names(pch2) <- rownames(d)

    # change color of US-born
    cols2[d$Years.in.US==0] <- "DarkGray"

    # try a beeswarm plot
    pdf(file="beeswarm_alphadiv_PD.pdf",useDingbats=F)
    beeswarm(d$alphadiv ~ d$BMI.Class, pwcol=cols2, pwpch=pch2, ylim=c(10,45),
        main="Biodiversity by BMI Class", ylab="PD Whole Tree", xlab="")
    bxplot(d$alphadiv ~ d$BMI.Class,add=T)

    legend(pch=c(19,17,19,19,19), col=c("black","black","DarkGray","white", "white"), x="topright", cex=.8, legend=c("Hmong", "Karen", "US-born Hmong", "", ""))
    # add gradient image for age colors
    legend_image <- as.raster(rbPal(10)[10:1])
    text(y=41, x = c(3.01,3.265,3.53), labels = c("0","Yrs in US","40"), cex=c(.6,.7,.6))
    rasterImage(legend_image, 3.0, 40.5, 3.05, 47.5, angle=270)

    dev.off()

}