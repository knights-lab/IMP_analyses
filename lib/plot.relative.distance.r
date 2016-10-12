# plots relative distance of foreign-born samples to US born samples
plot.relative.distance <- function(dm, imp_samples, US_born_samples, filename)
{
    # compare IMP against baseline of US-born Hmong
    IMP_to_US <- NULL
    for (i in 1:length(imp_samples))
    {
        IMP_to_US[i] <- mean(as.numeric(dm[imp_samples[i], US_born_samples]))
    }
    names(IMP_to_US) <- imp_samples
    
    mean_data <- cbind(Distance.to.US=IMP_to_US, map[imp_samples,], stringsAsFactors=F)

    # let's reset Years in US for the US-born to -1 so we can plot it too
    mean_data$Years.in.US[mean_data$Years.in.US==0] <- -1

    lookup <- c(19,17) # point type
    names(lookup) <- sort(unique(mean_data$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    final.pch <- lookup[as.character(mean_data$Ethnicity)] 

    rbPal <- colorRampPalette(c('yellow','blue'))
    final.cols <- rbPal(10)[as.numeric(cut(mean_data$BMI, breaks = 10))]

    pdf(file=filename,useDingbats=F)
    plot(mean_data$Years.in.US, mean_data$Distance.to.US, xlab="Years in US", ylab="Distance to US-born Hmong (Unweighted Unifrac)", main="Foreign-born Hmong and Karen \nbecome similar to US-born over time"
        , col=final.cols, pch=final.pch)
    m <- lm(mean_data$Distance.to.US~mean_data$Years.in.US)
    abline(m)

    legend(pch=c(19,17), col=c("black","black","white","white"), x="topright", cex=.8, legend=c("Hmong", "Karen", "", ""))
    # add gradient image for age colors
    legend_image <- as.raster(rbPal(10)[10:1])
    text(y=.553, x = c(14.2,15,15.8), labels = c("17","BMI","39"), cex=c(.6,.8,.6))
    rasterImage(legend_image, 14.1, .55, 14.3, .569, angle=270)

    # plot Rsquared and p value
    s <- summary(m)
    rsquared <- s$r.squared
    pval <- coef(s)[2,4]
    text(x=14.2, y=c(.47,.474), labels=c(paste("P =",round(pval,4)), paste("R2 =",round(rsquared,4))), cex=.7)
    dev.off()
}