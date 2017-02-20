require(faraway)
# plot BMI and Waist-Height over years
plot.body.trends <- function(map)
{
    # remove Thailand samples
    # remove US-born samples
    map <- map[map$Years.in.US > 0 & !is.na(map$Years.in.US),]

    plookup <- c(19,17) # point type
    names(plookup) <- sort(unique(map$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    pch2 <- plookup[as.character(map$Ethnicity)] 
    names(pch2) <- rownames(map)

   pdf(file="BodyTrends_immigrants.pdf",useDingbats=F)
   par(mfrow=c(2,2))
    
 
    m <- lm(map$BMI ~ map$Years.in.US)
    s <- summary(m)
    rsquared <- s$r.squared
    pval <- coef(s)["map$Years.in.US",4]

    plot(map$BMI ~ map$Years.in.US, col=alpha("black",.5), ylab="BMI", xlab="Years in the US", pch=pch2,
                main=paste("P =",round(pval,6),"R2 =",round(rsquared,4)))
                
    abline(m, col="gray")


    m <- lm(map$Waist.Height.Ratio ~ map$Years.in.US)
    s <- summary(m)
    rsquared <- s$r.squared
    pval <- coef(s)["map$Years.in.US",4]



    plot(map$Waist.Height.Ratio ~ map$Years.in.US, col=alpha("black",.5), ylab="Waist-to-Height Ratio", xlab="Years in the US", pch=pch2,
                main=paste("P =",round(pval,6),"R2 =",round(rsquared,4)))
                
    abline(m, col="gray")

    dev.off()
}