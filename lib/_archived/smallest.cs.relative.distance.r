# can't remember what this does....
# plot relative distance of cross sectional samples to a set of reference samples
find.smallest.relative.CS <- function(map, dm, main, query_samples, ref_samples, outputfn, ylab, saveplot=T)
{
    map0 <- map
    map <- map[query_samples,]

#     rel.distance.ref <- NULL
#     for (i in 1:length(ref_samples))
#     {
#         rel.distance.ref[i] <- mean(as.numeric(dm[ref_samples[i], ref_samples]))
#     }


    rel.distance <- NULL
    for (i in 1:length(query_samples))
    {
        rel.distance[i] <- mean(as.numeric(dm[query_samples[i], ref_samples]))
    }
    names(rel.distance) <- query_samples
    
    # let's reset Years in US for the Thailand samples as 0 for plotting purposes
    map$Years.in.US[is.na(map$Years.in.US)] <- -1

    lookup <- c(19,17) # point type
    names(lookup) <- sort(unique(map$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    final.pch <- lookup[as.character(map$Ethnicity)] 

#    cols <- "black"
    col.lookup <- c("#1a9641","#fecc5c","#d7191c")
    names(col.lookup) <- sort(unique(map$BMI.Class)) # let's Hmong to solid filled circle, Karen to filled triangle
    cols <- col.lookup[as.character(map$BMI.Class)] 

    if(saveplot)
        pdf(file=outputfn,useDingbats=F)
        
   plot(map$Years.in.US, rel.distance, xlab="Years in US", ylab=ylab, main=main, pch=final.pch, col=alpha(cols,.4))
   m <- lm(rel.distance~map$Years.in.US)
    abline(m)
    
    # plot Rsquared and p value
    s <- summary(m)
    rsquared <- s$r.squared
    pval <- coef(s)[2,4]
    legend(x="topright",legend=paste("P =",round(pval,5), "\nR2 =",round(rsquared,4)), bty ="n", pch=NA, cex=.7) 
    if(saveplot)
        dev.off()
}