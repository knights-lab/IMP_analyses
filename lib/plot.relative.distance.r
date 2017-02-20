require(faraway)

# plot relative distances of longitudinal samples to either a set of ref samples OR a specific reference timepoint
plot.relative.L <- function(map, dm, ref.sample.order, xlab, ylab, outputfn, ref_samples=NULL, saveplot=T)
{
    cols=c(brewer.pal(7,"Set3"))
    
    r <- range(map$Sample.Order)
    timepoints <- r[1]:r[2]
    
    all.rel.distance <- NULL
    subjects <- unique(map[map$Sample.Order==max(map$Sample.Order), "Subject.ID"])
    
    if(saveplot)
        pdf(file=outputfn,useDingbats=F)

    for(i in 1:length(subjects))
    {
        sub_map <- map[map$Subject.ID==subjects[i],]
        sub_map <- sub_map[order(sub_map$Sample.Order),] 
        
        query_samples <- rownames(sub_map)
        
        if(length(ref_samples)==0)
            ref_samples <- rownames(sub_map)[sub_map$Sample.Order == ref.sample.order]
    
        rel.distance <- NULL
        for (j in 1:length(query_samples))
        {
            rel.distance[j] <- mean(as.numeric(dm[query_samples[j], ref_samples]))
        }
        names(rel.distance) <- query_samples
        all.rel.distance <- c(all.rel.distance, rel.distance)
        
        if(i==1)
            plot(sub_map[query_samples,"Sample.Order"], rel.distance[query_samples], ylim=c(0,1), 
            main="Relative Distance to Self", col=cols[i], xlab=xlab, ylab=ylab, lwd=2, type="l")
        else
            lines(sub_map[query_samples,"Sample.Order"], rel.distance[query_samples], col=cols[i], lwd=2) 
    }
    if(saveplot)
        dev.off()
}

# plot relative distance of cross sectional samples to a set of reference samples
plot.relative.CS <- function(map, dm, main, query_samples, ref_samples, outputfn, ylab, saveplot=T)
{
    map0 <- map
    map <- map[query_samples,]

    rel.distance.ref <- NULL
    for (i in 1:length(ref_samples))
    {
        rel.distance.ref[i] <- mean(as.numeric(dm[ref_samples[i], ref_samples]))
    }


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

    cols <- "black"
#     col.lookup <- c("#1a9641","#fecc5c","#d7191c")
#     names(col.lookup) <- sort(unique(map$BMI.Class)) # let's Hmong to solid filled circle, Karen to filled triangle
#     cols <- col.lookup[as.character(map$BMI.Class)] 

    if(saveplot)
        pdf(file=outputfn,useDingbats=F)
        
    plot(map$Years.in.US, rel.distance, xlab="Years in US", ylab=ylab, main=main, pch=final.pch, col=alpha(cols,.4))
    m <- lm(rel.distance~map$Years.in.US)
    abline(m)
    
#    cols2 <- col.lookup[as.character(map0[ref_samples,"BMI.Class"])] 
#    points(rep(42,length(rel.distance.ref)), rel.distance.ref, col=alpha(cols2,.8), pch=19)

#     legend(pch=c(19,17), col=c("black","black","white","white"), x="topright", cex=.8, legend=c("Hmong", "Karen", "", ""))
#     # add gradient image for age colors
#     legend_image <- as.raster(rbPal(10)[10:1])
#     text(y=.553, x = c(14.2,15,15.8), labels = c("17","BMI","39"), cex=c(.6,.8,.6))
#     rasterImage(legend_image, 14.1, .55, 14.3, .569, angle=270)

    # plot Rsquared and p value
    s <- summary(m)
    rsquared <- s$r.squared
    pval <- coef(s)[2,4]
#    text(x=max(map$Years.in.US)-3, y=max(rel.distance)-c(.01, .02), labels=c(paste("P =",round(pval,5)), paste("R2 =",round(rsquared,4))), cex=.7)
    legend(x="topright",legend=paste("P =",round(pval,5), "\nR2 =",round(rsquared,4)), bty ="n", pch=NA, cex=.7) 
    if(saveplot)
        dev.off()
}