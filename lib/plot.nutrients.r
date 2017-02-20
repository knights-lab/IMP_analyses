# load nutrients and merge with mapping
# map is a dataframe that has already been filtered for samples desired
# consider using Principal Components to break down nutrients instead of hand-picking?
# consider looking at % of calories?
require(faraway)
plot.nutrients <- function(map0, ethnicity, nutrientsfn, nutrient.vars=NULL)
{
    nutrients0 <- read.table(nutrientsfn, sep="\t", header=T, row=1)

    samples <- intersect(rownames(nutrients0), map0$Old.Participant.ID.x)
    map <- map0[map$Old.Participant.ID.x %in% samples,]
    nutrients <- nutrients0[samples,]
    
    # subset map 
    ids <- rownames(map[map$Ethnicity %in% ethnicity,])
    map <- map[ids,]
    nutrients <- nutrients[ids,]
    
    lookup <- c("green","orange","red")
    names(lookup) <- levels(map$BMI.Class) 
    cols <- lookup[as.character(map$BMI.Class)] 
    names(cols) <- rownames(map)

    plookup <- c(19,17) # point type
    names(plookup) <- sort(unique(map$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
    pch2 <- plookup[as.character(map$Ethnicity)] 
    names(pch2) <- rownames(map)


    #print(table(map[,c("Ethnicity","BMI.Class")]))
    
    if(length(nutrient.vars)==0) # if not set, set it to all the columns in nutrients table
    {
        nutrient.vars <- colnames(nutrients)[c(1,3,4,5,6,7)]
    }
    
    cov1 <- "Years.in.US" # "Waist.Height.Ratio" #"BMI" "Years.in.US"
    fn <- paste("Nutrient",ethnicity,cov1,sep="_")
#    cat("Linear model interactions and main effects\n\n",file=paste(fn,".txt",sep=""), append=F)
#    pdf(file=paste(fn,".pdf",sep=""),useDingbats=F)
    par(mfrow=c(3,3))
    for(i in 1:length(nutrient.vars))
    {
        nutrient <- nutrient.vars[i]
        f <- nutrients[,nutrient] ~ map$Years.in.US * map$Age * map$BMI * map$Ethnicity #map$Waist.Height.Ratio
        # remove 3-way interactions
        #f <- nutrients[,nutrient] ~ map$Years.in.US + map$Age + map$BMI + map$Years.in.US:map$Age + map$Years.in.US:map$BMI + map$Age:map$BMI
        m <- lm(f)
        s <- summary(m)
        rsquared <- s$r.squared

        

        # plot half normal plots
#         coeffs <- coef(m)
#         halfnorm(coeffs[-1], nlab=3, labs=names(coeffs[-1]), main=nutrient)
#         cat(nutrient,file=paste(fn,".txt",sep=""), append=T)
#         write.table(coef(s),file=paste(fn,".txt",sep=""), append=T, sep="\t", quote=F)
        
        # plot scatter plots of associations
        pval <- coef(s)[paste("map$",cov1,sep=""),4]
        plot(nutrients[,nutrient] ~ map[,cov1], col=alpha(cols,.5), ylab=nutrient, xlab=cov1, pch=pch2,
                main=nutrient)
        abline(lm(nutrients[,nutrient] ~ map[,cov1]), col="gray")
        mtext(text=paste("P =",round(pval,4),"R2 =",round(rsquared,4)), side=3, cex=.7)
        # legend for special plotting purposes only
        # legend(pch=c(19), col=alpha(lookup,.5), x="topright", cex=.8, legend=c("Lean", "Overweight", "Obese"))


    }
#    dev.off()
}

