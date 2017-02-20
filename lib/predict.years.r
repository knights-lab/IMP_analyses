source("lib/rf.cross.validation.r")

predict.years <- function(map, otu)
{
     map[map$Sample.Group %in% c("KarenThai","HmongThai"),"Years.in.US"] <- 0

    # let's turn our groups into factors for classification
#     firstgen <- rownames(map)[map$Sample.Group %in% c("Hmong1st","Karen1st")]
#     firstgen.years <- cut(map[firstgen,"Years.in.US"], breaks=8)
#     names(firstgen.years) <- firstgen
#     
#     secondgen <- rownames(map)[map$Sample.Group=="Hmong2nd"]
#     secondgen.years <- rep("US-born",length(secondgen))
#     names(secondgen.years) <- secondgen
#     
#     thai <- rownames(map)[map$Sample.Group %in% c("KarenThai","HmongThai")]
#     thai.years <- rep("Thailand",length(thai))
#     names(thai.years) <- thai
# 
#     residence <- factor(c(firstgen.years, secondgen.years, thai.years), levels=c("Thailand", as.character(1:8), "US-born")) 
#     levels(residence)[2:9] <- c("0-5","5-10","10-15","15-20","20-25","25-30","30-35","35-40")
# 
#     map[,"residence"] <- residence[rownames(map)]

    otu <- otu[rownames(map),]
    
    res <- rf.cross.validation(x=otu, y=map$Years.in.US, nfolds=-1, regression=T)  
    
    plot(res$predicted, res$y, col=alpha("black",.5))

    
#     for(i in 1:n)
#     {
#  #       res_c[[i]] <- rf.cross.validation(x=otu, y=map$residence, nfolds=-1, regression=F)  
#         res[[i]] <- rf.cross.validation(x=otu, y=map$Years.in.US, nfolds=-1, regression=T)  
#     }
# 
#     rmse <- lapply(res, "[[", "rmse")
#     
#     saveRDS(res, file="rf_results_n50.rds")
# 
#     # read back using: readRDS(file="rf_results_n50.rds")
# 
#     # plot all 50 results
#     y <- map$Years.in.US
#     plot(y,predicted[[1]],main="predicted vs actual, n=50", xlab="Actual Years in US", ylab="Predicted Years in US", xlim=c(0,40),ylim=c(0,40) )
#     temp <- lapply(predicted, function(x) points(y,x))
#     
# 
#     predicted <- lapply(res, "[[", "predicted")
#     
#     clookup <- c("green","orange","red") 
#     names(clookup) <- sort(unique(map$BMI.Class)) # let's Hmong to solid filled circle, Karen to filled triangle
#     cols <- clookup[as.character(map$BMI.Class)] 
#     names(cols) <- rownames(map)
# 
#     lookup <- c(21,24) # point type
#     names(lookup) <- sort(unique(map$Ethnicity)) # let's Hmong to solid filled circle, Karen to filled triangle
#     pch2 <- lookup[as.character(map$Ethnicity)] 
#     names(pch2) <- rownames(map)
# 
# 
#     plot(predicted[[1]], map$Years.in.US, col=alpha(cols,.6), ylim=c(0,10), xlim=c(0,10), pch=pch2)

    return(res)
#     print(paste0("RMSE: ", mean(rmse)))
#     print(paste0("R-Squared: ", mean(rsq)))
}