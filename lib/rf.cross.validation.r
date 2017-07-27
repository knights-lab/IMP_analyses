# runs cross-validation 
# if nfolds > length(y) or nfolds==-1, uses leave-one-out cross-validation
# ...: additional parameters for train.fun
#
# value:
# y: true values
# predicted: cv predicted values
# probabilities: cv predicted class probabilities (or NULL if unavailable)
# confusion.matrix: confusion matrix (true x predicted)
# nfolds: nfolds, use -1 for leave-one-out cv
# params: list of additional parameters
# importances: importances of features as predictors
# regression: does both regression or classification (RMSE and r_squared will be meaningless here)
require(randomForest)
"rf.cross.validation" <- function(x, y, nfolds=10, verbose=FALSE, regression=FALSE, ...){
    if(regression==FALSE)
   	{
		if(class(y) != 'factor') stop('y must be factor for classification\n')
		y <- droplevels(y)
   	}
   
    if(nfolds==-1) nfolds <- length(y)
    #folds <- sample(rep(1:nfolds,ceiling(length(y)/nfolds)))
	folds <- balanced.folds(y, nfolds)
    
    result <- list()
    result$y <- y
    result$predicted <- result$y

    # K-fold cross-validation
    for(fold in sort(unique(folds))){
        if(verbose) cat(sprintf('Fold %d...\n',fold))
        foldix <- which(folds==fold)
        model <- randomForest(x[-foldix,], result$y[-foldix], importance=TRUE, do.trace=verbose, ...)
        newx <- x[foldix,]
        if(length(foldix)==1) newx <- matrix(newx,nrow=1)
        result$predicted[foldix] <- predict(model, newx)
    }

    result$rmse <- sqrt(mean((y - result$predicted)**2)) #standard error of the estimate 
    # to get 95% confidence, multiply rmse by 1.96
	SS_res <- sum((y - result$predicted)**2)
	SS_tot <- sum((y - mean(y))**2)
    result$r_squared <- 1-(SS_res/SS_tot)
	result$nfolds <- nfolds
    result$params <- list(...)
    return(result)    
}

