source(paste0(LIBDIR,"rf.cross.validation.r"))
source(paste0(LIBDIR,"balanced.folds.r"))

predict.run <- function()
{
    maphmong <- map[hmong_firstgen_cs,]
    mapkaren <- map[karen_firstgen_cs,]
    otu0 <- taxa_L7
    # save workspace
    # load workspace on teraminx, then
    # run these next lines on teraminx
    p.hmong <- predict.years(map0=maphmong, otu0=otu0, n=99)
    p.karen <- predict.years(map0=mapkaren, otu0=otu0, n=99)
}

predict.years <- function(map0, otu0, n=99)
{
    food_wuf_dm <- read.table(paste(food_datadir, "wuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
    food_uwuf_dm <- read.table(paste(food_datadir,"uwuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)


    firstgen <- c(karen_firstgen_cs, hmong_firstgen_cs)
    write.table(firstgen, file="firstgen.txt", quote=F, row.names=F, col.names=F)
    write.table(map, file="map.txt", quote=F, sep="\t")
    food.pcoa <- plot.pcoa(map[firstgen,], dm=food_wuf_dm, plot.title="food.pcoa", show.stats=F, save.pc=T, axis2=10)
    mb.pcoa <- plot.pcoa(map[firstgen,], dm=wuf_dm, plot.title="mb.pcoa", show.stats=F, save.pc=T, axis2=10)

    firstgen <- read.table("firstgen.txt", sep="\t", head=F, stringsAsFactors=F)[,1]
    food.pc <- read.table("food.pcoa-PC.txt", sep="\t", head=T, row=1)
    mb.pc <- read.table("mb.pcoa-PC.txt", sep="\t", head=T, row=1)
    map <- read.table("map.txt", sep="\t", head=T, row=1, stringsAsFactors=F)

    food.pc <- food.pc[firstgen,]
    colnames(food.pc) <- paste0("Food.PC",1:10)
    mb.pc <- mb.pc[firstgen,] 
    colnames(mb.pc) <- paste0("MB.PC",1:10)
    
    y <- map[firstgen,"Years.in.US"]
    
    rf.food <- rf.cross.validation(x=food.pc, y=y, nfolds=-1, regression=T)  
    rf.food$rmse
    
    rf.mb <- rf.cross.validation(x=mb.pc, y=y, nfolds=-1, regression=T)  
    rf.mb$rmse
    
    rf.food.mb <- rf.cross.validation(x=cbind(mb.pc,food.pc), y=y, nfolds=-1, regression=T)  
    rf.food.mb$rmse
    
   #  res.shuffled <- list()
    #     rmse <- NULL
    #     for(i in 1:n)
    #     {
    #         shuffled.y <- sample(y)
    #         res.shuffled[[i]] <- rf.cross.validation(x=otu0, y=shuffled.y, nfolds=-1, regression=T)  
    #         rmse[i] <- res.shuffled[[i]]$rmse
    # #        print(res.shuffled[[i]]$rmse)
    #     }
    #     return(list(actual=res.actual$rmse, shuffled=rmse))
    #return(mean(res.actual$rmse < rmse))
}
