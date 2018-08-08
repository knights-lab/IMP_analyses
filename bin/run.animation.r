# iMovie Notes
# 1. create new blank movie# 
# 2. import all files and drag to timeline
# 3. select all files then Window/Show Adjustments Bar -- under Crop, select "Fit"
# 4. with all files selected, click "i" icon and set duration to .1
# 5. all files selected, click Window/Movie Properties/Settings and adjust times to minimum possible
# 6. File/Share/File --> save as mp4
# 7. Create new movie, import mp4 from step 6, and drag to timeline
# 8. Select all, then Show Adjustments Bar - click on turtle icon and speed up movie.
# 9. File/Share/File --> save as mp4

library('vegan')
library('car')
library("plyr")
source('bin/load.r')
source('lib/animate.R')
dir.create("images")

# save these sample names for plotting background guys
non_longitudinal <- rownames(map)[map$Sub.Study == "CS"]
map0 <- map
dm <- bc_dm

## super important to include the entire DM here!
ddm <- as.dist(dm)
pc <- cmdscale(ddm,2)

pc <- -1*pc

col.lookup <- c(brewer.pal(8,"Set3"))
# drop yellow - hard to see!
col.lookup <- col.lookup[-2]

#map_L <- map[map$Subject.ID == "IMP.000",]
map_L <- map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",]

names(col.lookup) <- sort(as.character(unique(map_L[map_L$Sample.Order==max(map_L$Sample.Order), "Subject.ID"])))
cols <- col.lookup[as.character(map_L$Subject.ID)]

m1 <- sort(as.character(unique(map_L[map_L$Sample.Order==max(map_L$Sample.Order), "Subject.ID"])))

m1.samples <- rownames(map_L)[map_L$Sample.Order == 1]
m1.cols <- col.lookup[as.character(map_L[m1.samples, "Subject.ID"])]

for(step in seq(0,max(map_L$Sample.Order),.025)){
    cat(step,'\n')
   png(sprintf('images/step-%0.5f.png',step),width=800,height=800)
    par(mar=c(5.1, 4.1, 4.1, 9.1))
    
    plot(pc, type="n", xlab="PC1", ylab="PC2")

    # plot everyone else too 
    points(pc[non_longitudinal,1], pc[non_longitudinal, 2], 
         col=alpha("black", .2),
         pch=19, cex=1)

    # plot start points
    points(pc[m1.samples,1], pc[m1.samples,2], 
         col=m1.cols,
         pch=19,
         cex=1.5)
    
    # plot longitudinal samples
    plot.animation.frame(pc[rownames(map_L),], map_L$Subject.ID, map_L$Sample.Order, 
                         step=step, cols,
                         ellipse.alpha=.1, alpha.decay= 0.4) ## Modified Colors and Alpha Decay Level
        
    dev.off()
}

#source("code/cloud.sizes.r")

