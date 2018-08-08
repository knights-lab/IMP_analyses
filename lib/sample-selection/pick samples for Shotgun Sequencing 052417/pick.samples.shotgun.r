# pick samples for deep shotgun

# load all source files, data, and prep any variables for reuse
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/bin/load.r")

# exclude these, did NOT pass QC
# CS.171
# CS.082
# CS.083
# CS.010

# only CS.171 is an issue for our purposes 
#map[c("CS.171","CS.082","CS.083","CS.010"),]

dm <- wuf_dm[cs,cs] # best to use unifrac here
ddm <- as.dist(dm)
pc <- cmdscale(ddm,2)

map0 <- map[cs,]

bd <- betadisper(ddm, map0$Sample.Group)
centroid.dist <- bd$distances # distances of all samples to their group centroids!

hmong1 <- names(sort(centroid.dist[rownames(map0)[map0$Sub.Study=="CS" & map0$Years.in.US > 29 & map0$Sample.Group=="Hmong1st"]]))

bp <- plot.b.p.ratio(map0[cs,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="temp.b.p.ratio.pdf")
### comment out these lines because we don't want to pick high BP ratios anymore (try to be unbiased!)
# pick 1st gen samples if > 30x more bacteroides than prevotella, and those who have been in the US > 29 years (30 cut off gives us only n=14)
#hmong1 <- intersect(names(bp[bp > 30]), rownames(map0[map0$Sample.Group=="Hmong1st" & map0$Sub.Study=="CS" & map0$Years.in.US > 29,]))
# let's manually drop one sample (CS.284 has diabetes and is overweight, which we have too many of)
#hmong1 <- hmong1[-which(hmong1=="CS.284")]

hmongthai <- names(sort(centroid.dist[rownames(map0)[map0$Sample.Group=="HmongThai"]])[1:15])

d <- data.frame(x = pc[,1], y = pc[,2], group=map0$Sample.Group, distance=substring(as.character(centroid.dist[rownames(map0)]),1,4))
# set the levels of Sample.Group so that it's the same every time
d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))
group.cols <- get.group.colors(groups=as.character(levels(d$group)), alpha.val=.8)
p <- ggplot(data=d, aes(x, y)) + geom_point(colour=alpha("gray",.5), size=2) +
    scale_color_manual(values=group.cols) + #sets the color palette of the fill
    stat_ellipse(data=d, aes(colour=group), show.legend=T, type="t", level=.6)

# label points by distances
#p <- p + geom_text(data = d[hmong1, ], aes(label=distance) ,hjust=0, vjust=0, size=3)
# label points by sample name
p <- p + geom_text(data = d[hmong1, ], aes(label=hmong1) ,hjust=0, vjust=0, size=3)

# plot top 15 selected for each group
p <- p + geom_point(data=d[hmong1,], colour="black")
p <- p + geom_point(data=d[hmongthai,], colour="red")

# plot BP of the samples we picked
bp2 <- plot.b.p.ratio(map0[hmong1,], taxa, bug1=bacteroides, bug2=prevotella, outputfn="hmong1.b.p.ratio.pdf")

o.hmong.thai <- map0[hmongthai,c("Age","BMI")]
o.hmong.1 <- map0[hmong1,c("Age","BMI")]
o.hmong.thai <- o.hmong.thai[order(o.hmong.thai$Age,o.hmong.thai$BMI),]
o.hmong.1 <- o.hmong.1[order(o.hmong.1$Age,o.hmong.1$BMI),]

o.control <- map_all[map_all$Sample.Group=="Control", c("Age","BMI")]
o.control <- o.control[order(o.control$Age,o.control$BMI),]

# manually pick Control samples based on Age and BMI of Hmong1 and HmongThai
selected <- c("CS.177", "CS.146", "CS.369", "CS.393", "CS.262", "CS.181", "CS.391", "CS.268", "CS.166","CS.256","CS.326","CS.178","CS.255","CS.145","CS.370")

# compare selected samples across all three groups for Age and BMI
selected.control <- o.control[selected,]
plot(selected.control,type="n")
lines(selected.control[order(selected.control$Age),])
lines(o.hmong.1,col="blue")
lines(o.hmong.thai,col="red")


# now check that all samples are in the projects they should be in 
umgc <- read.table("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/sequences_052417/pick samples for Shotgun Sequencing/umgc.samples.txt", sep="\t", head=T)
merge(o.hmong.1, umgc, by.x=0, by.y="Sample.ID")
merge(o.hmong.thai, umgc, by.x=0, by.y="Sample.ID")

## NOTE, that some changes were made to the final list based on low DNA content in QC results -- 10/16/17