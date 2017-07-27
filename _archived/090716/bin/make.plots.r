setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.alphadiv.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.relative.distance.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/fraction.hits.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/load.data.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.b.p.ratio.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.group.distances.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.nutrients.r")

datadir_99rare <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/99/rare10002"
datadir_97 <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/97"
datadir <- datadir_97

otufn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/97/taxon_summaries/taxa_s2_L6.txt"


# load data
wuf_dm <- read.table(paste(datadir, "weighted_unifrac_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
uwuf_dm <- read.table(paste(datadir,"unweighted_unifrac_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
bc_dm <- read.table(paste(datadir,"bray_curtis_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
# let's try uwuf first
dm <- uwuf_dm

alphadiv0 <- read.table(paste(datadir,"alpha.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")

# load mapping and otu/taxa file - automatically normalizes, filters, and merges 
ret <- load.data(mapfile="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mapping.txt", otufile=otufn)
map <- ret$map
taxa <- ret$otu

# write out relative abundance OTU file out 

# recode metadata
map$Age <- as.numeric(map$Age)
map$Years.in.US <- as.numeric(map$Years.in.US)
# originally NA
map$Years.in.US[is.na(map$Years.in.US)] <- 0
# refactor BMI.Class to order by BMI instead of alphabetical
map$BMI.Class <- factor(map$BMI.Class, levels=c("Normal", "Overweight", "Obese")) 
# let's interpret Karenni as Karen - they are essentially the same (n=2)
map$Ethnicity[map$Ethnicity=="Karenni"] <- "Karen"
# let's remap the residence years for plotting later
Residence.Class <- character(nrow(map))
Residence.Class[map$Years.in.US != 0 & map$Years.in.US < 5] <- "Short"
Residence.Class[map$Years.in.US >=5 & map$Years.in.US < 10] <- "Medium"
Residence.Class[map$Years.in.US >=10 & map$Years.in.US < 15] <- "Long"
Residence.Class[map$Years.in.US >=15 ] <- "VeryLong"
Residence.Class[map$Years.in.US == 0 ] <- "USborn"    
Residence.Class <- factor(Residence.Class, levels=c("Short", "Medium", "Long", "VeryLong", "USborn")) 
map <- data.frame(map, Residence.Class)
# add waist-to-heigh ratio
map$Waist.Height.Ratio <- map$Waist/map$Height

write.table(map, "newmap.txt",sep="\t",quote=F)

# only work with samples that are in both the mapping AND the dm/alpha files
valid_samples <- intersect(rownames(dm), rownames(map))
map <- map[valid_samples,]
dm <- dm[valid_samples, valid_samples]

# subset samples
asian_born <- rownames(map)[map$Ethnicity %in% c("Hmong", "Karen") & map$Years.in.US != 0]
asian_born_20 <- rownames(map)[map$Ethnicity %in% c("Hmong", "Karen") & map$Years.in.US != 0 & map$Years.in.US < 20] 
US_born <- rownames(map)[map$Ethnicity %in% c("Hmong") & map$Years.in.US == 0] 

###### make plots #########

# plot fraction of sequences hitting GG
# fractionfn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/fraction_hitting_GG99.txt"
fractionfn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/fraction_hitting_GG97.txt"
plot.fraction.hits(fractionfn, map)

# plot distance of Asian-born to US-born, over first 20 years in US
plot.relative.distance(dm, asian_born_20, US_born, filename="relative_distance_20yrs.pdf")

plot.relative.distance(dm, asian_born, US_born, filename="relative_distance.pdf")

# plot alpha diversity
library(beeswarm)
plot.alphadiv(cbind(map[valid_samples,], alphadiv0[valid_samples,"PD_whole_tree"], stringsAsFactors=F))

# plot intra-inter group variabilities
plot.group.distances(map, dm, between.groups.fn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/between.groups.txt")
# plot for ISU only
plot.group.distances(map, dm, between.groups.fn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/between.groups.ISU.txt")

# bacteroides to prevotella
plot.b.p.ratio(map, taxa)

# plot nutrients 
# try limiting to first 10 years

plot.nutrients(map0=map[map$Years.in.US > 0 & map$Years.in.US < 15,], ethnicity=c("Hmong","Karen"), nutrient.vars=c("Added.Sugars.in.Grams",
    "Total.Calories", "Protein.in.Grams"),
    nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")
plot.nutrients(map0=map[map$Years.in.US > 0 & map$Years.in.US < 15,], ethnicity=c("Hmong","Karen"),
    nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")

plot.nutrients(map0=map[map$Years.in.US > 0,], ethnicity="Hmong", nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")
plot.nutrients(map0=map[map$Years.in.US > 0 & map$Years.in.US <=10,], ethnicity="Hmong", nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")
plot.nutrients(map0=map, ethnicity="Karen", nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")

plot.body.trends(map)

require(vegan)
plot.pcoa <- function()
{
    maporig <- map

    map <- map[map$Years.in.US > 0,]
    valid_samples <- intersect(rownames(dm), rownames(map))
    map <- map[valid_samples,]
    dm <- dm[valid_samples, valid_samples]

    ddm <- as.dist(dm)
    # plot clusters
    # get PCoA scores first
    pc <- cmdscale(ddm,2)
    
    rbPal <- colorRampPalette(alpha(c('wheat', "red")))
    cols2 <- rbPal(10)[as.numeric(cut(map$Years.in.US, breaks = 10))]
    names(cols2) <- rownames(map)
    
    plot(pc[,1], pc[,2], col=cols2, xlab="PC1",ylab="PC2",pch=19, main="Unweighted Unifrac")


}

