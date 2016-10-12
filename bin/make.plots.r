setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.alphadiv.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.relative.distance.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/fraction.hits.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/load.data.r")

datadir <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/rare10002"

# load data
wuf_dm <- read.table(paste(datadir,"weighted_unifrac_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
uwuf_dm <- read.table(paste(datadir,"unweighted_unifrac_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
bc_dm <- read.table(paste(datadir,"bray_curtis_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
# let's try uwuf first
dm <- uwuf_dm

alphadiv0 <- read.table(paste(datadir,"alpha_div.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")

# load mapping and otu/taxa file - automatically normalizes, filters, and merges 
ret <- load.data(mapfile="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mapping.txt",
        otufile="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/taxon_summaries/taxa_s2_L6.txt")
map <- ret$map
taxa <- ret$otu

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
plot.fraction.hits()

# plot distance of Asian-born to US-born, over first 20 years in US
plot.relative.distance(dm, asian_born_20, US_born, filename="relative_distance_20yrs.pdf")

# plot alpha diversity
plot.alphadiv(cbind(map[valid_samples,], alphadiv0[valid_samples,"PD_whole_tree"], stringsAsFactors=F))

# plot intra-inter group variabilities
plot.group.distances(map, dm, between.groups.fn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/between.groups.txt")

#
plot.b.p.ratio(map, taxa)

