source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/alphadiv.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/boxplot.by.group.x.bmi.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/relative.distance.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/fraction.hits.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/load.data.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/b.p.ratio.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/group.distances.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/nutrients.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/pcoa.color.by.nutrients.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/alphadiv.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/pcoa.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/body.trends.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/predict.years.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/utils.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/nutrient.mb.heatmap.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/taxa.summary.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/differential.taxa.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/heatmap.r")


datadir_gg97 <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/gg97"
datadir_dada2 <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/dada2"
datadir_refseq <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/refseq"

datadir <- datadir_refseq

otu_fn <- paste(datadir,"final_otu.txt",sep="/") # embalmer generated OTU file does not have taxonomy included
taxa_L6_fn <- paste(datadir,"taxatable_L6.txt",sep="/")
taxa_L2_fn <- paste(datadir,"taxatable_L2.txt",sep="/")
taxafile_L6_rare_fn <- paste(datadir,"taxatable_rare_L6.txt",sep="/")
taxafile <- taxa_L6_fn 

wuf_dm_fn <- "weighted_unifrac_dm.txt"
uwuf_dm_fn <- "unweighted_unifrac_dm.txt"
bc_dm_fn <- "bray_curtis_dm.txt"
alpha_fn <- "alpha.txt"


if(use.rare){ # consider just loading everything
    wuf_dm_fn <- "weighted_unifrac_dm_rare.txt"
    uwuf_dm_fn <- "unweighted_unifrac_dm_rare.txt"
    bc_dm_fn <- "bray_curtis_dm_rare.txt"
    alpha_fn <- "alpha_rare.txt"
}

#mapfile <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mapping.txt"
# TEMPORARILY use the mapping file with TFSCS for sample names for KCK
# (these samples were chipped off for sequencing for FMT only and renamed by UMGC - they'll be resequenced later with the correct names)
orig_mapfile <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mapping.txt"
mapfile <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mapping_with_KCKFMT_renamed.txt"
nutrientsfn<-"/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/nutrients.txt"
foodgroupfn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/foodgroups.txt"

bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"

# load data
wuf_dm <- read.table(paste(datadir, wuf_dm_fn,sep="/"), sep="\t", quote="", row=1, head=T)
uwuf_dm <- read.table(paste(datadir, uwuf_dm_fn,sep="/"), sep="\t", quote="", row=1, head=T)
bc_dm <- read.table(paste(datadir, bc_dm_fn,sep="/"), sep="\t", quote="", row=1, head=T)

alphadiv0 <- read.table(paste(datadir, alpha_fn, sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")

# load food data files
food_wuf_dm <- read.table(paste(datadir, "../wuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
food_uwuf_dm <- read.table(paste(datadir,"../uwuf_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
food_euc_dm <- read.table(paste(datadir,"../euc_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
food_bc_dm <- read.table(paste(datadir,"../bc_food_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
food_otu_L3 <- read.table(paste(datadir,"../food.otu_L3.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment="", skip=1)
food_alpha <- read.table(paste(datadir,"../food.alpha.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")

t_food_otu <- t(food_otu_L3)
food_otu_L3 <- sweep(t_food_otu, 1, rowSums(t_food_otu), '/')

# load mapping and otu/taxa file - automatically normalizes, filters, and merges 
# let's turn off normalization efforts. The taxa file should already be in relative abundance.
ret <- load.data(mapfile, otufile=taxafile, normalize=F)
map <- ret$map
taxa <- ret$otu

taxa_L2 <- load.data(mapfile, otufile=taxa_L2_fn, normalize=F)$otu
taxa_L6_rare <- load.data(mapfile, otufile=taxafile_L6_rare_fn, normalize=F)$otu

otu <- read.table(otu_fn, sep="\t", quote="", row=1, head=T, comment="")
otu <- as.data.frame(t(otu))

# only work with samples that are in both the mapping AND the dm/alpha files
valid_samples <- intersect(rownames(bc_dm), rownames(map))
map <- map[valid_samples,]

# format important variables factor levels and/or remove levels (change to char)
map$BMI.Class <- as.character(map$BMI.Class)
map$BMI.Class[map$BMI.Class == "Normal"] <- "Lean" # replace "Normal" with "Lean"
map$BMI.Class <- factor(map$BMI.Class, levels=c("Lean", "Overweight", "Obese")) 
map$Subject.ID <- as.character(map$Subject.ID)
map$Sample.Group <- factor(map$Sample.Group, levels=c("KarenThai","HmongThai","Karen1st","Hmong1st","Hmong2nd","Control"))

# all single-timepoint samples across both countries
cs <- rownames(map)[is.na(map$Sample.Order) | map$Sample.Order==1]
# 1st generation single-timepoint only
firstgen_cs <- rownames(map)[map$Sample.Group %in% c("Hmong1st","Karen1st") & (is.na(map$Sample.Order) | map$Sample.Order==1)]

hmong_secondgen_cs <- rownames(map)[map$Sample.Group == "Hmong2nd"]
hmong_firstgen_cs <- rownames(map)[map$Sample.Group == "Hmong1st" & map$Subject.ID != "IMP.000"]
karen_firstgen_cs <- rownames(map)[map$Sample.Group == "Karen1st" & (is.na(map$Sample.Order) | map$Sample.Order==1)]
karenthai <- rownames(map)[map$Sample.Group=="KarenThai"]
hmongthai <- rownames(map)[map$Sample.Group=="HmongThai"]

# calculate % of life spent in the US column
map[,"Fraction.Life.in.US"] <- map$Years.in.US/map$Age
map[hmong_secondgen_cs,"Fraction.Life.in.US"] <- 1.0
map[c(karenthai,hmongthai),"Fraction.Life.in.US"] <- 0

# let's reset Years.in.US so that 2ndGen == 50 and Thai == 0
map[map$Years.in.US==0 & !is.na(map$Years.in.US),"Years.in.US"] <- 50
map[is.na(map$Years.in.US),"Years.in.US"] <- 0

# calculate a socioeconomics factor columns
map[,"Socioeconomic.Index"] <- mean(map[,c("Children.Free.Lunch","Medical.Assistance","Public.Housing")] == "N", na.rm=T)

dm <- bc_dm[cs,cs]
ddm <- as.dist(dm)
pc <- cmdscale(ddm,2)
# ** important here that subsetting needs to be done after PCs are generated

# lets reload mapping file containing ALL participants (even those not sequenced yet)
# note that all food-related data has already been indexed by SAMPLE ID
map_all <- read.table(orig_mapfile, sep="\t", header=T, row=1, check.names=F, comment.char="", as.is=T)
# always remove excluded participants
map_all <- map_all[map_all$Exclude!="Y",]
# remove any FMT samples (no need to analyze them)
map_all <- map_all[map_all$Sub.Study!="FMT",]
# remove IMP000 samples from food (non exist)
map_all <- map_all[map_all$Subject.ID != "IMP.000",]
cs_all <- rownames(map_all)[is.na(map_all$Sample.Order) | map_all$Sample.Order==1]
map_all$BMI.Class <- as.character(map_all$BMI.Class)
map_all$BMI.Class[map_all$BMI.Class == "Normal"] <- "Lean"
map_all$BMI.Class <- factor(map_all$BMI.Class, levels=c("Lean", "Overweight", "Obese")) 
map_all$Subject.ID <- as.character(map_all$Subject.ID)
map_all$Sample.Group <- factor(map_all$Sample.Group, levels=c("KarenThai","HmongThai","Karen1st","Hmong1st","Hmong2nd","Control"))

# _all includes all participants, including those that have NOT been sequenced. important for metadata correlations.
hmong_secondgen_cs_all <- rownames(map_all)[map_all$Sample.Group == "Hmong2nd"]
hmong_firstgen_cs_all <- rownames(map_all)[map_all$Sample.Group == "Hmong1st" & map_all$Subject.ID != "IMP.000"]
karen_firstgen_cs_all <- rownames(map_all)[map_all$Sample.Group == "Karen1st" & (is.na(map_all$Sample.Order) | map_all$Sample.Order==1)]
karenthai_all <- rownames(map_all)[map_all$Sample.Group=="KarenThai"]
hmongthai_all <- rownames(map_all)[map_all$Sample.Group=="HmongThai"]
controls_all <- rownames(map_all)[map_all$Sample.Group=="Control"]

map_all[map_all$Years.in.US==0 & !is.na(map_all$Years.in.US),"Years.in.US"] <- 50
map_all[is.na(map_all$Years.in.US),"Years.in.US"] <- 0


# nutrient file should already be mapped to sample ids
nutrients <- read.table(nutrientsfn, sep="\t", header=T, check.names=F, as.is=T,row=1)
nutrients[,"% of Calories from Total Sugars"] <- nutrients[,"Total Sugars in Grams"] * 4 / nutrients[,"Total Calories"]
nutrients[,"g Fiber per 1000 Calories"] <- (nutrients[,"Dietary Fiber in Grams",] / nutrients[,"Total Calories"]) * 1000

# let's add some diet metadata into the mapping file 
x <- merge(map_all, nutrients[,c("SuperTracker.DATE", "Diet.Month", "Diet.ID", "Diet.Date"), drop=F], by=0)
rownames(x) <- x[,1]
x <- x[,-1]
map_all <- x
 
foodgroups <- read.table(foodgroupfn, sep="\t", header=T, check.names=F, as.is=T)
fgroup <- paste0(foodgroups$FoodType, " in ", foodgroups$PortionUnit)
foodgroups <- data.frame(foodgroups[, c("Sample.ID","Amount")],fgroup)
foodgroups <- reshape(foodgroups,direction="wide", idvar="Sample.ID", timevar="fgroup")

rownames(foodgroups) <- foodgroups$Sample.ID
foodgroups <- foodgroups[,-which(colnames(foodgroups)=="Sample.ID")]

# calculate Days since arrival for sample and diet dates
map_all$Sample.Day.Since.Arrival <- as.numeric(as.Date(map_all$Sample.Date, format="%m/%d/%y") - as.Date(map_all$Arrival.in.US, format="%m/%d/%y"))
map_all$Diet.Day.Since.Arrival <- as.numeric(as.Date(map_all$Diet.Date, format="%m/%d/%y") - as.Date(map_all$Arrival.in.US, format="%m/%d/%y") )
map$Sample.Day.Since.Arrival <- as.numeric(as.Date(map$Sample.Date, format="%m/%d/%y") - as.Date(map$Arrival.in.US, format="%m/%d/%y"))




