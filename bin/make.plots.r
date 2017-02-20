setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.alphadiv.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.relative.distance.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/fraction.hits.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/load.data.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.b.p.ratio.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.group.distances.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.nutrients.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.alphadiv.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.pcoa.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/plot.body.trends.r")
source("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/predict.years.r")

datadir <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data"
otufn <- paste(datadir,"final_otu_L6.txt",sep="/")

# load data
wuf_dm <- read.table(paste(datadir, "weighted_unifrac_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
uwuf_dm <- read.table(paste(datadir,"unweighted_unifrac_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
bc_dm <- read.table(paste(datadir,"bray_curtis_dm.txt",sep="/"), sep="\t", quote="", row=1, head=T)
# let's try uwuf first
dm <- uwuf_dm
#dm <- wuf_dm

alphadiv0 <- read.table(paste(datadir,"alpha.txt",sep="/"), sep="\t", quote="", row=1, head=T, comment.char="")

# load mapping and otu/taxa file - automatically normalizes, filters, and merges 
# let's turn off normalization efforts. The taxa file should already be in relative abundance.
ret <- load.data(mapfile=paste(datadir,"mapping.txt",sep="/"), otufile=otufn, normalize=F)
map <- ret$map
taxa <- ret$otu

# write out relative abundance OTU file out 

# only work with samples that are in both the mapping AND the dm/alpha files
valid_samples <- intersect(rownames(dm), rownames(map))
map <- map[valid_samples,]
dm <- dm[valid_samples, valid_samples]

# required: refactor BMI.Class to order by BMI instead of alphabetical
map$BMI.Class <- factor(map$BMI.Class, levels=c("Normal", "Overweight", "Obese")) 

# include 1st sample of longitudinal samples in any CS analyses
cs_samples <- rownames(map)[is.na(map$Sample.Order) | map$Sample.Order==1]

################ make plots ###################

plot.alphadiv(cbind(map[cs_samples,], alphadiv0[cs_samples,"PD_whole_tree"], stringsAsFactors=F))

plot.body.trends(map[cs_samples,])

plot.pcoa(map[cs_samples,], wuf_dm, "Weighted Unifrac")
plot.pcoa(map[cs_samples,], uwuf_dm, "Unweighted Unifrac")
plot.pcoa(map[cs_samples,], bc_dm, "Bray Curtis")

# plot relative distance of all samples to US-born
non_second_gen_cs <- rownames(map)[map$Sample.Group != "Hmong2nd" & (is.na(map$Sample.Order) | map$Sample.Order==1)]
second_gen_cs <- rownames(map)[map$Sample.Group == "Hmong2nd"]
hmong_firstgen_cs <- rownames(map)[map$Sample.Group == "Hmong1st" & map$Subject.ID != "IMP.000"]
karen_firstgen_cs <- rownames(map)[map$Ethnicity=="Karen" & (is.na(map$Sample.Order) | map$Sample.Order==1)]
karenthai <- rownames(map)[map$Sample.Group=="KarenThai"]
hmongthai <- rownames(map)[map$Sample.Group=="HmongThai"]

plot.relative.CS(map=map, dm=uwuf_dm, query_samples=non_second_gen_cs, ref_samples=second_gen_cs, filename="relative_distance_uwuf.pdf", ylab="Distance to US-born Hmong")
plot.relative.CS(map=map, dm=wuf_dm, query_samples=non_second_gen_cs, ref_samples=second_gen_cs, filename="relative_distance_wuf.pdf", ylab="Distance to US-born Hmong")

    # try Karen only - distance to KarenThai
    pdf("Karen_CS_to_Thai.pdf",useDingbats=F, width=7, height=2.5)
    par(mfrow=c(1,3))
    ylab <- "Distance to Karen in Thailand"
    plot.relative.CS(map=map, dm=uwuf_dm, main="Unweighted Unifrac", query_samples=karen_firstgen_cs, ref_samples=karenthai, ylab=ylab,saveplot=F)
    plot.relative.CS(map=map, dm=wuf_dm, main="Weighted Unifrac", query_samples=karen_firstgen_cs, ref_samples=karenthai, ylab=ylab,saveplot=F)
    plot.relative.CS(map=map, dm=bc_dm, main="Bray Curtis", query_samples=karen_firstgen_cs, ref_samples=karenthai, ylab=ylab,saveplot=F)
    dev.off()
    
    # Hmong only
    pdf("Hmong_CS_to_Thai.pdf",useDingbats=F)
    par(mfrow=c(1,3))
    plot.relative.CS(map=map, dm=uwuf_dm, query_samples=hmong_firstgen_cs, ref_samples=hmongthai, filename="hmong_relative_distance_uwuf.pdf", ylab="Distance to Hmong in Thailand",saveplot=F)
    plot.relative.CS(map=map, dm=wuf_dm, query_samples=hmong_firstgen_cs, ref_samples=hmongthai, filename="hmong_relative_distance_wuf.pdf", ylab="Distance to Hmong in Thailand",saveplot=F)
    plot.relative.CS(map=map, dm=bc_dm, query_samples=hmong_firstgen_cs, ref_samples=hmongthai, filename="hmong_relative_distance_bc.pdf", ylab="Distance to Hmong in Thailand",saveplot=F)
    dev.off()

    # plot relative distances of longitudinal 
    pdf("Karen_L_to_M1.pdf",useDingbats=F)
    par(mfrow=c(1,3))
    plot.relative.L(map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], dm=uwuf_dm,
                            ref_samples.order = 1, xlab="Month in the US", ylab="Distance to Self at Month 1", outputfn="Long_to_M1_uwuf.pdf",saveplot=F)
    plot.relative.L(map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], dm=wuf_dm,
                            ref.sample.order = 1, xlab="Month in the US",  ylab="Distance to Self at Month 1", outputfn="Long_to_M1_wuf.pdf",saveplot=F)
    plot.relative.L(map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], dm=bc_dm,
                            ref.sample.order = 1, xlab="Month in the US",  ylab="Distance to Self at Month 1", outputfn="Long_to_M1_bc.pdf",saveplot=F)
    dev.off()
    
    # plot relative distances of longitudinal to Karen in Thailand
    pdf("Karen_L_to_Thai.pdf",useDingbats=F)
    par(mfrow=c(1,3))
    plot.relative.L(map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], dm=uwuf_dm,
                            ref.sample.order = 1, xlab="Month in the US", ylab="Distance to Karen in Thailand", outputfn="Long_to_Thai_uwuf.pdf",
                            ref_samples=karenthai,saveplot=F)
    plot.relative.L(map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], dm=wuf_dm,
                            ref.sample.order = 1, xlab="Month in the US", ylab="Distance to Karen in Thailand", outputfn="Long_to_Thai_wuf.pdf",
                            ref_samples=karenthai,saveplot=F)
    plot.relative.L(map = map[map$Sub.Study == "L" & map$Subject.ID != "IMP.000",], dm=bc_dm,
                            ref.sample.order = 1, xlab="Month in the US", ylab="Distance to Karen in Thailand", outputfn="Long_to_Thai_bc.pdf",
                            ref_samples=karenthai,saveplot=F)
    dev.off()

    # plot relative distance to self
    pdf("Self_to_Day1.pdf", useDingbats=F)
    plot.relative.L(map = map[map$Subject.ID == "IMP.000",], dm=uwuf_dm,
                            ref.sample.order = 1, xlab="Day", ylab="Distance to Self at Day 1", outputfn="Self_to_Day1_uwuf.pdf",saveplot=F)
    plot.relative.L(map = map[map$Subject.ID == "IMP.000",], dm=wuf_dm,
                            ref.sample.order = 1, xlab="Day",  ylab="Distance to Self at Day 1", outputfn="Self_to_Day1_wuf.pdf",saveplot=F)
    plot.relative.L(map = map[map$Subject.ID == "IMP.000",], dm=bc_dm,
                            ref.sample.order = 1, xlab="Day",  ylab="Distance to Self at Day 1", outputfn="Self_to_Day1_bc.pdf",saveplot=F)
    dev.off()
    
# bacteroides to prevotella
plot.b.p.ratio(map[cs_samples,], taxa[cs_samples,])

prediction <- predict.years(map[map$Ethnicity=="Karen" & (is.na(map$Sample.Order) | map$Sample.Order==1),], taxa)

######
# plot fraction of sequences hitting GG
# fractionfn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/fraction_hitting_GG97.txt"
# plot.fraction.hits(fractionfn, map)

# plot intra-inter group variabilities
plot.group.distances(map, dm, between.groups.fn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/between.groups.txt")

plot.nutrients(map0=map[map$Years.in.US > 0 & map$Years.in.US < 15,], ethnicity=c("Hmong","Karen"), nutrient.vars=c("Added.Sugars.in.Grams",
    "Total.Calories", "Protein.in.Grams"),
    nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")
plot.nutrients(map0=map[map$Years.in.US > 0 & map$Years.in.US < 15,], ethnicity=c("Hmong","Karen"),
    nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")

plot.nutrients(map0=map[map$Years.in.US > 0,], ethnicity="Hmong", nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")
plot.nutrients(map0=map[map$Years.in.US > 0 & map$Years.in.US <=10,], ethnicity="Hmong", nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")
plot.nutrients(map0=map, ethnicity="Karen", nutrientsfn="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/Nutrients_101716.txt")

