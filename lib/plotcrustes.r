# https://github.com/RRShieldsCutler/plotCrustes/edit/master/plotcrustes_Food_v_MB_current.R
# 5/15/18
# Generate a more logical plot for procrustes comparisons
### By Robin Shields-Cutler
### Feburary 2018

library(ape)
library(vegan)
library(ggplot2)

setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")

# Read in beta diversity tables
unw_B <- read.delim('data/denovo/unweighted_unifrac_dm.txt',
                    header=1, row.names = 1, check.names = F)
unw_A <- read.delim('data/food/uwuf_food_dm.txt',
                    header=1, row.names = 1, check.names = F)
# Metadata needs three columns - 
# 1. sample IDs (e.g. "P21_pre_stool")
# 2. the unifying participant/unit ID (e.g."P21") 
# 3. the binary metadata group (e.g. "pre")
meta <- read.delim('SampleIDs-CS.txt',
                   header=F, check.names = F)
colnames(meta) <- 'participant'
meta_food <- meta

meta$data_type <- 'microbiome'
meta$sampleID <- unlist(lapply(X=meta$participant, FUN = function(xx) paste0(xx, '_microbiome')))
meta_food$data_type <- 'food'
meta_food$sampleID <- unlist(lapply(X=meta_food$participant, FUN = function(xx) paste0(xx, '_food')))

meta <- rbind(meta, meta_food)

#meta <- tibble::column_to_rownames(meta, var = 'sampleID')
rownames(meta) <- meta$sampleID
meta <- meta[,-which(colnames(meta)=="sampleID")]

thing <- 'participant'
group <- 'data_type'
groups <- c('microbiome','food')  # Names of the two distance matrices/categories
metaA <- meta[meta[,group] == groups[2],]  # Split the metadata into the two groups
metaB <- meta[meta[,group] == groups[1],]

metaA <- metaA[order(metaA[,thing]),]  # Sort the metadata by unifying ID
metaB <- metaB[order(metaB[,thing]),]

# CRITICAL:
# Ensure that the original distance matrices are in the same order by participant
# Uses the order generated from the metadata dataframe
rownames(unw_B) <- unlist(lapply(X=rownames(unw_B), FUN = function(xx) paste0(xx, '_microbiome'))) 
colnames(unw_B) <- unlist(lapply(X=colnames(unw_B), FUN = function(xx) paste0(xx, '_microbiome')))
rownames(unw_A) <- unlist(lapply(X=rownames(unw_A), FUN = function(xx) paste0(xx, '_food')))
colnames(unw_A) <- unlist(lapply(X=colnames(unw_A), FUN = function(xx) paste0(xx, '_food')))

unw_A <- unw_A[rownames(metaA),rownames(metaA)]
unw_B <- unw_B[rownames(metaB),rownames(metaB)]

# Get the principal coordinates
pcoa_A <- pcoa(unw_A)$vectors
for(c in 1:ncol(pcoa_A)){
  colnames(pcoa_A)[c] <- paste0("PC",c)
}
pcoa_B <- pcoa(unw_B)$vectors
for(c in 1:ncol(pcoa_B)){
  colnames(pcoa_B)[c] <- paste0("PC",c)
}

set.seed(11)
crusty <- procrustes(pcoa_A, pcoa_B, symmetric = T)  # Run Procrustes
set.seed(11)
crust_test_p <- protest(pcoa_A, pcoa_B, permutations = how(nperm = 999))$signif
A_crust <- data.frame(crusty$X)  # Recover the first group's coordinates
B_crust <- data.frame(crusty$Yrot)  # Recover the second group's coordinates
colnames(B_crust) <- colnames(A_crust)
ncoords = as.numeric(ncol(A_crust))
A_crust <- merge(A_crust, metaA, by=0)
B_crust <- merge(B_crust, metaB, by=0)
rownames(A_crust) <- A_crust[,1]; A_crust[,1] <- NULL  # Merge makes rownames column
rownames(B_crust) <- B_crust[,1]; B_crust[,1] <- NULL
sample_ids <- A_crust[,thing]  # Get all the unifying participant IDs

real_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))  # Initialize the dataframe
colnames(real_dist) <- c('sampleID_A', 'sampleID_B', 'distance')
# Loop through each participant to get the multidimensional distance between their rotated points
for (i in 1:length(sample_ids)) {
  ix <- as.character(sample_ids[i])
  A_ix <- A_crust[A_crust[,thing] == ix, 1:ncoords]  # Keep all the PC axes
  B_ix <- B_crust[B_crust[,thing] == ix, 1:ncoords]
  AB_mat <- rbind(A_ix, B_ix)
  AB_dist <- matrix(dist(AB_mat, method = 'euclidean'))  # Calculates the distance
  real_dist[i,1] <- as.character(rownames(A_ix))  # Fill in the dataframe
  real_dist[i,2] <- as.character(rownames(B_ix))
  real_dist[i,3] <- as.numeric(AB_dist[1])
}

A_crust$data_type <- 'food'; B_crust$data_type <- 'microbiome'
A_crust$realperm <- 'real'; B_crust$realperm <- 'real'

# The "Procrustes Distance"
reaL_pro_dist <- sqrt(sum(real_dist$distance^2))
cat(reaL_pro_dist)


### Multiple permutations ####

# Repeat the permutation "j" times
real_dist$real_perm <- 'true_distance'
dist_plot <- real_dist
PCOA_plot_many <- rbind(A_crust, B_crust)  # Start with the real data, then add the permutations
procrust_pvals <- vector(mode='numeric', length=9)
procrust_dists <- vector(mode='numeric', length=9)
for (j in 1:9) { 
  metaBp <- metaB
  set.seed(j)
  metaBp[,thing] <- sample(x = metaBp[,thing], size = length(metaBp[,thing]), replace = F)
  metaBp <- metaBp[order(metaBp[,thing]),]
  unw_Bp <- unw_B[rownames(metaBp),rownames(metaBp)]
  pcoa_Bp <- pcoa(unw_Bp)$vectors
  for(c in 1:ncol(pcoa_Bp)){
    colnames(pcoa_Bp)[c] <- paste0("PC",c)
  }
  set.seed(j)
  crusty_p <- procrustes(pcoa_A, pcoa_Bp, symmetric = T)
  set.seed(j)
  crust_test_perm_p <- protest(pcoa_A, pcoa_Bp, permutations = how(nperm = 999))$signif
  cat(paste0(crust_test_perm_p, '\n'))
  procrust_pvals[j] <- crust_test_perm_p
  Ap_crust <- data.frame(crusty_p$X)
  Bp_crust <- data.frame(crusty_p$Yrot)
  colnames(Bp_crust) <- colnames(Ap_crust)
  Ap_crust <- merge(Ap_crust, metaA, by=0)
  Bp_crust <- merge(Bp_crust, metaBp, by=0)
  rownames(Ap_crust) <- Ap_crust[,1]; Ap_crust[,1] <- NULL
  rownames(Bp_crust) <- Bp_crust[,1]; Bp_crust[,1] <- NULL
  sample_ids <- Ap_crust[,thing]
  
  perm_dist <- data.frame(matrix(nrow = length(sample_ids), ncol = 3))
  colnames(perm_dist) <- c('sampleID_A', 'sampleID_B', 'distance')
  for (i in 1:length(sample_ids)) {
    ix <- as.character(sample_ids[i])
    # ixp <- as.character(sample_id_perm[i])
    Ap_ix <- Ap_crust[Ap_crust[,thing] == ix, 1:ncoords]
    Bp_ix <- Bp_crust[Bp_crust[,thing] == ix, 1:ncoords]
    ABp_mat <- rbind(Ap_ix, Bp_ix)
    ABp_dist <- matrix(dist(ABp_mat, method = 'euclidean'))
    perm_dist[i,1] <- as.character(rownames(Ap_ix))
    perm_dist[i,2] <- as.character(rownames(Bp_ix))
    perm_dist[i,3] <- as.numeric(ABp_dist[1])
  }
  perm_dist_plot_2 <- perm_dist
  perm_dist_plot_2$real_perm <- paste0('permuted',j)  # Keep track of the permuted data
  dist_plot <- rbind(dist_plot, perm_dist_plot_2)
  
  # Add the permuted data to the existing runs
  Ap_crust$data_type <- 'food'; Bp_crust$data_type <- 'microbiome'
  Ap_crust$realperm <- paste0('permuted',j); Bp_crust$realperm <- paste0('permuted',j)
  PCOA_plot_many <- rbind(PCOA_plot_many, Ap_crust, Bp_crust)
  
  # Save the procrustes distance of the overall dataset
  procrust_dist_perm <- sqrt(sum(perm_dist$distance^2))
  procrust_dists[j] <- procrust_dist_perm
}

# What's the median permutation?
cat(procrust_pvals)
cat(procrust_dists)
median(procrust_dists)

### Plot the results ####

cols <- c("#A300FF", "#FF0000", "#00A696",  "#FE42AD", "#FBB400", "#2E1915")
names(cols) <- c("KarenThai", "HmongThai", "Karen1st",  "Hmong1st",  "Hmong2nd",  "Control" )

# Get more metadata for coloring
meta.big <- read.delim('data/mapping.txt', header=1, sep = '\t', check.names = F)
meta.big <- meta.big[meta.big$'#SampleID' %in% as.character(sample_ids), ]

dist_plot$sampleID <- unlist(lapply(X=dist_plot$sampleID_B,
                                    FUN=function(xx) gsub(pattern = '_microbiome', replacement = '', x = xx, fixed = T)))
dist_plot_groups <- merge(dist_plot, meta.big[,c('#SampleID','Sample.Group')], by.x = 'sampleID', by.y='#SampleID')
### PJ
temp <- meta.big[,c('#SampleID','Sample.Group')]
temp[,'#SampleID'] <- as.character(temp[,'#SampleID'])
dist_plot_groups2 <- merge(dist_plot, temp, by.x = 'sampleID', by.y='#SampleID')
###
dist_plot_groups$Sample.Group <- factor(dist_plot_groups$Sample.Group,
                                        levels = c("KarenThai", "HmongThai", "Karen1st",  "Hmong1st",  "Hmong2nd",  "Control"),
                                        ordered = T)

# Boxplot with scatter ALL PERMUTATIONS
pbox <- ggplot(dist_plot, aes(x=real_perm, y=distance, group=real_perm)) +
  geom_boxplot(outlier.colour = 'white') + geom_jitter(width = 0.2, alpha=0.20) + theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        axis.text = element_text(color='black'))
# To save
ggsave(pbox, filename = 'distance_boxplots_10x.png', width = 10, height = 6, dpi = 300)

# Set of typical procrustes plots ALL PERMUTATIONS
pscat <- ggplot(PCOA_plot_many) + geom_point(aes(x=PC1, y=PC2, color=data_type), size=0.7, alpha=0.6) +
  geom_line(aes(x=PC1, y=PC2, group=participant), alpha=0.1) +
  theme_classic() + facet_grid(. ~ realperm) +
  theme(axis.text = element_text(color='black', size = 6)) + coord_fixed(ratio = 1)
# To save
ggsave(pscat, filename = 'permute10x_groups.png', width = 15, height = 3, dpi = 300)


#### Set of two procrustes plots WITH GROUPS ####
PCOA_plot_many.groups <- merge(PCOA_plot_many, meta.big[,c('#SampleID','Sample.Group')], by.x = 'participant', by.y='#SampleID')
PCOA_plot_many.groups$Sample.Group <- factor(PCOA_plot_many.groups$Sample.Group,
                                        levels = c("KarenThai", "HmongThai", "Karen1st",  "Hmong1st",  "Hmong2nd",  "Control"),
                                        ordered = T)

# Enter the permuted set that features the median procrustes distance
PCOA_plot_many.groups.sub <- PCOA_plot_many.groups[PCOA_plot_many.groups$realperm == 'permuted9' | PCOA_plot_many.groups$realperm == 'real', ]
PCOA_plot_many.groups.sub$realperm[PCOA_plot_many.groups.sub$realperm == 'permuted9'] <- 'Permuted'
PCOA_plot_many.groups.sub$realperm[PCOA_plot_many.groups.sub$realperm == 'real'] <- 'Original'
PCOA_plot_many.groups.sub$realperm <- factor(x = PCOA_plot_many.groups.sub$realperm, levels = c('Original','Permuted'), ordered = T)

pscat.groups <- ggplot(PCOA_plot_many.groups.sub) +
  geom_line(aes(x=PC1, y=PC2, group=participant, color=Sample.Group), alpha=0.3) +
  geom_point(aes(x=PC1, y=PC2, color=Sample.Group, shape=data_type), size=1.2, alpha=0.7) +
  theme_classic() + facet_grid(. ~ realperm) +
  theme(axis.text = element_text(color='black', size = 6),
        axis.title = element_text(size=6, color='black')) +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = cols) + theme(legend.title = element_blank())
# To save
ggsave(pscat.groups, filename = '../IMP/results/procrustes_vs_med_permuted_group_colored_v3.png', width = 8.5, height = 4, dpi = 300)


dist_plot_groups.sub <- dist_plot_groups[dist_plot_groups$real_perm == 'permuted9' | dist_plot_groups$real_perm == 'true_distance', ]
dist_plot_groups.sub$real_perm[dist_plot_groups.sub$real_perm == 'permuted9'] <- 'permuted'
dist_plot_groups.sub$real_perm[dist_plot_groups.sub$real_perm == 'true_distance'] <- 'original'
dist_plot_groups.sub$real_perm <- factor(x = dist_plot_groups.sub$real_perm, levels = c('original','permuted'), ordered = T)


pbox.group <- ggplot(dist_plot_groups.sub, aes(x=real_perm, y=distance, group=real_perm)) +
  geom_boxplot(outlier.colour = 'white') +
  geom_jitter(width = 0.28, alpha=0.45) +
  theme_classic() +
  scale_color_manual(values=cols) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color='black', size=10),
        axis.title.y = element_text(color='black',size=12))
# To save
ggsave(pbox.group, filename = '../IMP/results/distance_boxplots_vs_med_permutation_group_v2.png', width = 3.5, height = 5, dpi = 300)

shapiro.test(dist_plot_groups.sub[dist_plot_groups.sub$real_perm=='original', 'distance'])
wilcox.test(dist_plot_groups.sub$distance ~ dist_plot_groups.sub$real_perm, paired=F)$p.value
