#############
#Procrustes of transcripts and beta div uni
beta <-  as.dist(uni[samples2, samples2])
samples3 <- mapping[samples2,"transcriptomic"]
trans1 <- ileum_transcript[,samples3]
colnames(trans1) <- colnames(beta)
trans <- vegdist(t(trans1))

PCOA1 <- pcoa(beta)$vectors
for(c in 1:ncol(PCOA1)){
  colnames(PCOA1)[c] <- paste("PC",c, sep="")
}

PCOA2 <- pcoa(trans)$vectors
for(c in 1:ncol(PCOA2)){
  colnames(PCOA2)[c] <- paste("PC",c, sep="")
}

pro <- procrustes(PCOA1, PCOA2)
pro2 <- protest(PCOA1, PCOA2, permutations = how(nperm = 999))
beta_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)
rownames(trans_pro) <- samples3
colnames(trans_pro) <- colnames(beta_pro)
trans_pro$treatment <- mapping[samples2, "treatment4"]
trans_pro$type <- "transcript"
beta_pro$treatment <- mapping[samples2, "treatment4"]
beta_pro$type <- "otu"
beta_pro$sample <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)
trans_pro$sample <- c(1,2,3,4,5,6,7,8,9,10,11,12,13)

PCOA <- rbind(beta_pro, trans_pro)
PCOA$SampleID <- rownames(PCOA)

ggplot(PCOA) +
  geom_point(size = 4, alpha=0.75, aes_string(x = "PC1", y = "PC2", color = "treatment", shape="type")) +
  theme_bw() +
  geom_line(aes(x= PC1, y=PC2, group=sample, color=treatment)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        legend.title = element_blank(),
        legend.key.size = unit(0.2, "in"),
        legend.text = element_text(size=5),
        legend.position = 'bottom',
        axis.text = element_text(size=5),
        axis.title = element_text(size=8))