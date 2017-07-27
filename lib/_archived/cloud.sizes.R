library('vegan')
library('beeswarm')
ALPHA <- .05

dir.create("dist_compare/")

# subset distance to keep only pups and dams
bd_pup <- as.matrix(bd_pup)[rownames(m_pup),rownames(m_pup)]
bd_dam <- as.matrix(bd_dam)[rownames(m_dam),rownames(m_dam)]

m_pup$MouseID <- as.character(m_pup$MouseID)
m_dam$MouseID <- as.character(m_dam$MouseID)

# get all within-mouse distances
mouse_number <- as.character(unique(m_pup$MouseID))
wds.pup <- numeric(length(mouse_number))
names(wds.pup) <- sprintf(mouse_number)
for(i in 1:length(mouse_number)){
    mouse <- mouse_number[i]
    # within-patient distances
    mouse.ix <- rownames(m_pup)[which(m_pup$MouseID == mouse)]
    if(length(mouse.ix) == 1){
      wds.pup[i] <- NA
    } else {
        # calculates all within-cloud distances
        tmp <- bd_pup[mouse.ix,mouse.ix]
        tmp <- tmp[upper.tri(tmp)]
        wds.pup[i] <- mean(tmp)
        
        # calculates only distance to previous timepoint
        # tmp <- sapply(1:(sum(patient.ix)-1), function(ixx) dgn[which(patient.ix)[ixx],which(patient.ix)[ixx+1]])
        # wds.pup[i] <- mean(tmp)
    }
}

# get all within-mouse distances (dam)
mouse_number <- as.character(unique(m_dam$MouseID))
wds.dam <- numeric(length(mouse_number))
names(wds.dam) <- sprintf(mouse_number)
for(i in 1:length(mouse_number)){
  mouse <- mouse_number[i]
  # within-patient distances
  mouse.ix <- rownames(m_dam)[which(m_dam$MouseID == mouse)]
  if(length(mouse.ix) == 1){
    wds.dam[i] <- NA
  } else {
    # calculates all within-cloud distances
    tmp <- bd_dam[mouse.ix,mouse.ix]
    tmp <- tmp[upper.tri(tmp)]
    wds.dam[i] <- mean(tmp)
    
    # calculates only distance to previous timepoint
    # tmp <- sapply(1:(sum(patient.ix)-1), function(ixx) dgn[which(patient.ix)[ixx],which(patient.ix)[ixx+1]])
    # wds.pup[i] <- mean(tmp)
  }
}

# test distances
test.xs <- list(pup=wds.pup, dam=wds.dam)
maps <- list(m_pup, m_dam)
# different group comparisons
test.ixs <- list('c57-con v. c57-stat'=c(c57.con, c57.stat),
                 'c57-con v. il10-con'=c(c57.con, il10.con),
                 'il10-con v. il10-stat'=c(il10.con, il10.stat),
                 'il10-stat v. c57-stat'=c(il10.stat, c57.stat)
                )
# compare.to <- list('HC v. IBS'='Healthy',
#                    'HC v. IBSD'='Healthy',
#                    'HC v. IBSC'='Healthy',
#                    'IBSC v. IBSD'='IBS-C'
# )

# run all combinations of tests

for(i in 1:length(test.xs)){
    x.name <- names(test.xs)[i]
    test.x <- test.xs[[i]]
    p_vals <- c()
    for(j in 1:length(test.ixs)){
        test.name <- names(test.ixs)[j]
        #compare.to.name <- compare.to[[test.name]]
        test.ix <- test.ixs[[j]]
        tmp_map <- data.frame(maps[i])
        tmp_map <- tmp_map[!duplicated(tmp_map$MouseID),]
        to_order <- names(test.x[which(names(test.x) %in% test.ix)])
        tmp_map <- tmp_map[match(to_order,tmp_map$MouseID),]
        tt <- t.test(test.x[which(names(test.x) %in% test.ix)] ~ tmp_map$TG_Mouse) 
        p_vals <- c(p_vals,tt$p.value)
        names(p_vals)[j] <- test.name
        
        if(tt$p.value < ALPHA) {
            # sink("dist_compare/sig_stats_file.txt", append=TRUE)
            # cat("\n\n")
            # cat(sprintf('t-test variability, %s, %s: ', x.name, test.name))
            # cat('p=',round(tt$p.value,4), ', t statistic=',round(tt$statistic,4),'\n',sep='')
            # sink()
            pdf(sprintf('dist_compare/variability_%s_%s.pdf',x.name, 
                        gsub(' ','_',test.name)),
                        width=4,
                        height=4, 
                        useDingbats=FALSE)
            beeswarm(test.x[which(names(test.x) %in% test.ix)] ~ droplevels(tmp_map$TG_Mouse),
                     xlab=paste(test.name,sep=''), cex.axis = 0.75,
                     ylab='mean within-mouse variability')
            bxplot(test.x[which(names(test.x) %in% test.ix)] ~ droplevels(tmp_map$TG_Mouse),add=TRUE)
            dev.off()
        } else {
          # sink("dist_compare/insig_stats_file.txt", append=TRUE)
          # cat("\n\n")
          # cat(sprintf('t-test variability, %s, %s: ', x.name, test.name))
          # cat('p=',round(tt$p.value,4), ', t statistic=',round(tt$statistic,4),'\n',sep='')
          # sink()
          pdf(sprintf('dist_compare/insig_variability_%s_%s.pdf',x.name, 
                      gsub(' ','_',test.name)),
              width=4,
              height=4, 
              useDingbats=FALSE)
          beeswarm(test.x[which(names(test.x) %in% test.ix)] ~ droplevels(tmp_map$TG_Mouse),
                   xlab=paste(test.name,sep=''), cex.axis = 0.75,
                   ylab='mean within-mouse variability')
          bxplot(test.x[which(names(test.x) %in% test.ix)] ~ droplevels(tmp_map$TG_Mouse),add=TRUE)
          dev.off()
        }
    }
    p_vals <- p.adjust(p_vals,method="fdr")
    sink("dist_compare/adjusted_pvals.txt", append=TRUE)
    cat(paste("\n", x.name, "\n"))
    print(p_vals)
    sink()
}
