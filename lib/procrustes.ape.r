# test out procrustes with ape

    pcoa1 <- pcoa(as.matrix(food_uwuf_dm[cs,cs]))
    pcoa2 <- pcoa(as.matrix(uwuf_dm[cs,cs]))
    
    pc1 <- pcoa1$vectors
    pc2 <- pcoa2$vectors
    pc2 <- pc2[rownames(pc1),]
    
    proc <- procrustes(pc1,pc2, symmetric=TRUE)

    proc1 <- data.frame(PC1=proc$X[,1],PC2=proc$X[,2], type="Food", Sample.Group=map[rownames(pc1),"Sample.Group"], Sample.ID=rownames(pc1), stringsAsFactors=F)
    proc2 <- data.frame(PC1=proc$Yrot[,1],PC2=proc$Yrot[,2], type="MB", Sample.Group=map[rownames(pc1),"Sample.Group"], Sample.ID=rownames(pc1), stringsAsFactors=F)
    proc.df <- rbind(proc1, proc2)

    cols <- get.group.colors.distinct()
    sizes <- c(.5, .75); names(sizes) <- c("MB","Food")
    
    p <- ggplot(proc.df) +
      geom_point(alpha=.7, aes(x = PC1, y = PC2, color = Sample.Group, shape=type, size=type)) +
      geom_line(aes(x = PC1, y = PC2, group=Sample.ID, color=Sample.Group), alpha=.2) +
      scale_color_manual(name="Groups", values=cols) + 
      scale_size_manual(name="Groups", values=sizes,guide="none") + 
      xlab("PC1") + ylab("PC2") + theme(legend.title = element_blank())
          
    save_plot("procrustes-ape1.pdf", p, useDingbats=FALSE, base_aspect_ratio = 1.3 )