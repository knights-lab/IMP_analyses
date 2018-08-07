require(hexbin)
# procrustes, for associating food and microbiome distances

# do NOT name this plot.procrustes (overrides base)
# rotate.ix = which dm to rotate
# dm.list = list of distance matrices, named
draw.procrustes <- function(map0, plot.title, max.axis = 2, dm1, dm2, dm1.name, dm2.name)
{   
    ret <- prep.dm(map0, dm=dm1)
    ddm1 <- ret$ddm # vegdist
    dm1 <- ret$dm
    map1 <- ret$map0
    ppc1 <- cmdscale(ddm1,max.axis, eig=TRUE) # return eigenvals for calculating % var explained
    pc1 <- ppc1$points

    ret <- prep.dm(map0, dm=dm2)
    ddm2 <- ret$ddm # vegdist
    dm2 <- ret$dm
    map2 <- ret$map0
    ppc2 <- cmdscale(ddm2,max.axis, eig=TRUE) # return eigenvals for calculating % var explained
    pc2 <- ppc2$points

    if(nrow(map1) != length(intersect(rownames(map1), rownames(map2))) | nrow(map2) != length(intersect(rownames(map1), rownames(map2))))
        print("Error: PCs contain different samples")

    # reorder samples in both PCs to be the same
    samples.ord <- rownames(pc1)
    pc2 <- pc2[samples.ord,]

    proc <- procrustes(X=pc1, Y=pc2, symmetric=TRUE) # Y is rotated (in this plot, its on the inside)
    # protest is always symmetric
########
    pro.t <- protest(pc1, pc2, permutations = how(nperm = 10))
    
    # find correlation between Food and MB distances
    valid <- intersect(rownames(dm1),rownames(dm2))
    man <- mantel(as.dist(dm1[valid,valid]), as.dist(dm2[valid,valid]))
    
    label <- paste0("Procrustes m2=", signif(proc$ss,2), " P=", signif(pro.t$signif,2), "\nMantel r = ", signif(man$statistic,2), " P=", man$signif)

    proc1 <- data.frame(PC1=proc$X[,1],PC2=proc$X[,2], type=dm1.name, Sample.Group=map1[samples.ord,"Sample.Group"], Sample.ID=samples.ord, stringsAsFactors=F)
    proc2 <- data.frame(PC1=proc$Yrot[,1],PC2=proc$Yrot[,2], type=dm2.name, Sample.Group=map1[samples.ord,"Sample.Group"], Sample.ID=samples.ord, stringsAsFactors=F)
    proc.df <- rbind(proc1, proc2)

    cols <- get.group.colors.distinct()
    sizes <- c(.5, .75); names(sizes) <- c(dm1.name, dm2.name)
    
    p <- ggplot(proc.df) +
      geom_point(alpha=.7, aes(x = PC1, y = PC2, color = Sample.Group, shape=type, size=type)) +
      geom_line(aes(x = PC1, y = PC2, group=Sample.ID, color=Sample.Group), alpha=.2) +
      scale_color_manual(name="Groups", values=cols) + 
      scale_size_manual(name="Groups", values=sizes,guide="none") + 
      xlab("PC1") + ylab("PC2") + theme(legend.title = element_blank())
      # add lines below for horizontal legend on bottom
      #, legend.position="bottom", legend.key = element_rect(size = 5), legend.key.size = unit(1.5, 'lines')) + 
      # guides(colour = guide_legend(nrow = 1)) 
     
     leg <- get_legend(p)
     
     p <- ggdraw(p) + draw_figure_label(label, size=8, position="bottom.right") # + theme(legend.position='none')
     
    # save_plot(paste0(plot.title, ".legend.pdf"), leg, base_aspect_ratio=5)
    save_plot(paste0(plot.title,".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 ) 
        
}

## ---> Try using this for the above to show the scatter better as a density!

# Default plot 
# c + geom_hex()
# Change the number of bins
# c + geom_hex(bins = 10)
# see other examples here: https://goo.gl/mLTNvX


# same as above, but instead of plotting procrustes plot, do a direct scatter of the two distance matrices
plot.two.dms <- function(map0, otu1=NULL, otu2=NULL, method1="euclidean", method2="euclidean", plot.title, axis1 = 1, axis2 = 2, dm1, dm2, flip.axis=NULL, xlab, ylab)
{   
    ret <- prep.dm(map0, otu1, dm=dm1, method=method1)
    ddm1 <- ret$ddm # vegdist
    dm1 <- ret$dm
    map1 <- ret$map0
    ppc1 <- cmdscale(ddm1,max(axis1,axis2), eig=TRUE) # return eigenvals for calculating % var explained
    pc1 <- ppc1$points

    ret <- prep.dm(map0, otu2, dm=dm2, method=method2)
    ddm2 <- ret$ddm # vegdist
    dm2 <- ret$dm
    map2 <- ret$map0
    ppc2 <- cmdscale(ddm2,max(axis1,axis2), eig=TRUE) # return eigenvals for calculating % var explained
    pc2 <- ppc2$points

    if(nrow(map1) != length(intersect(rownames(map1), rownames(map2))) | nrow(map1) != length(intersect(rownames(map1), rownames(map2))))
        print("Error: PCs contain different samples")
    
    # find correlation between Food and MB distances
    valid <- intersect(rownames(dm1),rownames(dm2))
    man <- mantel(as.dist(dm1[valid,valid]), as.dist(dm2[valid,valid]), permutations=999)

    label <- paste0("Mantel r = ", signif(man$statistic,2), " P=", man$signif)

 #   random.samples <- sample(valid)
 #    dm1.m <- as.matrix(dm1[random.samples,random.samples])
#     random.samples2 <- sample(valid)
#     dm2.m <- as.matrix(dm2[random.samples2,random.samples2])

     dm1.m <- as.matrix(dm1[valid,valid])
     dm2.m <- as.matrix(dm2[valid,valid])

    ggdata <- data.frame(x=dm1.m[upper.tri(dm1.m)], y=dm2.m[upper.tri(dm2.m)])

    p <- ggplot(ggdata, aes(x=x,y=y)) +
      geom_hex() + 
      xlab(xlab) + ylab(ylab) + theme(legend.title = element_blank())
     
     p <- ggdraw(p) + draw_figure_label(label, size=8, position="bottom.right")
     
    save_plot(paste0(plot.title,".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 ) 
        
}