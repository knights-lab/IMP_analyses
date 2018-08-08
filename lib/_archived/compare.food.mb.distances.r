# get correlation by transforming matrix to a vector
 food.sample.ids <- intersect(rownames(food_wuf_dm[cs,]), rownames(wuf_dm[cs,]))
    food_wuf_dm0 <- food_wuf_dm[food.sample.ids,food.sample.ids]
    wuf_dm0 <- wuf_dm[food.sample.ids,food.sample.ids]
    
    cor.test(unlist(wuf_dm0), unlist(food_wuf_dm0))
    
# Pearson's product-moment correlation
# 
# data:  unlist(wuf_dm0) and unlist(food_wuf_dm0)
# t = 46.994, df = 91807, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1469419 0.1595751
# sample estimates:
#       cor 
# 0.1532647 

# compare distance matrices directly
mantel(as.dist(food_wuf_dm0), as.dist(wuf_dm0))

# Mantel statistic based on Pearson's product-moment correlation 
# 
# Call:
# mantel(xdis = as.dist(food_wuf_dm0), ydis = as.dist(wuf_dm0)) 
# 
# Mantel statistic r: 0.1154 
#       Significance: 0.003 
# 
# Upper quantiles of permutations (null model):
#    90%    95%  97.5%    99% 
# 0.0454 0.0627 0.0731 0.0916 
# Permutation: free
# Number of permutations: 999

# generate a heatmap of mb distances vs food distances

plot.nutrient.mb.heatmap <- function(nutrients, taxa, map0, fn)
{
    # samples <- intersect(rownames(nutrients), rownames(map0))
#     map0 <- map0[samples,]
#     nutrients0 <- nutrients[samples,]
#     taxa0 <- taxa[samples,]
#     
    #nutrients0 <- nutrients0[,nutrient.vars]
   
    m <- cbind(unlist(wuf_dm0), unlist(food_wuf_dm0))
    colnames(m) <- c("MB", "FU")
#     rcorr.m <- rcorr(m, type="spearman")
#     corr.m <- rcorr.m$r
#     # let's filter and reorder our samples here first before melting 
#     rcorr.t1 <- rcorr(as.matrix(taxa0), type="s")
#     rcorr.t <- rcorr.t1$r
#     t.ordered.labels <- reorder_cormat(rcorr.t)
#      # now filter the full correlation matrix before melting. don't bother with nutrients, just taxa.
#     cormat <- corr.m[colnames(nutrients0), t.ordered.labels]
     
    cordata <- melt(m, na.rm=T)    
    
# Create a ggheatmap
ggheatmap <- ggplot(cordata, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white") + xlab("") + ylab("") + 
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
    name="Spearman's\nCorrelation") +
  theme_minimal()+ # minimal theme
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 5, hjust = 1), axis.text.y = element_text(size=5)) +
 coord_fixed()

    save_plot(fn, ggheatmap, useDingbats=FALSE, base_aspect_ratio = 1.3 )

}