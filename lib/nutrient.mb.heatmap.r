library(ggplot2)
library(reshape2)

distfun <- function(x) { as.dist((1-cor(t(x), method="spearman"))/2)}

plot.nutrient.mb.heatmap <- function(nutrients, taxa, map0, fn)
{
    samples <- intersect(rownames(nutrients), rownames(map0))
    map0 <- map0[samples,]
    nutrients0 <- nutrients[samples,]
    taxa0 <- taxa[samples,]
    
    nutrient.vars <- c(
   "Total Calories",
   "% of Calories from Total Sugars",
    "% of Calories from Added Sugars",
  "Fruits in Cups",

"Vegetables in Cups",
"Dietary Fiber in Grams",
 
 "Total Fat (% Calories Eaten )",
 "% of Calories from Saturated Fat",
 "% of Calories from Protein")
   
    nutrients0 <- nutrients0[,nutrient.vars]
   
    m <- as.matrix(cbind(nutrients0,taxa0))
    rcorr.m <- rcorr(m, type="spearman")
    corr.m <- rcorr.m$r
    # let's filter and reorder our samples here first before melting 
    rcorr.t1 <- rcorr(as.matrix(taxa0), type="s")
    rcorr.t <- rcorr.t1$r
    t.ordered.labels <- reorder_cormat(rcorr.t)
     # now filter the full correlation matrix before melting. don't bother with nutrients, just taxa.
    cormat <- corr.m[colnames(nutrients0), t.ordered.labels]
     
    cordata <- melt(cormat, na.rm=T)    
    
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


reorder_cormat <- function(cormat){
# Use correlation between variables as distance
dd <- as.dist((1-cormat)/2)
hc <- hclust(dd)
rownames(cormat[hc$order, hc$order])
}