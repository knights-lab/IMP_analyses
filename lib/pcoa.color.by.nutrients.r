# plots a PCOA and then generates multiple plots and colors by nutrient consumption to see if there are any differences in groups

plot.pcoa.nutrient <- function(map0, pc0, nutrients0, ethnicity, fn)
{
    nutrient.vars <- c(
        "Total Calories",
        "% of Calories from Total Sugars",
        "% of Calories from Added Sugars",
        "Fruits in Cups",
        "Vegetables in Cups",
        "g Fiber per 1000 Calories",
        "% of Calories from Total Fat",
        "% of Calories from Saturated Fat",
        "% of Calories from Protein",
        "% of Calories from Carbohydrate",
        "Grains in Ounces",
        "Dairy in Cups")

    valid_samples <- intersect(rownames(nutrients0), rownames(map0))
    map0 <- map0[valid_samples,]
    map0 <- map0[map0$Ethnicity %in% ethnicity,]
    
    # now grab only samples that are present in the dm
    valid_samples <- intersect(rownames(map0), rownames(pc0))
    
    map0 <- map0[valid_samples,]
    nutrients0 <- nutrients0[valid_samples,nutrient.vars]
    pc0 <- pc0[valid_samples,]
    
    
    # 17:triangle = Thai, 19:Circle = 1st, 15:square = 2nd
    lookup <- c(17,17,19,19,15) # point type
    names(lookup) <- c("HmongThai", "KarenThai", "Hmong1st", "Karen1st", "Hmong2nd")
    pch <- lookup[as.character(map0$Sample.Group)] 

    # color by nutrient
    nutrients0 <- log10(nutrients0 + 1)

    p <- NULL
    for(i in 1:length(nutrient.vars))
    {
        d <- data.frame(x = pc0[,1], y = pc0[,2], group=map0$Sample.Group, nutrient=nutrients0[,nutrient.vars[i]])

        p[[i]] <- ggplot(data=d, aes(x, y)) + geom_point(aes(color=nutrient), shape=pch) +
            scale_color_gradient(low=alpha("yellow",.5), high=alpha("red",.5)) +
            xlab("PC1") + ylab("PC2") + ggtitle(nutrient.vars[i]) 
#            stat_ellipse(data=d, aes(colour=group), show.legend=F, type="t", level=.6)

        ix <- map0$Sample.Group == "Hmong2nd"
        p[[i]] <- p[[i]] + annotate("text", x = pc0[ix,1], y = pc0[ix,2], label = rownames(map0)[ix], size=1.3)


        # let's test all nutrients for separation and print it to a file
        #adonis.years <- adonis(dm[firstgen,firstgen] ~ map0[firstgen,"Years.in.US"])
        

    }
    
    multiplot <- plot_grid(plotlist=p, ncol=3, nrow=4)
    save_plot(fn, multiplot,
        ncol = 3, # we're saving a grid plot of 2 columns
        nrow = 4, # and 2 rows
        # each individual subplot should have an aspect ratio of 1.3
        base_aspect_ratio = 1.3
    )

}