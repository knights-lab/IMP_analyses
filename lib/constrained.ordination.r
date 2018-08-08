require(cowplot)
require(ggplot2)
require(faraway)
require(vegan)


# save.plot.rda(map0=data.frame(map[valid_cs,], food.pc[valid_cs,]), otu0=taxa_clr_L7[valid_cs,], groupvar="Sample.Group", 
#     group.cols=c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"), plot.title="Euclidean - CLR - Constrained Vars", 
#     env.vars=c("Years.in.US","BMI","Age", "Ethnicity"))    
# 
# save.plot.dbrda(map0=data.frame(map[valid_cs,], food.pc[valid_cs,]), dm0=wuf_dm[valid_cs,], groupvar="Sample.Group", 
#     group.cols=c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"), plot.title="WUF - Constrained Vars", 
#     env.vars=c("Years.in.US","BMI","Age", "Ethnicity"))    

# run
# save.plot.pca(map0=data.frame(map[valid_cs,], food.pc[valid_cs,]), otu0=taxa_clr_L7[valid_cs,], groupvar="Sample.Group", 
#     group.cols=c("#e9a3c9", "#fee08b", "#c51b7d", "#80cdc1", "#018571"), plot.title="PCA - Constrained Vars")    


save.plot.rda <- function(map0, otu0, groupvar, group.cols, method="euclidean", plot.title, env.vars)
{
    p <- ordinate.linear(map0=map0, otu0=otu0, groupvar=groupvar, group.cols=group.cols, method=method, env.vars=env.vars)
    save_plot(paste0(plot.title,".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}

save.plot.dbrda <- function(map0, dm0, groupvar, group.cols, plot.title, env.vars)
{
    p <- ordinate.dist(map0=map0, dm0=dm0, groupvar=groupvar, group.cols=group.cols, env.vars=env.vars)
    save_plot(paste0(plot.title, ".pdf"), p, useDingbats=FALSE, base_aspect_ratio = 1.3 )
}

save.plot.pca <- function(map0, otu0, groupvar, group.cols, method="euclidean", plot.title)
{
    save.plot.rda(map0=map0, otu0=otu0, groupvar=groupvar, group.cols=group.cols, method=method, plot.title=plot.title, env.vars=NULL)
}

# add.samples.dm: allows for additional samples to be included in dm that are not in mapping
prep.dm.otu <- function(map0, otu0, method="euclidean", add.samples.dm=NULL)
{
    valid_samples <- intersect(rownames(otu0), rownames(map0))
    map0 <- map0[valid_samples,]
    if(!is.null(add.samples.dm)) 
        valid_samples <- c(valid_samples, add.samples.dm)
    otu0 <- otu0[valid_samples,]        
    ddm <- vegdist(otu0, method=method)
    dm <- as.matrix(ddm) # for use in stats later
    return(list(map0=map0, ddm=ddm, dm=dm))
}

prep.dm.dist <- function(map0, dm, add.samples.dm=NULL)
{
    valid_samples <- intersect(rownames(dm), rownames(map0))
    map0 <- map0[valid_samples,]
    if(!is.null(add.samples.dm)) 
        valid_samples <- c(valid_samples, add.samples.dm)
    dm <- dm[valid_samples, valid_samples]
    ddm <- as.dist(dm)
    return(list(map0=map0, ddm=ddm, dm=dm))
}

# basic RDA with euc distance with or without vectors
ordinate.linear <- function(map0, otu0, groupvar, group.cols, method="euclidean", plot.title, env.vars=NULL)
{
    ret <- prep.dm.otu(map0=map0, otu0=otu0, method=method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0
    
    if(is.null(env.vars)) rrda <- rda(otu0)
    else rrda <- rda(as.formula(paste("otu0 ~ ", paste(env.vars, collapse=" + "))), map0)

    p <- plot.ordination(rrda=rrda, group=map0[,groupvar], group.cols=group.cols, env.vars=env.vars)

    return(p)
}

# distance based ordination
ordinate.dist <- function(map0, dm0, groupvar, group.cols, plot.title, env.vars=NULL)
{
    ret <- prep.dm.dist(map0, dm0)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0

    rrda <- dbrda(as.formula(paste("ddm ~ ", paste(env.vars, collapse=" + "))), map0)

    p <- plot.ordination(rrda=rrda, group=map0[,groupvar], group.cols=group.cols, env.vars=env.vars)

    return(p)
}

# plot constrained ordination either by passing in a pre-calculated distance based matrix OR an OTU table
plot.ordination <- function(rrda, group, group.cols, env.vars=NULL)
{   
    rrda.plot <- plot(rrda)
    
    pc <- scores(rrda)$sites
    
    d <- data.frame(x = pc[,1], y = pc[,2], group=group)
    d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))
    group.cols <- alpha(group.cols,.8)

    # calculate % variation explained by each PC
    percent_var <- signif(eigenvals(rrda)/sum(eigenvals(rrda)), 4)*100
    
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=group), size=2) +
        xlab(paste0("Axis 1 [",percent_var[1],"%]")) +
        ylab(paste0("Axis 2 [",percent_var[2],"%]")) + 
        scale_color_manual(values=group.cols) + 
        stat_ellipse(data=d, aes(colour=group), show.legend=F, type="t", level=.6)

    if(!is.null(env.vars)) 
    {
        print(anova(rrda, by="terms")) # significant terms here contribute most to constrained model
    
        # total variation explained by environment vars (sum 1st and 2nd axes only)
        label <- paste0("Variation explained: ", signif(rrda$CCA$tot.chi/rrda$tot.chi, 4)*100, "%")
    
        vector.multiplier <- attributes(rrda.plot$biplot)$arrow.mul # for redrawing length of arrow

        fit_env_df <- as.data.frame(rrda.plot$biplot*vector.multiplier)
        fit_env_df <- cbind(fit_env_df, Variable = rownames(rrda.plot$biplot))
        colnames(fit_env_df)[1:2] <- c("Dim1","Dim2")

        text_fit_env_df <- fit_env_df
        text_fit_env_df[,2] <- text_fit_env_df[,2]*1.2 # let's shift the text out a little
    
        p <- p + coord_fixed() + ## need aspect ratio of 1 in order to redraw vectors
        geom_segment(data = fit_env_df,
                   aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
                   arrow = arrow(length = unit(0.25, "cm")), colour = alpha("black",.5)) +
        geom_text(data = text_fit_env_df, aes(x = Dim1, y = Dim2, label = Variable),
                size = 3)
        
        p <- ggdraw(p) + draw_figure_label(label, size=8, position="bottom.right")
    }
    return(p)
}