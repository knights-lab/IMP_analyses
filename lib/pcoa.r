require(ggpubr) 
require(vegan)
require(faraway)
require(RColorBrewer)
require(cowplot)

# label.samples = vector of samples to label by name
# used for plotting samples for selecting Hmong Mouse Donors
plot.pcoa.by <- function(map0, otu0, dm=NULL, fill.by="Years.in.US", shape.by=NULL, label.samples=NULL, flip.axis=NULL, fill.color="#40004B", low.fill.color="#FDD023")
{
    ret <- prep.dm(map0, otu0, dm, method="euclidean")
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0
    ppc <- cmdscale(ddm,2, eig=TRUE) # return eigenvals for calculating % var explained
    pc <- ppc$points
    
    if(!is.null(flip.axis))  pc[,flip.axis] <- -1*pc[,flip.axis] # sometimes we'll want to flip certain axes so orientations are the same
    
    d <- data.frame(x = pc[,1], y = pc[,2], 
                    fill.by=map0[,fill.by], shape.by=map0[,shape.by], group=map0$Sample.Group, label=rownames(map0))

    p <- ggplot(data=d, aes(x, y)) + xlab("PC1") + ylab("PC2")

    if(!is.null(shape.by))
        p <- p + geom_point(aes(fill=fill.by, shape=shape.by), color=alpha("black",.2), size=2) + scale_shape_manual(name=shape.by, values=c(21,24)) + scale_fill_gradient(name=fill.by, low="white", high=fill.color)
    else
        p <- p + geom_point(aes(color=fill.by), size=1.5) + scale_color_gradient(name="", low=low.fill.color, high=fill.color)


    if(!is.null(label.samples))
        p <- p + geom_text(aes(label=ifelse(label %in% label.samples, as.character(label), '')), hjust=0, vjust=0)

    return(p)
}

# dm is optional and only for things that have been generated outside of R (unifrac)
# env.vars = specify continuous vars to plot additional environment variables 
# this function constrains the ordination by a list of environment variables passed in
plot.constrained.ordination <- function(map0, otu0, method="euclidean", plot.title, dm0=NULL, env.vars=NULL)
{   
    ret <- prep.dm(map0, otu0, dm0, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0

    # RDA: distance is euclidean (CLR transformed) and relationship with env.vars is expected to be linear (otherwise, consider CCA)
    if(is.null(dm0) & !is.null(env.vars))
        rrda <- rda(as.formula(paste("otu0 ~ ", paste(env.vars, collapse=" + "))), map0)    # basic RDA with euc distance 
    else if(is.null(dm0) & is.null(env.vars))
        rrda <- rda(otu0)                                                                   # basic RDA with euc distance - full model
    else if(!is.null(dm0) & !is.null(env.vars))
        rrda <- dbrda(as.formula(paste("ddm ~ ", paste(env.vars, collapse=" + "))), map0)   # distance based RDA (e.g. constrained unifrac)
    else if(!is.null(dm0) & is.null(env.vars))
        rrda <- dbrda(ddm)                                                                  # distance based RDA - full model (e.g. unifrac)
    
    rrda.plot <- plot(rrda)
    
    pc <- scores(rrda)$sites
    
    d <- data.frame(x = pc[,1], y = pc[,2], Sample.Group=map0$Sample.Group)
    # set the levels of Sample.Group so that it's the same every time
    #d$group <- factor(d$group, levels=sort(as.character(unique(d$group))))

    groupnames <- as.character(unique(map0$Sample.Group))    
    cols <- get.group.colors(groups=groupnames) 
    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 


    # calculate % variation explained by each PC
    percent_var <- signif(eigenvals(rrda)/sum(eigenvals(rrda)), 4)*100
    
    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, size=Sample.Group, alpha=Sample.Group)) +
        xlab(paste0("PC1 [",percent_var[1],"%]")) +
        ylab(paste0("PC2 [",percent_var[2],"%]")) + 
        scale_color_manual(name="Groups", values=cols) + #sets the color palette of the fill
        scale_alpha_manual(name="Groups", values=alphas) +
        scale_shape_manual(name="Groups", values=shapes) +
        scale_size_manual(name="Groups", values=sizes) + 
        stat_ellipse(data=d, aes(colour=Sample.Group), show.legend=F, type="t", level=.6)
    if(!is.null(env.vars)) # if environment variables were passed in
    {
        #print(anova(rrda, by="terms")) # significant terms here contribute most to constrained model
        # total variation explained by environment vars (sum 1st and 2nd axes only)
        label <- paste0("Variation explained: ", signif(rrda$CCA$tot.chi/rrda$tot.chi, 4)*100, "%")
        vector.multiplier <- attributes(rrda.plot$biplot)$arrow.mul # for redrawing length of arrow

        fit_env_df <- as.data.frame(rrda.plot$biplot*vector.multiplier)
        fit_env_df <- cbind(fit_env_df, Variable = rownames(rrda.plot$biplot))
        colnames(fit_env_df)[1:2] <- c("Dim1","Dim2")

        text_fit_env_df <- fit_env_df
        text_fit_env_df[,2] <- text_fit_env_df[,2]*1.2 # let's shift the text out a little
        p <- p + theme(legend.position='none') + coord_fixed() + ## need aspect ratio of 1 in order to redraw vectors!
                    ylim(c(-4,4)) + 
                    xlim(c(-4,4)) +
                    geom_segment(data = fit_env_df,
                   aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
                   arrow = arrow(length = unit(0.25, "cm")), colour = alpha("black",.5)) +
                geom_label(data = text_fit_env_df, aes(x = Dim1, y = Dim2, label = Variable),
                size = 2.5)
        
        #p <- ggdraw(p) + draw_figure_label(label, size=8, position="bottom.right")
    }
    
    return(p)

}


# dm is optional and only for things that have been generated outside of R (unifrac)
# env.vars = specify continuous vars to plot additional environment variables 
# hide.groups = plot PCOA with all samples in map0, but only show certain samples
plot.pcoa <- function(map0, otu0, method="euclidean", plot.title, axis1 = 1, axis2 = 2, dm=NULL, env.vars=NULL, flip.axis=NULL, show.stats=TRUE, save.pc=FALSE, hide.groups=NULL)
{   
    ret <- prep.dm(map0, otu0, dm, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0
    ppc <- cmdscale(ddm,max(axis1,axis2), eig=TRUE) # return eigenvals for calculating % var explained
    pc <- ppc$points
        
    if(!is.null(flip.axis))  pc[,flip.axis] <- -1*pc[,flip.axis] # sometimes we'll want to flip certain axes so orientations are the same
    
    if(show.stats)
    {
        # adonis on residence class as factors for ALL samples
        adonis.sample.group <- adonis(dm ~ map0$Sample.Group)
        adonis.sample.group.pval <- adonis.sample.group$aov.tab[1,"Pr(>F)"]     
        adonis.sample.group.R2 <- adonis.sample.group$aov.tab[1,"R2"]
        print(adonis.sample.group)

        print(adonis(dm ~ map0$Ethnicity + map0$Birth.Continent + map0$Resident.Continent))

        anosim.result <- anosim(dm, grouping=map0$Sample.Group)
        estimate <- anosim.result$statistic
        pval <- anosim.result$signif
        print(anosim.result)

        # adonis on continuous years.in.us for 1st generation ONLY
        firstgen <- rownames(map0)[map0$Sample.Group %in% c("Hmong1st","Karen1st")]    
        adonis.years <- adonis(dm[firstgen,firstgen] ~ map0[firstgen,"Years.in.US"])
        print(adonis.years)
        adonis.years.pval <- adonis.years$aov.tab[1,"Pr(>F)"]     
        adonis.years.R2 <- adonis.years$aov.tab[1,"R2"]

        cortest <- cor.test(pc[,axis1],map0[,"Years.in.US"], method="spear")
        print(cortest)
    }                     
    # let's draw 95% standard error ellipses around Thai and 2nd gen samples
#    mod <- betadisper(ddm, map0$Sample.Group, type="centroid")
#    plot(mod, ellipse = TRUE, hull = FALSE) # 1 sd data ellipse
        
    # calculate variability explained
    percent_var <- signif(eigenvals(ppc)/sum(eigenvals(ppc)), 4)*100
        
    d <- data.frame(x = pc[,axis1], y = pc[,axis2], map0)
    
    groupnames <- as.character(unique(map0$Sample.Group))
    
    cols <- get.group.colors(groups=groupnames) 
    alphas <- get.group.alphas(groups=groupnames) 
    shapes <- get.group.shapes(groups=groupnames) 
    sizes <- get.group.sizes(groups=groupnames) 

    # hide some groups from the plot
    if(!is.null(hide.groups)) { 
        alphas[hide.groups] <- 0
    }

    p <- ggplot(data=d, aes(x, y)) + geom_point(aes(colour=Sample.Group, shape=Sample.Group, size=Sample.Group, alpha=Sample.Group), stroke=.6, fill=NA) + 
        xlab(paste0("PC",axis1, " [",percent_var[axis1],"%]")) + ylab(paste0("PC",axis2, " [",percent_var[axis2],"%]")) + 
        ggtitle(plot.title) +
        scale_color_manual(name="Groups", values=cols) + #sets the color palette of the fill
        scale_alpha_manual(name="Groups", values=alphas) +
        scale_shape_manual(name="Groups", values=shapes) +
        scale_size_manual(name="Groups", values=sizes) + 
        stat_ellipse(data=d[d$Group %in% c("Pre","2nd", "Control"),], aes(colour=Sample.Group, alpha=Sample.Group), show.legend=F, type="t", level=.6)

    label<-NULL
    if(!is.null(env.vars))
    {
        # Note, this will fit Env.Var ~ PC1 and PC2 -- and report the Multiple-R-squared of this fit (permutation based)
        # This is NOT the same as correlating Env.Var directly with PC1 or PC2 only (above)
        
        # Fit environmental variables
        fit_env <- envfit(ord=pc, env=map0[,env.vars])

        # need to plot these (invisibly) first then grab scaling and vectors out
        # try plotting PCOA with env variables
        oplot <- ordiplot(pc) # plots points only, plot is necessary for scaling vectors below
        # myplot <- plot(fit_env, p.max = 0.05) # plot significant environmental vectors
        vector.multiplier <- ordiArrowMul(scores(fit_env, choices=c(1,2), display="vectors"))
        fit_env_df <- as.data.frame(scores(fit_env, choices=c(1,2), display="vectors")*vector.multiplier)
#        fit_env_df <- cbind(fit_env_df, Variable = rownames(fit_env_df))
        fit_env_df <- cbind(fit_env_df, Variable = env.vars)  

        text_fit_env_df <- fit_env_df
        text_fit_env_df[,1] <- text_fit_env_df[,1]*.7 # let's shift the text out a little
        text_fit_env_df[,2] <- text_fit_env_df[,2]*1.3
                
        p <- p + coord_fixed() + ## need aspect ratio of 1!
        geom_segment(data = fit_env_df,
                   aes(x = 0, xend = Dim1, y = 0, yend = Dim2),
                   arrow = arrow(length = unit(0.25, "cm")), colour = "black") +
        geom_text(data = text_fit_env_df, aes(x = Dim1, y = Dim2, label = Variable),
                size = 4, fontface = "bold")
    }
    if(show.stats)
    {
        adonis.label <- paste0("1st-gen R2 = ", signif(adonis.years.R2, 2), ", P = ", signif(adonis.years.pval, 2), 
                                "\nGroup R2 = ", signif(adonis.sample.group.R2, 2), ", P = ", signif(adonis.sample.group.pval, 2))
        p <- ggdraw(p) + draw_figure_label(paste0(adonis.label), size=8, position="bottom.right")
    }
    if(save.pc)
    {
        write.table(pc, paste0(plot.title, "-PC.txt"), sep="\t", quote=F)
    }
    return(p)
}



# convex.hulls are shaded backgrounds around subject sample points
plot.pcoa.long <- function(map0, samples, otu0, method="euclidean", plot.title, dm=NULL, convex.hull=FALSE, flip.axis=NULL, num.clip.months=1, return.pc=FALSE)
{
    ret <- prep.dm(map0, otu0, dm, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0
    ppc <- cmdscale(ddm,2, eig=TRUE) # return eigenvals for calculating % var explained
    pc <- ppc$points

    if(!is.null(flip.axis))  pc[,flip.axis] <- -1*pc[,flip.axis] # sometimes we'll want to flip certain axes so orientations are the same
                    
    d <- data.frame(x = pc[,1], y = pc[,2], sample.id=rownames(pc), month=map0[rownames(pc),"Sample.Order"], subject=map0[rownames(pc),"Subject.ID"])

    mins <- NULL
    for(m in 0:(num.clip.months-1))
        mins <- rbind(mins, aggregate(d$month, list(d$subject), FUN=function(xx) min(xx)+m))
    colnames(mins) <- c("subject","month")            
    mins.d <- merge(mins, d, by=c("subject","month"))
    mins.d <- mins.d[!is.na(mins.d$month),]

    # find last x months, and average their response
    maxs <- NULL
    for(m in 0:(num.clip.months-1))
        maxs <- rbind(maxs, aggregate(d$month, list(d$subject), FUN=function(xx) max(xx)-m))
    colnames(maxs) <- c("subject","month")
    maxs.d <- merge(maxs, d, by=c("subject","month"))
    maxs.d <- maxs.d[!is.na(maxs.d$month),]

    # mark max and min samples appropriately
    d$Type <- "background"
    d[which(d$sample.id %in% mins.d$sample.id),"Type"] <- "start"
    d[which(d$sample.id %in% maxs.d$sample.id),"Type"] <- "end"
    d$Type <- factor(d$Type, levels=c("start","end"))

    p <- ggplot(data=d, aes(x, y)) + geom_point(color=alpha("gray",.5)) +
        xlab("PC1") + ylab("PC2")
    p <- p + geom_line(data=d[d$Type %in% c("start","end"),], aes(x, y, group=subject), color=alpha("black",.8)) + 
            geom_point(data=d[d$Type %in% c("start","end"),], aes(x, y, shape=Type, color=Type), fill="white", size=3) + 
            scale_shape_manual(name="legend", values=c(start=21,end=16)) + scale_color_manual(name="legend",values=c(start="black",end=alpha("#0B0BFD",.9))) 

    if(return.pc)
        return(d)
    else
        return(p)
}














#     if(length(cols) == 1) # single subject
#     {    
#         map.000 <- map0[samples,]
#         map.000$Sample.Day.Since.First.Sample <- as.numeric(as.Date(map.000$Sample.Date, format="%m/%d/%y") - as.Date("08/01/16", format="%m/%d/%y")) # hard code first sample date
#         map.000[,"travel.phase"] <- "Traveling"
#         map.000[map.000$Sample.Day.Since.First.Sample <= 2,"travel.phase"] <- "Pre"
#         map.000[map.000$Sample.Day.Since.First.Sample >= 28,"travel.phase"] <- "Post"
#         map.000$travel.phase <- factor(map.000$travel.phase, levels=c("Pre", "Traveling", "Post")) 
#         d2 <- data.frame(d[samples,], group2=map.000$travel.phase)
# 
#         p <- ggplot(data=d, aes(x, y)) +
#             geom_point(data = background.df, colour=alpha('gray',.2), shape=1, size=2) +
#             geom_point(data = d2, aes(colour=group2, shape=19), size=2) +
#             scale_shape_identity() +
#             xlab("PC1") + ylab("PC2") +
#             ggtitle(plot.title) + scale_color_manual(values=alpha(c("red","gray","blue"),.7)) 
# 
#     }
#     else{
