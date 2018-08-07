library('ggpubr') #needed for ggplot stat_chull 
require(vegan)
require(faraway)
require(RColorBrewer)
require(cowplot)

# subset mapping file to contain exact samples prior to calling this function
plot.mouse.pcoa.donors <- function(map0, otu0, method="euclidean", outputfn, dm=NULL, flip.axis=NULL, axis1 = 1, axis2 = 2, verbose=FALSE)
{
    # add donor samples to the map so that we can generate pcoa with all necessary samples
    donor.map <- map_donor[map_donor$Donor %in% as.character(unique(map0$Donor)),]
    map0 <- map0[,c("Donor", "Group.End","Group.Start"),drop=F]
    map0 <- rbind(map0, donor.map)
    
    ret <- prep.dm(map0, otu0, dm, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0

    ppc <- cmdscale(ddm,max(axis1,axis2), eig=TRUE) # return eigenvals for calculating % var explained
    pc <- ppc$points

    # calculate variability explained
    percent_var <- signif(eigenvals(ppc)/sum(eigenvals(ppc)), 4)*100

    if(!is.null(flip.axis))  pc[,flip.axis] <- -1*pc[,flip.axis] # sometimes we'll want to flip certain axes so orientations are the same
                    
    d <- data.frame(x = pc[,axis1], y = pc[,axis2], row.names=rownames(pc), map0)
    print(d)
    # plot all mouse points without the donors first
    p <- ggplot(data=d[d$Group.Start != "Donor",], aes(x, y)) +
            geom_point(size=2, aes(color=Group.End, shape=Group.End, fill=Group.End, alpha=Group.End), stroke=1.25) +
            scale_color_manual(name="Groups", values=GROUP.COLORS, labels=GROUP.NAMES.SHORT) + 
            scale_fill_manual(name="Groups", values=GROUP.FILLS, labels=GROUP.NAMES.SHORT) + 
            scale_alpha_manual(name="Groups", values=GROUP.ALPHAS, labels=GROUP.NAMES.SHORT) +
            scale_shape_manual(name="Groups", values=GROUP.SHAPES, labels=GROUP.NAMES.SHORT)
    # add donors separately    
    p <- p + geom_point(data=d[d$Group.Start == "Donor",], aes(color=Group.End), shape=8, size=3)
               
    p <- p + xlab(paste0("PC",axis1, " [",percent_var[axis1],"%]")) + ylab(paste0("PC",axis2, " [",percent_var[axis2],"%]"))
    if(verbose)
        p <- p + geom_text(aes(label=rownames(d)), check_overlap=TRUE, size=1)

    leg <- get_legend(p)
    p <- p + theme(legend.position="none")
    save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1 )
    save_plot(paste0("legend.",outputfn), leg, useDingbats=FALSE, base_aspect_ratio = 1)

}


# convex.hulls are shaded backgrounds around subject sample points
# type = "path", "endpoints", "allpoints"
plot.mouse.pcoa <- function(map0, otu0, method="euclidean", outputfn, dm=NULL, type="endpoints", flip.axis=NULL, axis1 = 1, axis2 = 2, verbose=FALSE)
{
    ret <- prep.dm(map0, otu0, dm, method)
    ddm <- ret$ddm
    dm <- ret$dm
    map0 <- ret$map0

    ppc <- cmdscale(ddm,max(axis1,axis2), eig=TRUE) # return eigenvals for calculating % var explained
    pc <- ppc$points

    # calculate variability explained
    percent_var <- signif(eigenvals(ppc)/sum(eigenvals(ppc)), 4)*100

    if(!is.null(flip.axis))  pc[,flip.axis] <- -1*pc[,flip.axis] # sometimes we'll want to flip certain axes so orientations are the same
                    
    d <- data.frame(x = pc[,axis1], y = pc[,axis2], row.names=rownames(pc), map0)
    
    # order by time for plotting path appropriately
    d <- d[order(d$Mouse.ID, d$Week),]

    if(type=="endpoints") # plot all points but highlight start and endpoints (meaning Week 2, 8, 10)
    {
        endpoints <- rbind(aggregate(Week ~ Mouse.ID, d, FUN=max),
                    aggregate(Week ~ Mouse.ID, d[d$Week != 0,], FUN=min))   # exclude Week 0
        d.endpoints <- merge(d, endpoints, by=c("Mouse.ID","Week"))
    
        p <- ggplot(data=d, aes(x, y)) +
            geom_point(alpha=.2, size=2, aes(color=Group.Start), shape=16) +
            scale_color_manual(name="Groups", values=GROUP.COLORS, labels=GROUP.NAMES.SHORT)

        p <- p + geom_line(data=d.endpoints, aes(x=x, y=y, group=Mouse.ID, color=Group.End), alpha=.4)

        p <- p + geom_point(d=d[d$timepoint == "endpoint",], size=2, aes(x=x, y=y, color=Group.End, alpha=Group.End, shape=Group.End, fill=Group.End), stroke=1.25) +
            scale_shape_manual(name="Groups", values=GROUP.SHAPES, labels=GROUP.NAMES.SHORT) +
            scale_alpha_manual(name="Groups", values=GROUP.ALPHAS, labels=GROUP.NAMES.SHORT) +
            scale_fill_manual(name="Groups", values=GROUP.FILLS, labels=GROUP.NAMES.SHORT) 

    } 
    else if (type=="path") # this is not very useful
    {
        # draw path between points
        p <- ggplot(data=d, aes(x, y)) +
        geom_path(aes(group=Mouse.ID, color=Group.End), alpha=.4, show.legend=F) + 
        #        stat_chull(aes(color=Group.Start), alpha=.4, geom="path",size=.5) +
        geom_point(data=d[d$timepoint == "endpoint",], size=2, aes(fill=Group.End, color=Group.End, shape=Group.End, alpha=Group.End), stroke=1.25) +
        scale_color_manual(name="Groups", values=GROUP.COLORS, labels=GROUP.NAMES.SHORT) + 
        scale_fill_manual(name="Groups", values=GROUP.FILLS, labels=GROUP.NAMES.SHORT) + 
        scale_alpha_manual(name="Groups", values=GROUP.ALPHAS, labels=GROUP.NAMES.SHORT) +
        scale_shape_manual(name="Groups", values=GROUP.SHAPES, labels=GROUP.NAMES.SHORT)
    }
    else if (type == "allpoints")
    {
        p <- ggplot(data=d, aes(x, y)) +
            geom_point(size=2, aes(color=Group.End, shape=Group.End, fill=Group.End, alpha=Group.End), stroke=1.25) +
            scale_color_manual(name="Groups", values=GROUP.COLORS, labels=GROUP.NAMES.SHORT) + 
            scale_fill_manual(name="Groups", values=GROUP.FILLS, labels=GROUP.NAMES.SHORT) + 
            scale_alpha_manual(name="Groups", values=GROUP.ALPHAS, labels=GROUP.NAMES.SHORT) +
            scale_shape_manual(name="Groups", values=GROUP.SHAPES, labels=GROUP.NAMES.SHORT)
  
    }
    p <- p + xlab(paste0("PC",axis1, " [",percent_var[axis1],"%]")) + ylab(paste0("PC",axis2, " [",percent_var[axis2],"%]"))
    if(verbose)
        p <- p + geom_text(aes(label=rownames(d)), check_overlap=TRUE, size=1)

    leg <- get_legend(p)
    p <- p + theme(legend.position="none")
    save_plot(outputfn, p, useDingbats=FALSE, base_aspect_ratio = 1 )
    save_plot(paste0("legend.",outputfn), leg, useDingbats=FALSE, base_aspect_ratio = 1)

}