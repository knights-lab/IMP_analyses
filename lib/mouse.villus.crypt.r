    calc.vc <- function(vc, type, return.summary=TRUE)
    {
        villus <- vc[vc$Type==type,]
        villus <- aggregate(Measurement ~ Mouse.ID + Group, villus, FUN=mean)
        v.mean <- aggregate(Measurement ~ Group, villus, FUN=mean)
        v.sd <- aggregate(Measurement ~ Group, villus, FUN=sd)
        v <- merge(v.mean, v.sd, by=1)
        colnames(v) <- c("Group", "Mean", "SD")
        if(return.summary) return(v)
        else return(villus)
    }   
     
    pairwise.lm <- function(groups, f, d)
    {
        pval <- NULL
        for(i in 1:length(groups))
        {
             pval <- c(pval, summary(lm(f, data=d[d$Group %in% groups[[i]],]))$coefficients[2,4])
             print(summary(lm(f, data=d[d$Group %in% groups[[i]],]))$coefficients)
        }
        return(p.adjust(pval, method="fdr"))
    }

    plot.vc <- function(villus_crypt)
    {
        compare.groups <- list(c("US.HighFiber","Thai.HighFiber"), c("US.LowFiber","Thai.LowFiber"), c("US.HighFiber","US.LowFiber"), c("Thai.HighFiber","Thai.LowFiber"))
        
        vl <- data.frame(calc.vc(villus_crypt, "villus (um)"), type="villus", stringsAsFactors=F)
        cd <- data.frame(calc.vc(villus_crypt, "crypt (um)"), type="crypt", stringsAsFactors=F)
        
        temp <- villus_crypt[villus_crypt$Type=="villus (um)", ]
        temp$Measurement <- villus_crypt[villus_crypt$Type=="villus (um)", "Measurement"] / villus_crypt[villus_crypt$Type=="crypt (um)", "Measurement"]
        temp$Type <- "VC Ratio"
        vc <- data.frame(calc.vc(temp, "VC Ratio", return.summary=FALSE), type="VC", stringsAsFactors=F)
        
        p2 <- ggplot(vc, aes(x=Group, y=Measurement, fill=Group)) + geom_boxplot(aes(fill=Group), colour="black", width=.9) +
                    scale_fill_manual(values=GROUP.COLORS) + ylab("") + ggtitle("Villus-Crypt Ratio") +
                    theme(axis.line.x = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
        save_plot("villus-crypt-ratio.pdf", p2, base_aspect_ratio = 1.3)

        adj.pval <- pairwise.lm(groups=compare.groups, f=as.formula("Measurement ~ Group"), d = aggregate(Measurement ~ Mouse.ID + Group, temp, FUN=mean))
        
        cd$Mean <- -1*cd$Mean
        d <- rbind(vl, cd)
        d$group2 <- paste0(d$Group, d$type)
    
        p <- ggplot(d, aes(x=Group, y=Mean, fill=Group, alpha=type)) + 
            geom_bar(stat="identity", position="identity") + geom_errorbar(aes(ymax=Mean+SD, ymin=Mean-SD), position = "identity", width = 0.25) +
            scale_fill_manual(values=GROUP.COLORS) + 
            scale_alpha_manual(name="", values=c(.5, .9)) +
            theme(axis.line.x = element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
            guides(alpha = guide_legend(order = 0), fill = guide_legend(order = 1)) +
            ylab("(um)") + xlab("") + ggtitle("Villus Height and Crypt Depth")
        p <- p + geom_hline(yintercept=0, color="black", size=.5)

        # play around with adding asterisks for significance    
        villus <- villus_crypt[villus_crypt$Type=="villus (um)",]
        villus <- aggregate(Measurement ~ Mouse.ID + Group, villus, FUN=mean)
        adj.pval <- pairwise.lm(groups=compare.groups, f=as.formula("Measurement ~ Group"), d = villus)

        #p <- p + stat_compare_means(comparisons = compare.groups)
    
        save_plot("villus-crypt.pdf", p, base_aspect_ratio = 1.3)
    }