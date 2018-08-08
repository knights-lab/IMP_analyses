require(scales)

### AESTHETICS
    SAMPLE.GROUP.NAMES <- c("KarenThai", "HmongThai", "Karen1st",  "Hmong1st",  "Hmong2nd",  "Control" )
    SAMPLE.GROUP.NAMES.SHORT <- c("KT", "HT", "K1",  "H1",  "H2",  "C" )
    names(SAMPLE.GROUP.NAMES.SHORT) <- SAMPLE.GROUP.NAMES
    
    SAMPLE.GROUP.COLS <- c("#018571", "#c51b7d", "#018571", "#c51b7d", "#c51b7d",  "black") # original

    SAMPLE.GROUP.COLS.DISTINCT <- c("#A300FF", "#FF0000", "#00A696",  "#FE42AD", "#FBB400", "#2E1915") # purple, red, teal, orange, blue, black

    SAMPLE.GROUP.SHAPES <- c(16, 16, 21, 21, 17, 17)    
    names(SAMPLE.GROUP.SHAPES) <- SAMPLE.GROUP.NAMES

    #SAMPLE.GROUP.SIZES <- c(1.5, 1.5, 1.3, 1.3, 1.6, 1.6)    
    SAMPLE.GROUP.SIZES <- c(2, 2, 2, 2, 2, 2)    
    names(SAMPLE.GROUP.SIZES) <- SAMPLE.GROUP.NAMES

    SAMPLE.GROUP.ALPHAS <- c(.7, .7, .7, .7, .5, .5)    
    names(SAMPLE.GROUP.ALPHAS) <- SAMPLE.GROUP.NAMES

    COLORS.BMI <- c(Lean="#FDD023",Overweight="#A78FA9",Obese="#644356") # FFDD91  DEE2ED

    get.bmi.colors <- function(alpha.val=1)
    {
        cols <- alpha(COLORS.BMI, alpha.val)
        names(cols) <- names(COLORS.BMI)
        return(cols)
    }
    
    get.group.colors <- function(groups=NULL, alpha.val=1)
    {
        cols <- SAMPLE.GROUP.COLS
    
        if(alpha.val != 1) cols <- alpha(SAMPLE.GROUP.COLS, alpha.val)

        names(cols) <- SAMPLE.GROUP.NAMES

        if(!is.null(groups))
            return(cols[groups])
        else
            return(cols)
    }

    get.group.colors.distinct <- function(groups=NULL, alpha.val=1)
    {
        cols <- SAMPLE.GROUP.COLS.DISTINCT
    
        if(alpha.val != 1) cols <- alpha(SAMPLE.GROUP.COLS.DISTINCT, alpha.val)
    
        names(cols) <- SAMPLE.GROUP.NAMES

        if(!is.null(groups))
            return(cols[groups])
        else
            return(cols)
    }

    get.group.alphas <- function(groups=NULL)
    {
        if(!is.null(groups))
            return(SAMPLE.GROUP.ALPHAS[groups])
        else
            return(SAMPLE.GROUP.ALPHAS)
    }

    get.group.sizes <- function(groups=NULL)
    {
        if(!is.null(groups))
            return(SAMPLE.GROUP.SIZES[groups])
        else
            return(SAMPLE.GROUP.SIZES)
    }
 
    get.group.shapes <- function(groups=NULL)
    {
        if(!is.null(groups))
            return(SAMPLE.GROUP.SHAPES[groups])
        else
            return(SAMPLE.GROUP.SHAPES)
    }

    check.normality <- function(x)
    {
        qqnorm(x)
        qqline(x)
        print(ks.test(x,"pnorm",mean(x),sqrt(var(x))))
        print(shapiro.test(x))
    }

    plot.scatter <- function(x,y,group=NULL, xlab, ylab, main, do.stats=TRUE)
    {
        ggdata <- data.frame(x=x, y=y)
        if(!is.null(group))
            ggdata <- cbind(ggdata, group)

        p <- ggplot(ggdata, aes(x=x, y=y)) + geom_point(alpha=.3, color="black", size=1.5) +
            xlab(xlab) + ylab(ylab) + ggtitle(main) + theme(legend.position='none') + 
            geom_smooth(method='lm', se=F)
    
        if(do.stats)
        {
            if(shapiro.test(ggdata$y)$p.value >= 0.05)
                cortest <- cor.test(ggdata$x, ggdata$y)
            else
                cortest <- cor.test(ggdata$x, ggdata$y, method="spear", exact=F)

            lab <- paste0(names(cortest$estimate), "=",signif(cortest$estimate, 2), "\nP=", signif(cortest$p.value,2))
    
            p <- ggdraw(p) + draw_figure_label(label=lab, size=8, position="bottom.right")                
        }
        return(p)
    }

    # plots first column of d by second column (x, y), with pval
    # useful for looking at interactions of 3-variables (used for BP-BMI-Years plot)
    # z.symmetric = use 3 colors for mean=0 z-values
    plot.XY <- function(d, z.symmetric=TRUE)
    {
        orig_cols <- colnames(d)
        if(ncol(d) == 2)
        {
            colnames(d) <- c("x","y")
            names(orig_cols) <- colnames(d)

            cortest <- cor.test(d$x, d$y, method="spear")
            pval <- signif(cortest$p.value, 2)
            rho <- signif(cortest$estimate,2)
            p <- ggplot(d, aes(x, y)) + geom_point() + ylab(orig_cols["y"]) + xlab(orig_cols["x"])
            p <- ggdraw(p) +  draw_figure_label(label=paste0("rho=", rho, "\np=",pval), size=8, position="bottom.right")                
        } 
        else # with interaction term
        {

            colnames(d) <- c("a","b","z")
            names(orig_cols) <- colnames(d)
        
            s <- summary(lm(z ~ a * b, data=d))
            pvals <- signif(s$coefficients[c("a","b"),"Pr(>|t|)"],2)
            pval <- paste(orig_cols[names(pvals)], ", P=", pvals, collapse="\n", sep="")
            p <- ggplot(d, aes(x=a, y=b)) + geom_point(aes(fill=z), shape=21, color=alpha("black",.2), size=2) +
                 ylab(orig_cols["b"]) + xlab(orig_cols["a"])

            if(z.symmetric)
                 p <- p + scale_fill_gradient2(low="#00441B", mid="white", high="#40004B", name=orig_cols["z"]) 
            else 
                p <- p + scale_fill_gradient(low="white", high="#40004B", name=orig_cols["z"]) 
                
            p <- ggdraw(p) +  draw_figure_label(label=pval, size=8, position="bottom.right")                
            print(s)        
        }
        return(p)        
    }


    # human study function for plotting by groups
    map.boxplot <- function(y, Group, main, facet.var=NULL, alpha=.1, add.pval=FALSE, plot.legend.only=FALSE, ylab="", strip.text.size=5, y.size=10, 
                            x.size=9, show.x=TRUE, group.vars.df=NULL)
    {
        legend.labels <- SAMPLE.GROUP.NAMES.SHORT[levels(Group)]
        if(show.x)
            x.labels <- legend.labels
        else
            x.labels <- NULL
            
        p <- map.boxplot.base(y=y, Group=Group, main=main, facet.var=facet.var, alpha=alpha, add.pval=add.pval, plot.legend.only=plot.legend.only, ylab=ylab, 
                            strip.text.size=strip.text.size, y.size=y.size, x.size=x.size,
                            fills=NA, cols=get.group.colors(), alphas=get.group.alphas(), shapes=get.group.shapes(), legend.labels=legend.labels, xlabels=x.labels, group.vars.df=group.vars.df) 
        return(p)
    }

    # base function for plotting groups by boxplot
    # fills = boxplot fill color (points will not be colored)
    # group.vars.df = extra variables used to make groups
    map.boxplot.base <- function(y, Group, facet.var=NULL, 
                        alpha=.1, add.pval=FALSE, 
                        ylab="", main, 
                        strip.text.size=5, title.size=12, y.size=10, x.size=9, 
                        hide.box=FALSE, cols, alphas, shapes, legend.labels, xlabels, fills=NA,
                        plot.legend.only=FALSE, group.vars.df=NULL)
    {
        d <- data.frame(y=y, Group=Group)
        if(!is.null(group.vars.df)) d <- cbind(d, group.vars.df)
        if(!is.null(facet.var)) d$facet.var <- facet.var

        p <- ggplot(d, aes(x=Group, y=y, group=Group)) + 
            scale_fill_manual(name="Groups", values=fills, labels=legend.labels) + 
            scale_color_manual(name = "Groups", values = cols, labels=legend.labels) +
            scale_alpha_manual(name="Groups", values=alphas, labels=legend.labels) +
            scale_shape_manual(name="Groups", values=shapes, labels=legend.labels) +
            ggtitle(main) + ylab(ylab) + xlab("") + 
            scale_x_discrete(labels=xlabels) + #
            theme(axis.text.y = element_text(size=y.size), axis.text.x = element_text(size=x.size),
            axis.title.x=element_blank(), axis.title.y=element_text(size=11),
            strip.background = element_blank(),
            strip.text.x = element_text(size = strip.text.size), 
            plot.title = element_text(size=title.size)) 

            
        if(hide.box){ # show median horizontal lines only, useful for immune plots
            d.summary <- aggregate(y ~ Group, median, data=d)
            p <- p + geom_quasirandom(dodge.width=.1, aes(shape=Group, color=Group), size=2.5, stroke=1) +
             geom_crossbar(data=d.summary, aes(ymin = y, ymax = y), size=.5, col="black", width = .75)
        }
        else{
            if(is.na(fills))
                p <- p + geom_quasirandom(dodge.width=.75, aes(shape=Group, color=Group, alpha=Group), size=2, stroke=1) + geom_boxplot(colour="black", width=.9, fill=NA)
            else
                p <- p + geom_boxplot(colour="black", width=.9, aes(fill=Group, alpha=Group)) + geom_quasirandom(dodge.width=.75, aes(shape=Group), color="black", size=2, stroke=1)
        }
    
        if(plot.legend.only)
          p <- get_legend(p)
        else
        {
            p <- p + theme(legend.position="none")    
            if(!is.null(facet.var))
                p <- p + facet_grid(. ~ facet.var)
            if(add.pval){
                if(is.null(group.vars.df))
                    p <- add.pvals.to.plot(p=p, model=aov(d$y ~ d$Group), group.var="d$Group", alpha=alpha)
                else
                {
                    f <- as.formula(paste0("y ~ ", paste(colnames(group.vars.df), collapse="+")))
                    p <- add.pvals.to.plot(p=p, model=aov(f, data=d), group.var=colnames(group.vars.df), alpha=alpha, method="two-way")
                }
            }
        }
        return(p)
    }


    "shorten.taxonomy" <- function(ids,delim=';'){
	ids <- gsub('[kpcofgs]__','',ids)
	newids <- ids
	ids <- strsplit(ids,delim)
	for(i in seq_along(ids)){
		n <- length(ids[[i]])
		j <- n
		while(ids[[i]][j] == 'NA' || ids[[i]][j] == '') j <- j - 1 # NA instead of "Other" here
		newids[i] <- ids[[i]][j]
	}
	return(newids)
}

### STATISTICS

    # x can be continuous or factor (vector)
    # allows for control variables to be included (passed in as a data.frame)
    # fits a linear model and assumes data is parametric (or has been transformed to be)
    test.features.parametric <- function(otu, x, controls=NULL, sig.level, paired=FALSE)
    {
        base.df <- x
        f <- "feature ~ x"
        if(!is.null(controls))
        {
            f <- paste0(f, " + ", paste(colnames(controls), collapse=" + "))
            base.df <- data.frame(x, controls)
        }
        ff <- as.formula(f)

        pvals <- apply(otu, 2, function(feature) {
                c <- summary(lm(ff, data.frame(feature, base.df)))$coefficients; return(c[2,4]);})
    
        adj.pvals <- p.adjust(pvals, "fdr")
    
        diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

        list(features=diff.features, adj.pvals=adj.pvals, pvals=pvals)
    }

    # allows only for group comparisons, therefore x must be a factor
    test.features.nonparametric <- function(otu, x, sig.level, paired=FALSE)
    {
        stopifnot(is.factor(x))

        if(length(unique(x))==2)
        {
            pvals <- apply(otu, 2, function(feature) 
                (wilcox.test(feature~x, data.frame(feature=feature, x=x), paired=paired))$p.value)
        }
        else
        {
            pvals <- apply(otu, 2, function(feature) 
            (kruskal.test(feature~x, data.frame(feature=feature, x=x)))$p.value)
        }
    
        adj.pvals <- p.adjust(pvals, "fdr")
    
        diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

        list(features=diff.features, adj.pvals=adj.pvals, pvals=pvals)
    }

    # allows only for two-group comparisons, therefore x must be a factor
    # samples1 and samples2 are paired samples ordered by subject
    test.features.paired <- function(otu, samples1, samples2, sig.level, parametric=FALSE)
    {
        if(parametric) pvals <- apply(otu, 2, function(feature) t.test(x=feature[samples1], y=feature[samples2], paired=T)$p.value)
        else pvals <- apply(otu, 2, function(feature) wilcox.test(x=feature[samples1], y=feature[samples2], paired=T, exact=F)$p.value)
    
        adj.pvals <- p.adjust(pvals, "fdr")
    
        diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

        list(features=diff.features, adj.pvals=adj.pvals, pvals=pvals)
    }

    # one-sample test, tests if values are == 0 (useful alternative to paired test for differences)
    test.features.one.sample <- function(otu, sig.level, parametric=FALSE)
    {
        if(parametric) pvals <- apply(otu, 2, function(feature) t.test(x=feature)$p.value)
        else pvals <- apply(otu, 2, function(feature) wilcox.test(x=feature, exact=F)$p.value)
    
        adj.pvals <- p.adjust(pvals, "fdr")
    
        diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]

        list(features=diff.features, adj.pvals=adj.pvals, pvals=pvals)
    }

    # x = +2-level factor, y = continuous value
    # simply runs wilcoxon test, OR manually loops through all possible combinations and corrects p-values for 3+ level group factors
    # better than kruskal so we can see direction
    test.groups <- function(y, x, parametric=TRUE)
    {
        comparisons <- combn(levels(x), 2)
    
        pvals <- NULL
        pnames <- NULL
        for(i in 1:ncol(comparisons))
        {
            this.y <- y[x %in% comparisons[,i]]
            this.x <- x[x %in% comparisons[,i]]

            if(parametric) m <- t.test(this.y ~ this.x)
            else m <- wilcox.test(this.y ~ this.x)
            pvals <- c(pvals, m$p.value)
            pnames <- c(pnames, paste0(comparisons[1,i], "-", comparisons[2,i]))
        }
        names(pvals) <- pnames
        return(pvals)
    }

    # if there are pairwise comparisons that are significant, add them to the plot
    # e.g. p <- add.pvals.to.plot(p, model=aov(ggdata$y ~ ggdata$Group), group.var="ggdata$Group", alpha=.10)
    add.pvals.to.plot <- function(p, model, group.var, alpha=.10, textsize=1.8, stepinc=.08, method="one-way")
    {
        if(method=="one-way")
        {
            tuke <- TukeyHSD(model)[[group.var]][,"p adj"] 
            sig.pvals <- signif(tuke[tuke < alpha],2)

            # if there are significant comparisons, show them
            if(length(sig.pvals) > 0) 
            {
                if(length(sig.pvals)==1) # comparing only a single comparison
                    p <- p + geom_signif(comparisons = list(c(1,2)), annotation=as.character(sig.pvals), color="black", tip_length=.01, textsize=textsize, step_increase=stepinc)
                else
                {
                    comparisons <- strsplit(names(sig.pvals),"-")
                    p <- p + geom_signif(comparisons = comparisons, annotation=as.character(sig.pvals), color="black", tip_length=.01, textsize=textsize, step_increase=stepinc)
                }
            }
            else
                print(tuke)
        }
        else
        {
            # for factorial designs, best to just print out p-values
            a <- Anova(model, type = "III")
            pvals <- a[group.var,"Pr(>F)"]
            names(pvals) <- rownames(a[group.var,])
            sig.pvals <- signif(pvals[pvals < alpha],2)
            
            if(length(sig.pvals) > 0) 
                p <- ggdraw(p) + draw_figure_label(label=paste(paste0(names(sig.pvals), " P=", sig.pvals), collapse="; "), size=7, position="bottom.right")               
            else
                print(pvals)
        }
        
        p
    }

### DISTANCES

    # depending on what's available (map, otutable, and distance matrix), generates appropriate DM and PCs
    # using samples that are currently in map0
    # add.samples.dm: allows for additional samples to be included in dm that you don't care to keep in the mapping
    prep.dm <- function(map0, otu0=NULL, dm=NULL, method="euclidean", add.samples.dm=NULL)
    {
        if(!is.null(dm)) # if dm is passed in, use it
        {    
            valid_samples <- intersect(rownames(dm), rownames(map0))
            map0 <- map0[valid_samples,]
            if(!is.null(add.samples.dm)) 
                valid_samples <- c(valid_samples, add.samples.dm)
            dm <- dm[valid_samples, valid_samples]
            ddm <- as.dist(dm)
        }
        else # if no dm, generate distance from OTU table
        {
            valid_samples <- intersect(rownames(otu0), rownames(map0))
            map0 <- map0[valid_samples,]
            if(!is.null(add.samples.dm)) 
                valid_samples <- c(valid_samples, add.samples.dm)
            otu0 <- otu0[valid_samples,]        
            ddm <- vegdist(otu0, method=method)
            dm <- as.matrix(ddm) # for use in stats later
        }
        return(list(map0=map0, ddm=ddm, dm=dm))
    }

    get.relative.distance <- function(query_samples, ref_samples, dm)
    {
        rel.distance <- NULL
        for (i in 1:length(query_samples))
        {
            rel.distance[i] <- mean(as.numeric(dm[query_samples[i], ref_samples]))
        }
        names(rel.distance) <- query_samples
        return(rel.distance)
    }

    get.within.dist<-function(dm, samples)
    {
        within <- numeric(length(samples))
        for(i in 1:length(samples))
        {
            within[i] <- mean(as.numeric(dm[samples[i], samples[-i]]))
        }
        names(within) <- samples
        return(within)
    }

    get.between.dist<-function(dm, samples, ref_samples)
    {
        between <- numeric(length(samples))
        for(i in 1:length(samples))
        {
            between[i] <- mean(as.numeric(dm[samples[i], ref_samples]))
        }
        names(between) <- samples
        return(between)
    }

