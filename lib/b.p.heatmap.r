# plots a bubble heatmap of all bacteroides and prevotella strains, filtered by various settings
library(tidyr)
library(car)

# per = all, group, person
b.p.heatmap <- function(mapfile, otufile="data/shotgun/humarine-plus/otutable.txt", covfile="data/shotgun/humarine-plus/bpgenomes-coverage-per-sample.txt", 
                                    rescale="standard", outputfn, max.features=NULL, min.coverage=NULL, per="all",
                                    show.samplenames=FALSE, baseheight=4, basewidth=15, sample.groups=NULL)
{
    ret <- load.data(mapfile, otufile, normalize=F)
    otu0 <- ret$otu
    map0 <- ret$map
    
    cov.all <- read.table(covfile, comment="", sep="\t", head=T, as.is=T)
    if(!is.null(sample.groups)){
        map0 <- map0[map0$Sample.Group %in% sample.groups,]
        otu0 <- otu0[rownames(map0),]
    }

    valid.samples <- intersect(rownames(otu0), cov.all$Sample.ID)
    otu0 <- otu0[valid.samples,]
    map0 <- map0[valid.samples,]
    cov.all <- cov.all[cov.all$Sample.ID %in% valid.samples,]
    valid.otus <- intersect(colnames(otu0), cov.all$Genome.ID)
    otu0 <- otu0[,valid.otus]
    cov.all <- cov.all[cov.all$Genome.ID %in% valid.otus,]

    # normalize coverage by total sample depth
    #cov.all.norm <- merge(cov.all, data.frame(Counts=rowSums(otu0)), by.x="Sample.ID", by.y=0)
    #cov.all$Coverage <- cov.all.norm$Coverage/cov.all.norm$Counts

    covdf <- spread(cov.all, Genome.ID, Coverage)
    rownames(covdf) <- covdf$Sample.ID
    covdf <- covdf[,!(colnames(covdf) =="Sample.ID")]
    
    # turn NAs into zeros
    covdf[is.na(covdf)] <- 0
        
    cov.all <- merge(cov.all, map0[,"Sample.Group",drop=F], by.x="Sample.ID", by.y=0)
    #print(head(cov.all))
    
    if(per=="all")
    {
        # show abundance heatmap by only the highest coverage genomes overall
        cov.overall <- aggregate(Coverage ~ Genome.ID, data=cov.all, mean)    
        cov.overall$Genus <- ""
        cov.overall$Genus[grep("Bacteroides",cov.overall$Genome.ID)] <- "Bacteroides"
        cov.overall$Genus[grep("Prevotella",cov.overall$Genome.ID)] <- "Prevotella"
    
        b <- cov.overall[cov.overall$Genus=="Bacteroides",]
        p <- cov.overall[cov.overall$Genus=="Prevotella",]    

        if(!is.null(max.features)){
            topb <- b[order(b$Coverage,decreasing=T),][1:max.features,]
            topp <- p[order(p$Coverage,decreasing=T),][1:max.features,]    
            top.bugs <- c(topb$Genome.ID, topp$Genome.ID)
        } else if (!is.null(min.coverage)) {
            topb <- b[b$Coverage > min.coverage.overall,]
            topp <- p[p$Coverage > min.coverage.overall,]
            top.bugs <- c(topb$Genome.ID, topp$Genome.ID)
        }        
    } 
    else
    {
            if(per=="person") grouping.var <- "Sample.ID"
            else grouping.var <- "Sample.Group"
            
             # filter by min cov per group
            # let's take a look at highest coverage bugs in each sample group
            cov.overall <- aggregate(as.formula(paste0("Coverage ~ Genome.ID + ",grouping.var)), data=cov.all, mean)    
            cov.overall$Genus <- ""
            cov.overall$Genus[grep("Bacteroides",cov.overall$Genome.ID)] <- "Bacteroides"
            cov.overall$Genus[grep("Prevotella",cov.overall$Genome.ID)] <- "Prevotella"

            b <- cov.overall[cov.overall$Genus=="Bacteroides",]
            p <- cov.overall[cov.overall$Genus=="Prevotella",]    
        # allow filtering by BOTH coverage and number of features
        if (!is.null(min.coverage)) {           
            b <- b[b$Coverage > min.coverage,]
            p <- p[p$Coverage > min.coverage,]
            top.bugs <- c(b$Genome.ID, p$Genome.ID)
            top.bugs <- unique(top.bugs)
            topbpcov <- rbind(b,p)
        }

         if (!is.null(max.features)) {
            # order within group
            topb.all <- b[order(b[,grouping.var], b$Coverage, decreasing=T),]
            topp.all <- p[order(p[,grouping.var], p$Coverage, decreasing=T),]
            # now grab just the top X bugs
            topb <- aggregate(as.formula(paste0("Genome.ID ~ ", grouping.var)), data=topb.all, FUN=function(xx) xx[1:max.features])
            topp <- aggregate(as.formula(paste0("Genome.ID ~ ", grouping.var)), data=topp.all, FUN=function(xx) xx[1:max.features])
            top.bugs <- c(unique(as.vector(topb[,-1])), unique(as.vector(topp[,-1])))
            top.bugs <- top.bugs[!is.na(top.bugs)]

            # this code block for printing coverage to file only
            # note, this saves everything after the first column to be saved as a matrix! convert it properly after
            topb <- data.frame(as.character(topb[,1]), as.character(topb[,-1]), stringsAsFactors=F)
            topp <- data.frame(as.character(topp[,1]), as.character(topp[,-1]), stringsAsFactors=F)
            colnames(topb) <- c(grouping.var,"Genome.ID")
            colnames(topp) <- c(grouping.var,"Genome.ID")
            topb.all[,grouping.var] <- as.character(topb.all[,grouping.var])
            topp.all[,grouping.var] <- as.character(topp.all[,grouping.var])
            
            topbpcov <- rbind(merge(topb, topb.all, by=c(grouping.var,"Genome.ID")),
                        merge(topp, topp.all, by=c(grouping.var,"Genome.ID")))
        }
        # how many strains per person based on these constraints:
        print(table(topbpcov[,c(grouping.var,"Genus")]))
        # print the top coverage and genomes based on grouping var
        write.table(topbpcov, gsub("pdf","txt",outputfn), sep="\t",quote=F)

    }

    # do some custom wrangling of genome ids
    topotu <- otu0[,top.bugs]
    topcovdf <- covdf[,top.bugs]
    otu0.labels <- top.bugs

    otu0.labels <- gsub("-_.*", "", gsub(".*Prevotella", "Prevotella", otu0.labels))
    otu0.labels <- gsub("-_.*", "", gsub(".*Bacteroides", "Bacteroides", otu0.labels))
    colnames(topotu) <- otu0.labels
    colnames(topcovdf) <- otu0.labels

    if(rescale=="none") # if we're skipping rescaling by feature, let's just plot the normalized counts directly
        topotu <- sweep(topotu, 1, rowSums(topotu), '/')
    
    make.bubble.heatmap(otu0=topotu, map0=map0, outputfn=outputfn, 
                             min.prev=0.01, sig.level=1, show.rownames=T, rescale=rescale, as.bubble=T, cov.df=topcovdf, 
                             filter.mode="none", show.samplenames=show.samplenames, baseheight=baseheight, basewidth=basewidth)

}

### ARCHIVED - these don't look as good

#         # let's take a look at highest coverage bugs in each sample group
#         cov.group <- aggregate(Coverage ~ Genome.ID + Sample.Group, data=cov.all, mean)    
#         cov.group50 <- cov.group[cov.group$Coverage >= 50,]
#         top.bugs.group <- unique(cov.group50$Genome.ID)
#     
#         # keep only bugs with 90% coverage or higher per person
#         top.bugs.all <- unique(cov.all[cov.all$Coverage > .9,"Genome.ID"])
#         # there appears to be 3 strains that are not in the otu0 table but in the coverage table
#         top.bugs.all[which(!(top.bugs.all %in% colnames(otu0)))]
# 
#         valid <- intersect(colnames(otu0), colnames(covdf))
#         make.bubble.heatmap(otu0=otu0[,top.bugs.group], map0=map0, outputfn="heatmap.bubble.BP.top-per-group.pdf", 
#                              min.prev=0.01, sig.level=1, show.rownames=T, rescale=T, as.bubble=T, cov.df=covdf, filter.mode="none")
#         # Try all strains
#         make.bubble.heatmap(otu0=otu0[,valid], map0=map0, outputfn="heatmap.bubble.BP.ALL.pdf", 
#                              min.prev=0.01, sig.level=1, show.rownames=T, rescale=T, as.bubble=T, cov.df=covdf, filter.mode="none")
# 
#         # try normalizing coverage by counts and showing that
#         cov.norm <- covdf[rownames(otu0), c(topb,topp)] / otu0[,c(topb,topp)]
#         cov.norm[is.infinite(as.matrix(cov.norm))] <- 0
#         # heatmap of coverage normalized by counts (note that this will not distinguish between  highcov/lowcounts and lowcov/highcounts
#         make.heatmap.traditional(otu0=cov.norm, map0=map0, outputfn="heatmap.trad.BP.Cov-by-Counts.pdf", 
#                                  min.prev=0.01, sig.level=1, show.rownames=T, rescale=F)


    ### can we figure out a way to exclude low coverage bugs at high abundance?
#         counts <- unlist(as.list(otu0[rownames(covdf),valid]))
#         coverage <- unlist(as.list(covdf[,valid]))
#         # plot only things that have counts or coverage > 0
#         coverage.p <- coverage[!(coverage == 0 | counts == 0)]
#         counts.p <- counts[!(coverage == 0 | counts == 0)]
# 
#         # this looks weird
#         ggdata <- data.frame(counts=log10(counts.p), coverage=log10(coverage.p))    
#         ggdata$coverage <- ggdata$coverage/sd(ggdata$coverage)
#         m <- lm(coverage ~ counts, data=ggdata)
#         qqPlot(m) 
#         hist(residuals(m))
#             # residuals are basically linear --> our data can be assumed to be normal
#         ggplot(ggdata, aes(x=counts, y=coverage)) + geom_point(alpha=.1) + geom_smooth(method="loess", se = TRUE, color='black', fill="yellow", level=.95)

