library('vegan')

"vector.end" <- function(v1,v2,frac=.5) v1 + (v2-v1) * frac


# plots an animation timepoint 
# subject: vector of subject IDs
# timept: vector of integer timepoints; must be 0, 1, 2, ...
# step: plot all timepoints up to step (inclusive)
#       if step is non-integer, plots a partial line toward timepoint ceiling(step)
# cols: per-subject colors
# min.alpha: minimum level for transparency 
# alpha.decay: amount that alpha decays per timepoint
# align.subjects: if TRUE, ensures that each subject has timepoints 0, 1, 2, ...
"plot.animation.frame" <- function(pc, subject, timept, step, cols,
        min.alpha=0,alpha.decay=.5, align.subjects=TRUE,
        ellipse.alpha=.15){ ## Modified min alpha to be 0

    pc <- pc[,1:2]


    #plot(pc,type='n')
            
    subject <- factor(subject)

    if(align.subjects) {
        for(j in 1:length(levels(subject))){
            subj <- levels(subject)[j]
            ix <- subject == subj
            timept[ix] <- rank(timept[ix])-1
        }
    }

    # plot each timepoint successively, connected
    prev.ix <- NULL
    for(j in 1:length(levels(subject))){
        for(i in 0:floor(step)){
            subj <- levels(subject)[j]
            ix <- subject == subj & timept == i
            subj.col <- cols[ix]
            subj.col.alpha <- round(255 * max(min.alpha, 1 - alpha.decay * (step - i)))
            subj.col <- sprintf('%s%s',substr(subj.col,1,7), format(as.hexmode(subj.col.alpha),width=2))
            if(i==0) points(pc[ix,1],pc[ix,2],col=subj.col,pch=16,cex=1)
            # if(i==0) text(pc[ix,1],pc[ix,2],labels=paste(subj))
            if(i > 0 & any(ix)){
                prev.ix <- subject == subj & timept == (i-1)
                lines(pc[prev.ix | ix,1],pc[prev.ix | ix,2],col=subj.col,lwd=3)
            }
        }
    }

    frac <- step - floor(step)

    if(ceiling(step) > step){
        for(j in 1:length(levels(subject))){
            subj <- levels(subject)[j]
            ix1 <- subject == subj & timept == (ceiling(step) - 1)
            ix2 <- subject == subj & timept == ceiling(step)
            subj.col <- cols[ix1]
            if(any(ix2)){
                v <- vector.end(pc[ix1,], pc[ix2,],frac=frac)
                lines(c(pc[ix1,1],v[1]), c(pc[ix1,2], v[2]),col=subj.col, lwd=3)
            }
        }        
    }


    # draw transparent ellipse for each subject, only after past final timepoint
    # for that subject
    # for(j in 1:length(levels(subject))){
    #     subj <- levels(subject)[j]
    #     all.ix <- subject == subj & timept <= step
    #     pc.subj <- pc[all.ix,,drop=F]
    #     subj.col <- cols[which(all.ix)[1]]
    #     # fade out last ellipse, fade in new
    #     if(nrow(pc.subj) > 2 && !any(subject == subj & timept > step)){
    #         ellipse.alpha.j <- ellipse.alpha
    #         # fade in only at first
    #         if(any(subject == subj & timept == floor(step))){
    #             ellipse.alpha.j <- ellipse.alpha.j * frac
    #         }
    #         subj.col.alpha <- round(ellipse.alpha.j * 255)
    #         subj.col <- sprintf('%s%s',substr(subj.col,1,7), format(as.hexmode(subj.col.alpha),width=2))
    #         dataEllipse(pc.subj[,1], pc.subj[,2], plot.points=FALSE,center.cex=NA, levels=2 * pnorm(1) - 1,fill=TRUE,fill.alpha=ellipse.alpha.j,col=subj.col)
    #     }
    # }        
}
