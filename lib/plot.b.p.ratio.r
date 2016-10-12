plot.b.p.ratio <- function(map, taxa)
{
    bacteroides <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides"
    prevotella <- "k__Bacteria;p__Bacteroidetes;c__Bacteroidia;o__Bacteroidales;f__Prevotellaceae;g__Prevotella"

    samples <- rownames(map)
    
    US.born <- rownames(map[is.na(map$Years.in.US),])
    
    plot(map[samples,"Years.in.US"], taxa[samples,bacteroides], col=alpha("red", 0.5),pch=19, xlab="Years in US", ylab="Relative Abundance", xlim=c(0,41))

    points(map[samples,"Years.in.US"], taxa[samples,prevotella], col=alpha("green", 0.5),pch=19)

    abline(v=41.5,lty=2,col="gray28")
    points(rep(42,length(US.born)), taxa[US.born,prevotella], col=alpha("green", 0.8),pch=19)
    points(rep(42,length(US.born)), taxa[US.born,bacteroides], col=alpha("red", 0.8),pch=19)
    
    legend(x="topright", pch=19, col=alpha(c("red","green"), 0.5), legend=c("Bacteroides","Prevotella"),bg="white")

    mtext(side=1,line=1,text="US-\nborn", at=42)
}