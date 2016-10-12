# calculates the fraction of sequences per sample that successfully hit the Greengenes database at 99% using EMBALMER
# hits = embalmer_cap.b6f from running EMBALMER, these are the sequences that actually hit the ref database
# queries = filter_d5_headers2.txt - these are the original raw sequence headers (parsed from input file passed in with EMBALMER)
fraction.hits <- function(hits, queries)
{
    data_table <- table(hits[,1])
    queries_table <- table(queries[,1])
    fraction <- data_table/queries_table
    write.table(fraction,"fraction_hitting.txt",quote=F,row.name=F,col.name=F)
}

library(scales)
plot.fraction.hits <- function(valid_samples, map)
{
    fraction <- read.table("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/fraction_hitting_GG99.txt", row=1, head=F, quote="", sep=" ")
    imp_fraction <- fraction[valid_samples,1] # load in fraction of hits that hit GG 99 that was previously calculated
    cols<-rep("blue",nrow(map))
    cols[map$Years.in.US==0]<-"Red"
    pdf(file="fraction_hits.pdf",useDingbats=F)
    plot(map$Years.in.US, imp_fraction, ylab="Fraction of known sequences", xlab="Years in the US", col=alpha(cols, 0.5),pch=19)
    dev.off()
}