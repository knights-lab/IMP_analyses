# should only need to do this once
# format the mouse mapping file 

    datadir <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse"
    taxa_L2_fn <- paste(datadir,"taxatable_L2.txt",sep="/")
    mice_fn <- "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/mouse/All Mouse Data.txt"

    taxa <- t(read.table(file=taxa_L2_fn, sep="\t", header=T, as.is=T, quote="", row=1,comment=""))
    mice <- read.table(file=mice_fn, sep="\t", header=T, as.is=T, quote="")
    mice$Sample.ID <- paste0(mice$Mouse.ID, ".", gsub("/","\\.",mice$Date))

    missing.seq <- mice$Sample.ID[which(!(mice$Sample.ID %in% rownames(taxa)))] # in mice, not in taxa
    not.in.mice <- rownames(taxa)[!(rownames(taxa) %in% mice$Sample.ID)]        # in taxa, not in mice

    # these are off by 1 day, because pellets were collected 1 day before sac (date that all other data was collected)
    not.mice.df <- as.data.frame(matrix(unlist(strsplit(not.in.mice, split="\\.")), ncol=4, byrow=T), stringsAsFactors=F)
    colnames(not.mice.df) <- c("Mouse.ID","Month","Day","Year")
    not.mice.df[,-1] <- apply(not.mice.df[,-1],1:2,as.numeric)
    not.mice.df$Sac.Day <- not.mice.df$Day + 1
    not.mice.df$Sac.Sample.ID <- paste(not.mice.df$Mouse.ID, not.mice.df$Month, not.mice.df$Sac.Day, not.mice.df$Year, sep=".")
    not.mice.df$Seq.Sample.ID <- paste(not.mice.df$Mouse.ID, not.mice.df$Month, not.mice.df$Day, not.mice.df$Year, sep=".")

    # now check if new sample IDs are indeed sequenced
    if(length(intersect(not.mice.df$Sac.Sample.ID, missing.seq)) != length(not.mice.df$Sac.Sample.ID))
        print("ERROR: there are still some samples that are missing in taxa!") 

    # replace the Sac SampleIDs with the Seq Samples IDs
    mice <- merge(mice, not.mice.df[,c("Sac.Sample.ID","Seq.Sample.ID")], by.x="Sample.ID", by.y="Sac.Sample.ID", all.x=TRUE)
    mice[!is.na(mice$Seq.Sample.ID), "Sample.ID"] <- mice[!is.na(mice$Seq.Sample.ID), "Seq.Sample.ID"]
    mice <- mice[,-which(colnames(mice)=="Seq.Sample.ID")]
    
    # any sample IDs that remain that aren't sequenced means that they actually WERE NOT sequenced
    if(length(mice$Sample.ID[mice$Sample.ID %in% rownames(taxa)]) != length(rownames(taxa)))
        print("ERROR: some sequenced samples are still missing")
             
    map <- mice
    rownames(map) <- mice$Sample.ID
    
    map$Date <- as.Date(map$Date, format="%m/%d/%y")
    
    # get the start and end dates for the entire map, but exclude dead mice 
    map.dead <- map[map$Mouse.ID %in% c("M25","M26","M28"),]
    map.alive <- map[!(rownames(map) %in% rownames(map.dead)),]
    
    endpoint <- aggregate(Date ~ Mouse.ID, map.alive, FUN=max) # only count mice who didn't die in the middle
    baseline <- aggregate(Date ~ Mouse.ID, map, FUN=min)
    endpoint$Date <- as.Date(endpoint$Date)
    baseline$Date <- as.Date(baseline$Date)
    endpoint$Mouse.ID <- as.character(endpoint$Mouse.ID)
    baseline$Mouse.ID <- as.character(baseline$Mouse.ID)

    map <- merge(map, rbind(cbind(endpoint,timepoint="endpoint"), cbind(baseline, timepoint="baseline")), by=c("Mouse.ID","Date"), all.x=T)
    map$timepoint <- as.character(map$timepoint)
    map$timepoint[is.na(map$timepoint)] <- "middle"

    # transform dates into experiment week
    for(i in 1:nrow(baseline))
    {
      this.mouse <- baseline$Mouse.ID[i]
      map[map$Mouse.ID == this.mouse, "Week"] <- (as.numeric(map[map$Mouse.ID == this.mouse, "Date"] - baseline$Date[i]))/7
    }

    map$Group.Start <- paste(map$Donor.Type, map$Diet.Type, sep=".")

    map$Group.End <- paste0(map$Donor.Type, ".", map$Diet.Type, ifelse(map$Cohoused, ".Cohoused", ""))
                    
    rownames(map) <- map$Sample.ID
    map <- map[,colnames(map)!="Sample.ID"]
    write.table(map, file="mapping.txt", sep="\t", quote=F, row.names=T, col.names=NA)
