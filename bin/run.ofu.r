    oFu_fns <- c("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/OFU/ofuprofile_species_40.txt",
    "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/OFU/ofuprofile_species_80.txt",
    "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/OFU/ofuprofile_species_90.txt",
    "/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/data/OFU/ofuprofile_species_95.txt")

### plot differential OFUs by BMI class
    output.table <- NULL
    for(i in 1:length(oFu_list))
    {
        ofu <- oFu_list[[i]]
        ret <- plot.diff.taxa(lean.obese.cs.map, ofu, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend=paste0("BMI.all.", names(oFu_list)[i]), sig.level=.25)
    
        if(nrow(ret) > 0)    
            output.table <- rbind(output.table, cbind(ret, similarity=names(oFu_list)[i], group="All Groups"))
    }
    # let's look at 1st gen immigrants only
    submap <- lean.obese.cs.map[lean.obese.cs.map$Sample.Group %in% c("Karen1st", "Hmong1st"),]
    for(i in 1:length(oFu_list))
    {
        ofu <- oFu_list[[i]]
        ret <- plot.diff.taxa(submap, ofu, x.var="BMI.Class", 
            control.vars=c("Age","Years.in.US","Ethnicity"), outputfn.prepend=paste0("BMI.1stgen.", names(oFu_list)[i]), sig.level=.25)

        if(nrow(ret) > 0)
            output.table <- rbind(output.table, cbind(ret, similarity=names(oFu_list)[i], group="1stGen"))
    }
    write.table(output.table, file="IMP-OFU-lean-v-obese.txt", quote=F, row.names=F, sep="\t")
