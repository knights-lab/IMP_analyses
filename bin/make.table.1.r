# generate participant characteristics

datadirname <- "denovo"
LIBDIR="/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis/lib/"
setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/analysis")
source("bin/load.r")

map_orig <- map

# reset the placeholders that were put in for Years in US for other groups
map[!(map$Sample.Group %in% c("Hmong1st","Karen1st")),"Years.in.US"] <- NA

# aggregate tobacco usage
map$Tobacco <- as.character(map$Tobacco.Use)
map[which(map$Tobacco == "BN + Tobacco Daily"), "Tobacco"] <- "Daily"
map[which(map$Tobacco == "BN + Tobacco Occasionally"), "Tobacco"] <- "Occasionally"
map[is.na(map$Tobacco),"Tobacco"] <- "Never"
map[grep("BN", map$Tobacco), "Tobacco"] <- "Never" # any BN users without Tobacco, mark these as never tobacco users
map$Tobacco <- as.factor(map$Tobacco)

map$Betelnut <- as.factor(ifelse(grepl("BN", map$Tobacco.Use), "Y", "N"))

# order the education
map$Highest.Education <- factor(map$Highest.Education, levels=c("NONE","ESL","LHS","HS","C","G"))

vars.factor <- c("BMI.Class", "Alcohol.Use", "Tobacco", "Betelnut", "Highest.Education", "Type.Birth.Location", "Medical.Assistance", "Public.Housing", "Children.Free.Lunch")
vars.numeric <- c("Age", "Waist.Height.Ratio", "Years.in.US")

subject.df <- as.data.frame(table(map[cs,"Sample.Group"]))
colnames(subject.df) <- c("Sample.Group", "N")
rownames(subject.df) <- subject.df$Sample.Group
subject.df <- subject.df[,-1,drop=F]

# format as "mean (min-max)"
pval <- NULL
for(vn in vars.numeric)
{
    pval[vn] <- summary(aov(map[cs,vn] ~ map[cs,"Sample.Group"]))[[1]][[1,"Pr(>F)"]]
    
    temp <- data.frame(aggregate(map[cs,vn], by=list(map[cs,"Sample.Group"]), mean)[,2], aggregate(map[cs,vn], by=list(map[cs,"Sample.Group"]), range)[,2])
    subject.df[,paste(c("Mean","Min","Max"), vn, sep=".")] <- temp
    subject.df[,vn] <- paste0(signif(temp[,1],2), " (", signif(temp[,2],2), "-", signif(temp[,3],2), ")")
}

vars.us.only <-  c("Medical.Assistance","Public.Housing","Children.Free.Lunch")

# for each level, format as "N (%)"
for(vf in vars.factor)
{
    this.map <- map[cs,c("Sample.Group",vf)]
    if(vf %in% vars.us.only)
        pval[vf] <- chisq.test(this.map[!(this.map$Sample.Group %in% c("HmongThai","KarenThai")),"Sample.Group"], this.map[!(this.map$Sample.Group %in% c("HmongThai","KarenThai")),vf], , simulate.p.value=TRUE)$p.value
    else
        pval[vf] <- chisq.test(this.map[,"Sample.Group"], this.map[,vf], simulate.p.value=TRUE)$p.value
    
    var.level.counts <- reshape(as.data.frame(table(this.map)),direction="wide",timevar=vf,idvar="Sample.Group")[,-1]
    colnames(var.level.counts) <- paste(vf, colnames(var.level.counts),sep=".")
    subject.df[,colnames(var.level.counts)] <- var.level.counts
    var.level.percents <- signif(subject.df[,colnames(var.level.counts)]/subject.df$N,3)*100
    colnames(var.level.percents) <- gsub("Freq", "Percent", colnames(var.level.percents))
    subject.df[,colnames(var.level.percents)] <- var.level.percents

    var.levels <- levels(map[,vf])
    x <- lapply(var.levels, function(var.level) paste0(subject.df[, paste(vf,"Freq",var.level,sep=".")], " (", subject.df[, paste(vf,"Percent",var.level,sep=".")], ")"))

    subject.df[,gsub("Freq\\.","",colnames(var.level.counts))] <- x
}
# final edits
subject.df[c("HmongThai","KarenThai","Hmong2nd","Control"),"Years.in.US"] <- "NA"
subject.df[c("HmongThai","KarenThai"), paste(vars.us.only,"Y",sep=".")] <- "NA"

final.colnames <- colnames(subject.df)[-grep("Mean|Max|Min|Percent|Freq",colnames(subject.df))]

# keep only 1 value for booleans
final.colnames <- final.colnames[-grep("\\.N$", final.colnames)]

# save pvalues into table
table1 <- subject.df[,final.colnames]
table1["P", unlist(lapply(names(pval), function(p) min(grep(p, final.colnames))))] <- signif(pval,2)
table1["P", is.na(table1["P",])] <- ""

write.table(t(table1), file="Table1.Participants.txt", sep="\t", quote=F)


