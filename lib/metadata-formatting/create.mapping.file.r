create.mapping <- function()
{
    setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/metadata formatting")

    samples_metadata_fp <- "IMP metadata - Samples.txt"
    participants_metadata_fp <- "IMP metadata - Participants.txt"  

    samples <- read.table(samples_metadata_fp, sep="\t", header=T, as.is=T, stringsAsFactors=F)
    samples$Sample.ID <- as.character(samples$Sample.ID)
    samples$Subject.ID <- as.character(samples$Subject.ID)
    
    participants <- read.table(participants_metadata_fp, sep="\t", header=T, as.is=T, stringsAsFactors=F)
    participants$Subject.ID <- as.character(participants$Subject.ID)

    # merge sample and participant metadata
    map <- merge(samples, participants, by.x="Subject.ID", by.y="Subject.ID")

    # let's fix up some of the columns and make it QIIME formatted
    map <- map[,c(2,1,3:ncol(map))]
    colnames(map)[colnames(map)=="Sample.ID"] <- "#SampleID"
    colnames(map)[colnames(map)=="Birthdate..Year."] <- "Birth.Year"
    colnames(map)[colnames(map)=="Notes.x"] <- "Notes.Samples"
    colnames(map)[colnames(map)=="Notes.y"] <- "Notes.Participants"
    colnames(map)[colnames(map)=="Old.Participant.ID.x"] <- "Old.Participant.ID"
    colnames(map)[colnames(map)=="Recruitment.Date.x"] <- "Recruitment.Date"
    # remove redundant columns
    ix <- match(c("Old.Participant.ID.y","Recruitment.Date.y"), colnames(map))
    map <- map[,-ix]

    map$Age <- as.numeric(map$Age)
    map$Years.in.US <- as.numeric(map$Years.in.US)

    # store ordering for longitudinal samples
    map$Sample.Order <- map$Sample.Month
    imp000.dates <- map[map$Subject.ID == "IMP.000", "Sample.Date", drop=F] # rownames not set but defaults should still work
    imp000.dates <- imp000.dates[order(as.Date(imp000.dates[,"Sample.Date"], format="%m/%d/%Y")),,drop=F]
    imp000.dates <- cbind(imp000.dates, Sample.Order=1:nrow(imp000.dates))
    map[rownames(imp000.dates),"Sample.Order"] <- as.numeric(imp000.dates$Sample.Order)
    
    # let's interpret Karenni as Karen - they are essentially the same (n=2)
    map$Ethnicity[map$Ethnicity=="Karenni"] <- "Karen"

    map$Sample.Group <- NA
    map$Sample.Group[!is.na(map$Years.in.US) & map$Ethnicity == "Karen"] <- "Karen1st"
    map$Sample.Group[!is.na(map$Years.in.US) & map$Ethnicity == "Hmong"] <- "Hmong1st"
    map$Sample.Group[is.na(map$Years.in.US) & map$Recruitment.Location %in% c("Maela Camp")] <- "KarenThai"
    map$Sample.Group[is.na(map$Years.in.US) & map$Recruitment.Location %in% c("Khun Chiang Kian")] <- "HmongThai"
    map$Sample.Group[is.na(map$Years.in.US) & !(map$Recruitment.Location %in% c("Khun Chiang Kian","Maela Camp"))] <- "Hmong2nd"
    map$Sample.Group[map$Sub.Study=="HC"] <- "Control"
    map$Years.in.US[map$Sample.Group %in% c("Hmong2nd")] <- 0

    # add waist-to-height ratio
    map$Waist.Height.Ratio <- map$Waist/map$Height
        
    # remove NA's in the SampleID column (dropouts)
    map <- map[!is.na(map[,"#SampleID"]),]

    # do some checks on formatting of strings
    check <- c("Recruitment.Location", "Researcher", "Sub.Study", "Exclude", "Public.Housing", "Medical.Assistance", "Children.Free.Lunch", "Highest.Education", "Religion", "Type.Birth.Location", "Type.location.before.US", "Tobacco.Use", "Alcohol.Use", "Breastfed", "Years.Breastfed")
    unique.vals <- apply(map[, check], 2, unique)
    lapply(unique.vals, sort)

    write.table(map, "mapping.txt",sep="\t",quote=F,qmethod="double",row.names=F)
    
    print("After running this, always search for special characters: ' \" < > ?")
}

# this takes a participant-to-date file, nutrients-to-date, and samples file and merges them 
# final output should be a sampleID to nutrients 
map.supertrackerdates.to.samples <- function()
{
    setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/metadata formatting")
    
    #### !!! Make sure that dates are properly formatted to MM/DD/YY in both of these files before running!
    participants.date.fn <- "participants-to-date.txt"
    samples.metadata.fn <- "IMP metadata - Samples.txt"  
    
    participants.date <- read.table(participants.date.fn, sep="\t", header=T, stringsAsFactors=F, strip.white=T)
    samples <- read.table(samples.metadata.fn, sep="\t", header=T, as.is=T, strip.white=T, na.strings=c("NA",""))
    participants.date$Diet.ID <- as.character(participants.date$Diet.ID)
    samples$Subject.ID <- as.character(samples$Subject.ID)

    # immediately remove any samples that were never received
    samples <- samples[!is.na(samples$Sample.ID),]
    rownames(samples) <- samples$Sample.ID
    samples.diet <- samples
    samples.diet$Diet.ID <- rep(NA, nrow(samples.diet))

    # fill in CS Diet IDs and SuperTracker Dates
    samples.diet[samples.diet$Sub.Study %in% c("CS","HC"), "Diet.ID"] <- samples[samples$Sub.Study %in% c("CS","HC"), "Old.Participant.ID"]
    cs.diets <- merge(samples.diet, participants.date, by="Diet.ID")
    samples.diet$SuperTracker.DATE <- rep(NA, nrow(samples.diet))
    samples.diet[cs.diets$Sample.ID, "SuperTracker.DATE"] <- cs.diets$SuperTracker.DATE
        
    # L fecal samples and diet samples do NOT have the same IDs. Therefore we need to match on Old Participant ID AND Month Order
    l.participants.date <- participants.date[substring(participants.date$Diet.ID, 1,1) == "L",]
    l.participants.date$Old.Participant.ID <- substring(l.participants.date$Diet.ID, 1, 5)
    l.participants.date <- l.participants.date[-which(l.participants.date$Diet.ID == "L.011.M1"),] # drop L.011.M1 because it's considered CS now
    # create copy of Sample.Month to preserve it after merge 
    l.samples.temp <- samples[samples$Sub.Study == "L",]

    # manually override the T.FS.CS. with L. because these initial samples were recruited as FullStool samples and not L samples
    camp.ix <- grep("T\\.FS\\.CS", l.samples.temp$Old.Participant.ID)
    l.samples.temp[camp.ix,"Old.Participant.ID"] <- gsub("\\.M.*", "", l.samples.temp[camp.ix, "Sample.ID"])

    l.samples.temp$Sample.Month.Copy <- l.samples.temp$Sample.Month
    l.samples <- merge(l.participants.date, l.samples.temp, by.x=c("Old.Participant.ID", "Diet.Month"), all.x=T, by.y=c("Old.Participant.ID", "Sample.Month"))     
    # rename copy of Sample.Month back 
    c.ix <- which(colnames(l.samples)=="Sample.Month.Copy")
    colnames(l.samples)[c.ix] <- "Sample.Month"
    
    # NAs in l.samples are diets without matching samples. let's save these for later. 
    extra.diets <- l.samples[is.na(l.samples$Sample.ID),]
    l.samples <- l.samples[!is.na(l.samples$Sample.ID),]
    
    samples.diet$Diet.Date <- rep(NA, nrow(samples.diet))
    samples.diet$Diet.Month <- rep(NA, nrow(samples.diet))    
    
    # fill longitudinal data
    columnnames <- c("Diet.ID","SuperTracker.DATE","Diet.Date","Diet.Month")
    samples.diet[l.samples$Sample.ID, columnnames] <- 
        l.samples[, columnnames]
    
    # fill in diet date for cs samples (same as recruitment date)
    samples.diet[samples.diet$Sub.Study %in% c("HC","CS"), "Diet.Date"] <- samples.diet[samples.diet$Sub.Study %in% c("HC","CS"),"Recruitment.Date"]

    # manually add L.011.M1
    samples.diet["L.011.M1", "SuperTracker.DATE"] <- participants.date[which(participants.date$Diet.ID=="L.011.M1"),"SuperTracker.DATE"]
    samples.diet["L.011.M1", c("Diet.ID","Diet.Date")] <- c("L.011.M1",samples.diet["L.011.M1", "Recruitment.Date"])
    
    # add back in any additional diets without matching samples in case we want to do diet-only analyses later
    extra.subjects <- unique(samples.diet[samples.diet$Old.Participant.ID %in% unique(extra.diets$Old.Participant.ID),c("Old.Participant.ID","Subject.ID")])
    merged <- merge(extra.diets[,"Old.Participant.ID",drop=F], extra.subjects, by="Old.Participant.ID")
    extra.diets$Subject.ID <- merged$Subject.ID
    extra.diets$Sample.ID <- paste0(extra.diets$Old.Participant.ID, ".M", extra.diets$Diet.Month, ".NOSAMPLE") ## do this here so we don't have NA samples!
    samples.diet <- rbind(samples.diet, extra.diets[,colnames(samples.diet)])
    
    # do some checks here
#     samples.diet[is.na(samples.diet$Diet.Date),] # print anything that hasn't been assigned a diet id
#     x <- samples.diet[,c("Diet.ID", "Sample.ID")]
#     y <- participants.date
#     dim(x)
#     dim(y)
#     x$df <- "samples"
#     y$df <- "participants"
#     xy <- merge(x, y, by="Diet.ID", all.x=T, all.y=T)
#     write.table(xy, file="test.txt", sep="\t", quote=F, row.names=F) # take a look to see what's getting assigned

    # remove samples marked for exclusion that also have no diet IDs (any with diet IDs, keep in case just so we can do proper matches later)
    remove.ix <- which(!is.na(samples.diet$Exclude) & is.na(samples.diet$Diet.Date))
    samples.diet <- samples.diet[-remove.ix,]
    
    # remove extra columns we don't care about
    remove.cols <- which(colnames(samples.diet) %in% c("Recruitment.Date", "Sample.Date", "Recruitment.Location", "Researcher", "Sub.Study", "Exclude", "Notes"))
    samples.diet <- samples.diet[,-remove.cols]

    # NOTE we started with 692 MB samples - 22 excluded + 6 extra diets 
    # final samples.diet == 676 rows

    write.table(samples.diet, "sampleid-to-dietid.txt", sep="\t", row.names=F, quote=F, qmethod="double")

    # in excel, reformat date columns to be MM/DD/YY  

}


map.nutrient.data <- function()
{    
    setwd("/Users/pvangay/Dropbox/UMN/KnightsLab/IMP/ANALYSES/metadata formatting")
    
    map.fn <- "sampleid-to-dietid.txt"
    nutrients.date.fn <- "nutrients-to-date.txt"
    meals.date.fn <- "meals-to-date.txt"
    foodgroups.date.fn <- "foodgroups-to-date.txt"
    
    map <- read.table(map.fn, sep="\t", header=T, check.names=F, strip.white=T, as.is=T)
    nutrients.date <- read.table(nutrients.date.fn, sep="\t", header=T, check.names=F, strip.white=T, as.is=T)
    meals.date <- read.table(meals.date.fn, sep="\t", header=T, check.names=F, strip.white=T, as.is=T)
    foodgroups.date <- read.table(foodgroups.date.fn, sep="\t", header=T, check.names=F, strip.white=T, as.is=T)
    
    # remove any excluded samples from map to avoid duplicates!
    metadata <- read.table("mapping.txt", sep="\t", header=T, check.names=F, as.is=T, comment="",row=1)
    metadata <- metadata[metadata$Exclude != "Y",]
    map <- map[map$Sample.ID %in% rownames(metadata),]
    
    # now process all nutrients, meals, foodgroups 
    ix <- which(colnames(nutrients.date) == "Total Fat (% Calories Eaten )")
    colnames(nutrients.date)[ix] <- "% of Calories from Total Fat"
    
    # merge two -to-date files, fix some colnames 
    # note that some samples might be dropped after merge since samples are excluded
    nutrients <- merge(map, nutrients.date, by="SuperTracker.DATE")
    meals <- merge(map, meals.date, by="SuperTracker.DATE")
    foodgroups <- merge(map, foodgroups.date, by="SuperTracker.DATE")

    # check that any dates not assigned are just excluded samples    
#     missing.nutrients <- nutrients.date[which(!(nutrients.date$SuperTracker.DATE %in% final_nutrients$SuperTracker.DATE)),"SuperTracker.DATE"]
#     missing.meals <- meals.date[which(!(meals.date$SuperTracker.DATE %in% final_meals$SuperTracker.DATE)),"SuperTracker.DATE"]
#     missing.foodgroups <- foodgroups.date[which(!(foodgroups.date$SuperTracker.DATE %in% final_foodgroups$SuperTracker.DATE)),"SuperTracker.DATE"]
#     participants.date.fn <- "participants-to-date.txt"
#     participants.date <- read.table(participants.date.fn, sep="\t", header=T, stringsAsFactors=F, strip.white=T)
#     samples <- read.table(samples.metadata.fn, sep="\t", header=T, as.is=T, strip.white=T, na.strings=c("NA",""))
#     missing.ids <- participants.date[participants.date$SuperTracker.DATE %in% unique(c(missing.nutrients, missing.meals, missing.foodgroups)), "Diet.ID"]
#     print(samples[samples$Old.Participant.ID %in% missing.ids,]) # make sure these are all excluded 

    # reorder with sampleID first
    final_nutrients <- nutrients[,c(which(colnames(nutrients)=="Sample.ID"), which(colnames(nutrients)!="Sample.ID"))]
    final_foodgroups <- foodgroups[,c(which(colnames(foodgroups)=="Sample.ID"), which(colnames(foodgroups)!="Sample.ID"))]
        
    final_meals <- meals[,c(which(colnames(meals)=="Sample.ID"), which(colnames(meals)!="Sample.ID"))]
    # rename these columns to match ASA24 so that FoodTree generation will be easier
    colnames(final_meals)[which(colnames(final_meals)=="modcode")] <- "ModCode"
    colnames(final_meals)[which(colnames(final_meals)=="foodcode")] <- "FoodCode"
    colnames(final_meals)[which(colnames(final_meals)=="Food")] <- "Main.food.description"
    
    # ** note that sometimes R is wonky and doesn't write out stuff properly, might just have to rerun
    
    write.table(final_nutrients, "nutrients.txt", sep="\t", row.names=F, quote=F,qmethod="double")
    write.table(final_meals, "meals.txt", sep="\t", row.names=F, quote=F,qmethod="double")
    write.table(final_foodgroups, "foodgroups.txt", sep="\t", row.names=F, quote=F,qmethod="double")
 
   
}

