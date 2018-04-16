# useful to keep separate so we can use this for other maps without loading in otu data

        # remove all IMP.000 from this analyses
        map <- map[map$Subject.ID != "IMP.000",]

        # format important variables factor levels and/or remove levels (change to char)
        map$BMI.Class <- as.character(map$BMI.Class)
        map$BMI.Class[map$BMI.Class == "Normal"] <- "Lean" # replace "Normal" with "Lean"
        map$BMI.Class <- factor(map$BMI.Class, levels=c("Lean", "Overweight", "Obese")) 
        map$Subject.ID <- as.character(map$Subject.ID)
        map$Sample.Group <- factor(map$Sample.Group, levels=c("KarenThai","HmongThai","Karen1st","Hmong1st","Hmong2nd","Control"))

        # all single-timepoint samples across both countries
        cs <- rownames(map)[is.na(map$Sample.Order) | map$Sample.Order==1]
    
        firstgen_cs <- rownames(map)[rownames(map) %in% cs & map$Sample.Group %in% c("Hmong1st","Karen1st")]

        hmong_secondgen_cs <- rownames(map)[map$Sample.Group == "Hmong2nd"]
        hmong_firstgen_cs <- rownames(map)[map$Sample.Group == "Hmong1st"]
        karen_firstgen_cs <- rownames(map)[rownames(map) %in% cs & map$Sample.Group == "Karen1st"]
        karenthai <- rownames(map)[map$Sample.Group=="KarenThai"]
        hmongthai <- rownames(map)[map$Sample.Group=="HmongThai"]
        controls <- rownames(map)[map$Sample.Group=="Control"]
    
        # calculate % of life spent in the US column
        map[,"Fraction.Life.in.US"] <- map$Years.in.US/map$Age
        map[c(hmong_secondgen_cs,controls),"Fraction.Life.in.US"] <- 1.0
        map[c(karenthai,hmongthai),"Fraction.Life.in.US"] <- 0

        # let's reset Years.in.US so that 2ndGen == 50 Controls = 60 and Thai == 0
        map[map$Ethnicity=="Caucasian","Years.in.US"] <- 60
        map[hmong_secondgen_cs,"Years.in.US"] <- 50
        map[c(karenthai,hmongthai),"Years.in.US"] <- 0

        # additional columns for easy plotting
        map[, "Group"] <- rep("1st", nrow(map))
        map[c(karenthai,hmongthai), "Group"] <- "Pre"
        map[hmong_secondgen_cs, "Group"] <- "2nd"    
        map[controls, "Group"] <- "Control"    
        map$Group <- factor(map$Group, levels=c("Pre","1st","2nd"))