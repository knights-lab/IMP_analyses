
### Predict Years in US using different covariates 
#   1) Covariates only (Age, BMI, Waist, Breastfed, Ethnicity, Medical.Assistance, Public.Housing, Children.Free.Lunch, Highest.Education, Religion)
#   2) Covariates + MB    
#   3) Diet? Diet + MB?    
    
    mapcs <- map[cs,]
    
    # remove "IMP000" since food doesn't exist for it
    mapcs <- mapcs[mapcs$Subject.ID != "IMP.000",]
    
    vars <- c("Age", "BMI", "Waist.Height.Ratio", "Breastfed", "Ethnicity", "Medical.Assistance", "Public.Housing", "Children.Free.Lunch", "Highest.Education", "Religion")

    mapcs[,vars[-(1:3)]] <- apply(mapcs[,vars[-(1:3)]], 1:2, as.character)
    varset1 <- mapcs[,vars]
    varset1[is.na(varset1)] <- "NA"
    templist <- lapply(varset1[,-(1:3)], as.factor)
    
    # for variable sets here
    vars_only <- cbind(varset1[,1:3], do.call(cbind.data.frame, templist))

    samples <- rownames(mapcs)
    
    mb_only <- as.data.frame(taxa_clr_L7[samples,])
    
    vars_mb <- cbind(vars_only, mb_only)
    
    food_only <- food_otu_L3[samples,]
    
    # this takes extremely long -- save data from here and run this line on teraminx instead
    # preds <- compare.rf.vars(mapcs$Years.in.US, list(vars_only, mb_only, vars_mb, food_only), n=10)
    
    preds <- list()
    varlist <- list(vars_only, mb_only, vars_mb, food_only)
    for(i in 1:length(varlist))
    {
        preds[[i]] <- compare.rf.vars(mapcs$Years.in.US, varsdf=varlist[[i]], n=1)
    }
    
    preds.df <- do.call(cbind.data.frame, preds)
    colnames(preds.df) <- c("Vars.Only", "MB.Only", "Vars.MB", "Food.Only")
    preds.df <- cbind(Sample.ID = samples, Actual = mapcs$Years.in.US, preds.df)

    ggdata <- melt(preds.df, id.vars=c("Actual","Sample.ID"))
    
    pre <- rownames(mapcs[mapcs$Years.in.US == 0,])
    secondgen <- rownames(mapcs[mapcs$Years.in.US == 50,])
    firstgen <- rownames(mapcs[!(mapcs$Years.in.US %in% c(0,50)),])
    
    cols <- c("#5fa55a", "#01b4bc", "#fa8925", "#fa5457")
    cols2 <- alpha(c("#5fa55a", "#01b4bc", "#fa8925", "#fa5457"), .5)
    names(cols) <- c("Vars.Only", "MB.Only", "Vars.MB", "Food.Only")
    names(cols2) <- c("Vars.Only", "MB.Only", "Vars.MB", "Food.Only")
        
    p.first <- ggplot(ggdata[ggdata$Sample.ID %in% firstgen,], aes(x = Actual, y = value, group=variable, color = variable)) +  
           geom_point() + geom_line() + scale_color_manual(values=cols2) + theme(legend.position='none') + xlab("Actual") + ylab("Predicted") + ggtitle("1st-Gen")
               
    p.pre <- ggplot(ggdata[ggdata$Sample.ID %in% pre,], aes(x = variable, y = value, group=variable, color = variable)) +  
            geom_boxplot(aes(fill=variable), colour="black", width=.9) + 
            scale_fill_manual(values = cols) + 
            ggtitle("Pre") + ylab("Predicted") + xlab("") + 
            theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), strip.background = element_blank(),
            strip.text = element_blank()) + coord_flip() + geom_hline(yintercept = 0, color="red")  + theme(legend.position='none')


    p.second <- ggplot(ggdata[ggdata$Sample.ID %in% secondgen,], aes(x = variable, y = value, group=variable, color = variable)) +  
            geom_boxplot(aes(fill=variable), colour="black", width=.9) + 
            scale_fill_manual(values = cols) +
            ggtitle("2nd-Gen") + ylab("Predicted") + xlab("") + 
            theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), strip.background = element_blank(),
            strip.text = element_blank()) + coord_flip() + geom_hline(yintercept = 50, color="red")

    leg <- get_legend(p.second)
    p.second <- p.second + theme(legend.position='none')
    
    p.combined <- plot_grid(plot_grid(p.pre, p.first, p.second, nrow=1, rel_widths=c(1,3,1)), leg, nrow=1, rel_widths=c(4,1))
    
    save_plot("rf.vars.pdf", p.combined, base_aspect_ratio = 1.7)