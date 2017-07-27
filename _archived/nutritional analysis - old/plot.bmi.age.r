data <- read.table(file="data.txt", header=T, sep="\t", row=1)

data[data$Ethnicity == "Karenni", "Ethnicity"] <- "Karen" #code these basically the same

lookup <- c(16,17)
names(lookup) <- sort(unique(data$Ethnicity))
final.pch <- lookup[as.character(data$Ethnicity)] 

rbPal <- colorRampPalette(c('orange','blue'))
data$Col <- rbPal(10)[as.numeric(cut(data$Age,breaks = 10))]

# ylim=c(17, 33)
par(mfrow=c(1,3))
plot(data$Years.in.US, data$BMI, pch=final.pch, col=data$Col, ylab="BMI", xlab="Years in US", main="BMI")
lmfit <- lm(data$BMI~data$Years.in.US)
pval <- summary(lmfit)$coefficients[2,4]
abline(lmfit, col="red")
mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)

plot(data$Years.in.US, data$Waist, pch=final.pch, col=data$Col, ylim = c(25,42), ylab="Waist Circumference (in)", xlab="Years in US", main="Waist Cir.")
lmfit <- lm(data$Waist~data$Years.in.US)
pval <- summary(lmfit)$coefficients[2,4]
abline(lmfit, col="red")
mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)

plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '', main = '')
legend(col=c("black","black","white","white"), pch=lookup, x="topleft", legend=c("Hmong", "Karen", " ", "Age"), bty='n')

legend_image <- as.raster(rbPal(10)[10:1])
text(x=0, y = seq(min(data$Age),max(data$Age),l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, .06, .3, .2, .6)


#legend(col="black", pch=lookup, x="topright", legend=c("Hmong", "Karen"))

#legend_image <- as.raster(rbPal(10)[10:1])
#plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '', main = 'Age')
#text(x=0, y = seq(min(data$Age),max(data$Age),l=5), labels = seq(0,1,l=5))
#rasterImage(legend_image, 0, 0, .5, .5)

############ plot nutrient data

data <- read.table(file="data.txt", header=T, sep="\t", row=1)

fiber <- read.table(file="nutrients.txt", header=T, sep="\t", row=1)

fiber.data <- cbind(data[rownames(fiber),], fiber)

data <- fiber.data

data[data$Ethnicity == "Karenni", "Ethnicity"] <- "Karen" #code these basically the same

lookup <- c(16,17)
names(lookup) <- sort(unique(data$Ethnicity))
final.pch <- lookup[as.character(data$Ethnicity)] 

rbPal <- colorRampPalette(c('blue','red'))
data$Col <- rbPal(10)[as.numeric(cut(data$BMI,breaks = 10))]

pdf("nutrients_over_time.pdf", useDingbats=F)
par(mfrow=c(2,3))

plot(data$Years.in.US, data$Sugar, pch=final.pch, col=data$Col, ylab="Sugar(g) Consumption", xlab="Years in US", main="Sugar")
lmfit <- lm(data$Sugar~data$Years.in.US)
pval <- summary(lmfit)$coefficients[2,4]
abline(lmfit, col="red")
mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)

# plot(data$Years.in.US, data$Fiber, pch=final.pch, col=data$Col, ylab="Fiber(g) Consumption", xlab="Years in US", main="Fiber")
# lmfit <- lm(data$Fiber~data$Years.in.US)
# pval <- summary(lmfit)$coefficients[2,4]
# abline(lmfit, col="red")
# mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)
# 
# plot(data$Years.in.US, data$Protein, pch=final.pch, col=data$Col, ylab="Protein(g) Consumption", xlab="Years in US", main="Protein")
# lmfit <- lm(data$Protein~data$Years.in.US)
# pval <- summary(lmfit)$coefficients[2,4]
# abline(lmfit, col="red")
# mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)

plot(data$Years.in.US, data$Total.Fat, pch=final.pch, col=data$Col, ylab="Total Fat (% Calories)", xlab="Years in US", main="%Fat")
lmfit <- lm(data$Total.Fat~data$Years.in.US)
pval <- summary(lmfit)$coefficients[2,4]
abline(lmfit, col="red")
mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)

plot(data$Years.in.US, data$Calories, pch=final.pch, col=data$Col, ylab="Total Calories", xlab="Years in US",main="Calories")
lmfit <- lm(data$Calories~data$Years.in.US)
pval <- summary(lmfit)$coefficients[2,4]
abline(lmfit, col="red")
mtext(paste("P = ", format(round(pval, 3), nsmall = 3), sep=""), side=3)

plot(c(0,1),c(0,1), type = 'n', axes = F, xlab = '', ylab = '', main = '')
legend(col=c("black","black","white","white"), pch=lookup, x="topleft", legend=c("Hmong", "Karen", " ", "BMI"), bty='n')

legend_image <- as.raster(rbPal(10)[10:1])
text(x=0, y = seq(min(data$Age),max(data$Age),l=5), labels = seq(0,1,l=5))
rasterImage(legend_image, 0, .5, .2, .7)

dev.off()
