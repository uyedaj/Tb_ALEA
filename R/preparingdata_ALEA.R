setwd("~/repos/Tb_ALEA/")
require(bayou)
require(OUwie)
require(treeplyr)
#source("./R/utilityfunctions.R")

## Read in a tree
birdtree <- read.tree("../data/bird.tre")
fristoe <- read.csv("../data/fristoe.csv")
fristoe$genspec <- gsub(" ", "_", fristoe$Species)
tdbird <- make.treedata(birdtree, fristoe)
tdbird <- mutate(tdbird, lnBMR = log(BMR..mlO2.hour.,10), Tb=Tb...C., lnMass=log(Mass..g.,10))
plot(tdbird[['lnBMR']]/tdbird[['lnMass']], tdbird[['Tb']])
f1 <- phylolm(lnBMR~lnMass, data=tdbird$dat, phy=tdbird$phy, model="OUfixedRoot")
r1 <- resid(f1)
f2 <- phylolm(lnBMR ~ lnM + I(lnM^2), data=td$dat, phy=td$phy, model="OUfixedRoot")
r2 <- resid(f2)
tmp1 <- phylolm(r1~Tb, data=tdbird$dat, phy=tdbird$phy, model="OUfixedRoot")
tmp2 <- phylolm(r2~Tb, data=td$dat, phy=td$phy, model="OUfixedRoot")

plot(td$dat$Tb, r2,  xlab="Tb", ylab="Mass-Corrected BMR",pch=21, bg="green", xlim=c(30,45))
points(tdbird$dat$Tb, r1, pch=21, bg="blue")
abline(tmp1, col="blue", lwd=3, lty=2)
abline(tmp2, col="green", lwd=3, lty=2)

clarke <- read.csv("../data//Clarke.csv")
clarke[clarke=="*"] <- NA
clarke$Species <- gsub(" ", "_", clarke$Species)
tree <- read.nexus("../data/FritzTree.tre")
td <- make.treedata(tree[[1]], clarke)
td <- filter(td, !is.na(Tb), !is.na(BMR.W.))
td <- mutate(td, lnBMR=log(BMR.W.), lnM=log(Bm.g.))
td$phy <- multi2di(td$phy)
picBMR <- pic(resid(phylolm(td[['lnBMR']]~td[['lnM']], phy=td$phy)), phy=td$phy)
emat <- data.frame(node=td$phy$edge[,1], el=td$phy$edge.length)
bsum <- group_by(emat, node) %>% summarize(., mean=mean(el), min=min(el))
picTb <-  pic(resid(phylolm(td[['Tb']]~1, phy=td$phy)), phy=td$phy)
tmp <- data.frame(apicBMR= log(abs(picBMR)+0.01), apicTb=log(abs(picTb)+0.01))
tmp <- cbind(tmp, bsum)
tmp <- filter(tmp, apicBMR > -10, apicTb>-10)
plot(tmp$apicBMR, tmp$apicTb, pch=21, bg="white",
     xlab="ln absolute contrast (BMR)", ylab="ln absolute contrast (Tb)")
lm1 <- (lm(apicTb ~ apicBMR, data = tmp))


summary(lm1)
abline(lm1)
lm2 <- lm(apicTb ~ apicBMR, data=filter(tmp, mean<10))
abline(lm2)
lm2 <- lm(apicTb ~ mean, data = tmp)
lm3 <- lm(apicBMR ~ mean, data = tmp)
lm4 <- lm(resid(lm2)~resid(lm3))

plot(resid(lm3), resid(lm2))
abline(lm4)

nH <- nodeHeights(td$phy)
internal <- which(td$phy$edge[,2] > length(td$phy$tip.label))
plot(nH[internal,1], picBMR[-1])
td <- filter(td, Tb<45)

plot(td$dat$TempMean.C,td$dat$Tb, xlab="Ambient Temp", ylab="Body Temp", pch=21, bg="white")
.td <- filter(td, !is.na(Mass.g), !is.na(Tb), !is.na(TempMean.C))
.td <- mutate(.td, logBMR=log(BasalMetRate.mL02hr))
library(phylolm)
tmp <- phylolm(Tb ~ log(Mass.g), data=.td$dat, phy=.td$phy, model="OUfixedRoot")
plot(.td[['TempMean.C']], resid(tmp), xlab="Ambient Temp", ylab="Body Temp (mass corrected)", pch=21, bg="white")

.td <- filter(.td, !is.na(logBMR), is.finite(logBMR), logBMR <5)
tmp2 <- phylolm(logBMR ~ log(Mass.g), data=.td$dat, phy=.td$phy, model="OUfixedRoot")
plot(.td$dat$Tb, resid(tmp2) , xlab="BMR", ylab="Body Temp (mass corrected)", pch=21, bg="white")




## OLD Tb analyses with the ISIS/SP360 database
## Load in data
#load("./data/Tb.data2014.dta")
Tb.data <- read.csv("./data/panTb2014.csv")
isis <- read.csv("./data/NORMALS.csv")
Tb.isis <- subset(isis, TEST=="Body Temperature:")
pantheria <- read.csv("./data/PantheriaReduced.csv",header=TRUE)
Tb.data <- compileByPriority(Tb.data, pantheria)
Tb.isis$genspec <- unname(sapply(gsub(" ", "_", tolower(Tb.isis$SPECIES)), function(x) simpleCap(x)))
isis.pan <- Tb.isis$genspec[which((Tb.isis$genspec %in% pantheria$genspec) & !(Tb.isis$genspec %in% Tb.data$genspec))]
pan2add <- pantheria[match(isis.pan, pantheria$genspec),]
isis2add <- select(do.call(rbind, lapply(pan2add$genspec, function(x) filter(Tb.isis, genspec==x))), ISIS_MEAN, ISIS_SD, ISIS_N)
isis2add <- mutate(isis2add, "ME"= ifelse(ISIS_N==1, NA, sqrt(ISIS_SD^2/ISIS_N)), "Testicles"=NA)
isis2add <- select(isis2add, Testicles, ISIS_MEAN, ME, ISIS_N)
colnames(isis2add) <- c("Testicles", "Tb", "ME", "N")
pan2add <- cbind(pan2add, isis2add)
pan2add <- pan2add[,colnames(Tb.data)]
pan2add$Testicles <- NA
head(pan2add)

head(Tb.isis)
overlapping <- Tb.isis$genspec[which(Tb.isis$genspec %in% Tb.data$genspec)]
isis.overlapping <- Tb.isis[Tb.isis$genspec %in% overlapping,]
isis.overlapping <- filter(isis.overlapping, !is.na(SPECIES))
isis.overlapping <- select(isis.overlapping, genspec, ISIS_MEAN, ISIS_SD, ISIS_N, ISIS_MIN, ISIS_MAX)
Tb.overlapping <- Tb.data[Tb.data$genspec %in% overlapping,]
Tb.overlapping <- filter(Tb.overlapping, !is.na(genspec))
Tb.overlapping <- select(Tb.overlapping, genspec, Tb, ME, N)
Tb.overlapping <- arrange(Tb.overlapping, genspec)
isis.overlapping <- isis.overlapping[!duplicated(isis.overlapping$genspec),]
isis.overlapping <- arrange(isis.overlapping, genspec)
comb.overlap <- cbind(Tb.overlapping, isis.overlapping)
head(comb.overlap)
plot(comb.overlap$Tb, comb.overlap$ISIS_MEAN)
curve(1*x, add=TRUE)
cor.test(comb.overlap$Tb, comb.overlap$ISIS_MEAN)
replaceWithIsis <- function(i, comb.table, original.table){
  spp <- comb.table$genspec[i]
  origin.row <- which(original.table$genspec==spp)
  if(comb.table$ISIS_N > ifelse(is.na(original.table[origin.row, 'N']), 0,original.table[origin.row, 'N'])  | is.na(original.table[origin.row, 'N'])){
    original.table[origin.row, 'Tb'] <- comb.table$ISIS_MEAN[i]
    original.table[origin.row, 'ME'] <- sqrt((comb.table$ISIS_SD[i])^2/comb.table$ISIS_N[i])
    original.table[origin.row, 'N'] <- comb.table$ISIS_N[i]
  }
  return(original.table[origin.row,])
}
Tb.replaced <- lapply(1:nrow(comb.overlap),function(x) replaceWithIsis(x, comb.overlap, Tb.data))
Tb.replaced <- do.call(rbind, Tb.replaced)
Tb.data[which(Tb.data$genspec %in% Tb.replaced$genspec),] <- Tb.replaced
Tb.data <- arrange(Tb.data, genspec)
head(Tb.data)
Tb.data <- rbind(Tb.data, pan2add)


tree <- td$phy
dat <- td$dat
dat$ME[dat$ME==0] <- NA
MEr <- dat$ME
MEr[is.na(MEr) & !is.na(dat$N)] <- sqrt(1/dat$N[is.na(MEr) & !is.na(dat$N)])
MEr[is.na(MEr)] <- median(dat$ME, na.rm=TRUE)
dat$MEr <- MEr+0.05

phenogram(tree, setNames(dat$Tb, tree$tip.label), fsize=0.3)

plot(tree, show.tip.label=FALSE, type="fan")
tiplabels(pch=21, col=rainbow(105)[round(scale(dat$Tb)*10-min(scale(dat$Tb)*10),0)] ,bg=rainbow(105)[round(scale(dat$Tb)*10-min(scale(dat$Tb)*10),0)], cex=(scale(dat$Tb)+0.1-min(scale(dat$Tb)))/5)
median(dat$ME, na.rm=TRUE)
oufit <- fitContinuous(multi2di(tree), setNames(dat$Tb, tree$tip.label), SE=setNames(dat$MEr, tree$tip.label), model="OU")
bmfit <- fitContinuous(multi2di(tree), setNames(dat$Tb, tree$tip.label), SE=setNames(dat$MEr, tree$tip.label), model="BM")
## Phylogenetic half-life of the dataset
log(2)/oufit$opt$alpha
oufit$opt$aic
bmfit$opt$aic

