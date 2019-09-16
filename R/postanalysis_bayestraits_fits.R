setwd("~/repos/Tb_ALEA/")
require(bayou)
require(OUwie)
require(treeplyr)
require(phytools)
#source("./R/utilityfunctions.R")
Tbsim <- read.csv("output/trait_data")
bt1 <- read.table("output/Tb_BayesTraitvarrates_1.txt",header=TRUE, sep="\t")
bt2 <- read.table("output/Tb_BayesTraitvarrates_2.txt",header=TRUE, sep="\t")
bt3 <- read.table("output/Tb_BayesTraitvarrates_3.txt",header=TRUE, sep="\t")
pairs(cbind(bt1$Mean.Scalar, bt2$Mean.Scalar, bt3$Mean.Scalar))
tree <- read.nexus("output/tree")
#plot(tree$edge.length, bt1$Original.BL)
clarke <- read.csv("./data/Clarke.csv")
clarke[clarke=="*"] <- NA
clarke$Species <- gsub(" ", "_", clarke$Species)

taxalist <- lapply(2:nrow(bt1), function(x) strsplit(as.character(bt1$Taxa.List[x]), ",")[[1]])
bt1 <- bt1[-1,-22]
rtree <- tree
rtree$edge.length <- bt1$Original.BL*bt1$Mean.Scalar
plot(ladderize(rtree), cex=0.75)
plot(log(bt2$Original.BL), log(bt2$Mean.Scalar))


td <- make.treedata(rtree, Tbsim)
colnames(td$dat) <- "Tb"
nH <- nodeHeights(td$phy)
terminals <- which(td$phy$edge[,2] <= length(td$phy$tip.label))
td$dat$pathwiserates <- nH[terminals,2]
td$phy <- tree
td <- filter(td, !is.na(Tb))
plot(td[['pathwiserates']], td[['Tb']])
lm1 <- lm(td[['Tb']]~ td[['pathwiserates']])
lm2 <- phylolm::phylolm(Tb~pathwiserates, data=td$dat, phy=td$phy, model="BM")
abline(lm1)
abline(lm2)
