## Figure 2 script
setwd("~/repos/Tb_ALEA/R/")
require(bayou)
require(OUwie)
require(treeplyr)
require(phytools)
#source("./R/utilityfunctions.R")
Tbsim <- read.csv("../data/trait_data_old")
bt1 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_1.txt",header=TRUE, sep="\t")
bt2 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_2.txt",header=TRUE, sep="\t")
bt3 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_3.txt",header=TRUE, sep="\t")
pairs(cbind(bt1$Mean.Scalar, bt2$Mean.Scalar, bt3$Mean.Scalar))
tree <- read.nexus("../data/tree_old")
#tree <- read.tree("../data/finaltree_BBMV")
#plot(tree$edge.length, bt1$Original.BL)
clarke <- read.csv("../data/Clarke.csv")
clarke[clarke=="*"] <- NA
clarke$Species <- gsub(" ", "_", clarke$Species)

taxalist <- lapply(2:nrow(bt1), function(x) strsplit(as.character(bt1$Taxa.List[x]), ",")[[1]])
bt1 <- bt1[-1,-22]
rtree <- tree
rtree$edge.length <- bt1$Original.BL*bt1$Mean.Scalar
plot(ladderize(rtree), cex=0.75)
plot(log(bt2$Original.BL), log(bt2$Mean.Scalar))


td <- make.treedata(rtree, clarke) %>% dplyr::select(., Tb)
colnames(td$dat) <- "Tb"
nH <- nodeHeights(td$phy)
terminals <- which(td$phy$edge[,2] <= length(td$phy$tip.label))
td$dat$pathwiserates <- nH[terminals,2]
td$phy <- tree
td <- filter(td, !is.na(Tb))
plot(td[['pathwiserates']], td[['Tb']], ylim=c(15, 45), xlim=c(0,1200))
lm1 <- lm(td[['Tb']]~ td[['pathwiserates']])
lm2 <- phylolm::phylolm(Tb~pathwiserates, data=td$dat, phy=td$phy, model="BM")
abline(lm1)
abline(lm2, col="red")
real <- list(lm=list(ols=lm1, pgls=lm2), td=td)

empirical <- list(chain=bt1, pgls=lm2, ols=lm1, td=td)
saveRDS(empirical, file="../output/empirical_bayestraits.rds")
#########################
setwd("~/repos/Tb_ALEA/R/")
require(bayou)
require(OUwie)
require(treeplyr)
require(phytools)
#source("./R/utilityfunctions.R")
Tbsim <- read.csv("../data/trait_data_old")
bt1 <- read.table("../data/var_rates_bayes_traits/Tb_BayesTraitvarrates_1.txt",header=TRUE, sep="\t")
bt2 <- read.table("../data/var_rates_bayes_traits/Tb_BayesTraitvarrates_2.txt",header=TRUE, sep="\t")
bt3 <- read.table("../data/var_rates_bayes_traits/Tb_BayesTraitvarrates_3.txt",header=TRUE, sep="\t")
pairs(cbind(bt1$Mean.Scalar, bt2$Mean.Scalar, bt3$Mean.Scalar))
#tree <- read.nexus("../data/tree")
tree <- read.tree("../data/finaltree_BBMV")
#clarke <- read.csv("../data/Clarke.csv")
#clarke[clarke=="*"] <- NA
#clarke$Species <- gsub(" ", "_", clarke$Species)

taxalist <- lapply(2:nrow(bt1), function(x) strsplit(as.character(bt1$Taxa.List[x]), ",")[[1]])
bt1 <- bt1[-1,-22]
td <- make.treedata(tree, Tbsim)

plot(td$phy$edge.length[td$phy$edge.length!=0], bt1$Original.BL)

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
plot(td[['pathwiserates']], td[['Tb']], ylim=c(15, 45), xlim=c(0,1200))
lm1 <- lm(td[['Tb']]~ td[['pathwiserates']])
lm2 <- phylolm::phylolm(Tb~pathwiserates, data=td$dat, phy=td$phy, model="BM")
abline(lm1)
abline(lm2)
sim <- list(lm=list(ols=lm1, pgls=lm2), td=td)

bt_addRates <- function(bt1, td){
  taxalist <- lapply(2:nrow(bt1), function(x) strsplit(as.character(bt1$Taxa.List[x]), ",")[[1]])
  bt1 <- bt1[-1,-22]
  td$phy <- reorder(td$phy, "cladewise")
  rtree <- td$phy
  oldtree <- td$phy
  rtree$edge.length <- bt1$Original.BL*bt1$Mean.Scalar
  td <- make.treedata(rtree, data.frame(species=td$phy$tip.label,"Tb"=td$dat$Tb))
  nH <- nodeHeights(td$phy)
  terminals <- which(td$phy$edge[,2] <= length(td$phy$tip.label))
  td$dat$pathwiserates <- nH[terminals,2]
  td$phy <- oldtree
  td <- filter(td, !is.na(Tb))
  pgls <- phylolm::phylolm(Tb~pathwiserates, data=td$dat, phy=td$phy, model="BM")
  return(list(td=td, pgls=pgls))
}
bayou_addRates <- function(chain, td){
  chain <- set.burnin(chain, 0.3)
  postburn <- floor(0.3*length(chain$gen)):length(chain$gen)
  td$phy$edge.length[td$phy$edge.length==0] <- 0.00001
  cache <- bayou:::.prepare.ou.univariate(td$phy, td[['Tb']])
  pullBL <- function(i, chain){
    pars <- pull.pars(i, chain, model=model.auteur)
    map <- bayou:::.pars2map(pars, cache)
    transfBL <- map$segs*pars$sig2[map$theta]
    cache$phy$edge.length <- unname(tapply(transfBL, map$branch, sum))/cache$phy$edge.length
    return(cache$phy)
  }
  allBLs <- sapply(postburn, function(x) pullBL(x, chain)$edge.length)
  postmedBLs <- apply(allBLs, 1, mean)
  postmedBLs <- postmedBLs/mean(sapply(chain$sig2, function(x) x[1]))
  
  ## Get all branches with a posterior probability above a certain cutoff:
  L <- Lposterior(chain, td$phy)
  sb <- which(L$pp > 0.3)
  pp <- list(k=length(sb), ntheta=length(sb)+1, sb=sb, loc=rep(0, length(sb)), t2=2:(length(sb)+1))
  tr <- pars2simmap(pp, td$phy)
  tr$tree$edge.length <- postmedBLs * tr$tree$edge.length
  tipEdges <- which(td$phy$edge[,2]  <= length(td$phy$tip.label))
  nH <- nodeHeights(tr$tree)
  tipRates <- nH[tipEdges,2]
  tipRates <- setNames(tipRates, tr$tree$tip.label[tr$tree$edge[tipEdges,2]])
  
  td$dat$pathwiserates <- tipRates[td$phy$tip.label]
  pgls <- phylolm::phylolm(Tb~pathwiserates, data=td$dat, phy=td$phy, model="BM")
  return(list(td=td, pgls=pgls))
}

## Loading simulations back in:
registerDoParallel(cores=10)
res <- foreach(i=1:10) %dopar% {
  bt1 <- read.table(paste("../data/parsed_VarRates_", i, ".txt", sep=""),header=TRUE, sep="\t")
  tree <- read.tree(paste("../data/finaltree_BBMV_", i, sep=""))
  tree <- reorder(tree, "postorder")
  dat <- read.csv(paste("../data/trait_data_Tb_", i, sep=""))
  chain <- readRDS(paste("BBMVsim_chain_", i,".rds", sep=""))
  td <- make.treedata(tree, dat)
  colnames(td$dat) <- "Tb"
  td_bt <- bt_addRates(bt1, td)
  td_bayou <- bayou_addRates(chain, td)
  return(list(bt=td_bt, bayou=td_bayou))
}


empirical_bayestraits <- readRDS("../output/empirical_bayestraits.rds")
empirical_bayou <- readRDS("../output/empirical_bayou.rds")

#### Figure 2
pdf("../output/figure2.pdf", height=4, width=6)
par(mfrow=c(1,2), cex.lab=1.25, mgp=c(2,1,0), mar=c(5,4,3,0), cex.main=1.5)
col1 <- bayou:::makeTransparent("#2ecc71", 100)
col3 <- bayou:::makeTransparent("#3498db", 100)
col2 <- bayou:::makeTransparent("#27ae60", 250)
col4 <- bayou:::makeTransparent("#2980b9", 250)
plot(empirical_bayestraits$td[['pathwiserates']], empirical_bayestraits$td[['Tb']], xlab="Pathwise rate for Tb", ylab="Tb",
     ylim=c(15,45), xlim=c(0,1200), pch=21, bg=col1, col=col1, main="Empirical data")
points(empirical_bayou$td[['tipRates']], empirical_bayou$td[['Tb']], pch=21, bg=col3, col=col3)

abline(empirical_bayestraits$pgls, col=col2, lwd=2, lty=2)
abline(empirical_bayou$pgls, col=col4, lwd=2, lty=2)

par(mar=c(5,0, 3, 4))
root=25
bt_points <- lapply(res, function(x) data.frame(pathwiserates=x$bt$td[['pathwiserates']], Tb=x$bt$td[['Tb']]))
bayou_points <- lapply(res, function(x) data.frame(pathwiserates=x$bayou$td[['pathwiserates']], Tb=x$bayou$td[['Tb']]))
plot(do.call(rbind, bt_points), type="n", xlab="Pathwise rate for Tb", ylab="Tb",
     ylim=c(15,45), xlim=c(0,1200), pch=21, bg=col1, col=col1, main="Simulated Data", yaxt="n")
#points(do.call(rbind, bayou_points), pch=21, bg=col3, col=col3)
lapply(res[1:5], function(x) abline(x$bt$pgls, lwd=1, lty=2, col=col2))
lapply(res, function(x) abline(x$bayou$pgls, lwd=1, lty=2, col=col4))
abline(h=root, lty=2, lwd=3)
text(900,root-2, labels="Root Tb", cex=1)

random = seq(from=15, to=45, length =200)
bounds = c(min(random), max(random))

#reverse of thermoregulation curve of enzyme
x <- random + 273.15
k <- 12.617e-5
B0 <- 8
V6 <- (B0 *(exp(-2.8/k*(1/x - 1/311.15)))/(1+ exp(9/k*(1/311.90-1/x))))
#plot(V6, random, type="n", xlim=c(1, -4), ylim=c(15,45), ylab="Tb",
#     main="Macroevolutionary Landscape", yaxt="n", xaxt="n", xlab="Thermal Performance")
lines(V6*200, random, lwd=3, lty=3, col="gray80")
dev.off()

######Landscape
#Params


##########
tdnm <- real$td
tdnm <- filter(tdnm, !tdnm$phy$tip.label %in% c("Tachyglossus_aculeatus", "Zaglossus_bruijni", "Ornithorhynchus_anatinus"))
ll_fitnm<- BBMV::lnL_FPK(tdnm$phy, tdnm[['Tb']], Npts=50, a=NULL, b=NULL, c=NULL)
fitnm <- BBMV::find.mle_FPK(ll_fitnm)
#fitnm <- BBMV::MH_MCMC_FPK(tdnm$phy, tdnm[['Tb']], bounds=c(0,49))
asrnm <- ACE_FPK(fitnm)


ll_fit4 <- BBMV::lnL_FPK(real$td$phy, real$td[['Tb']], Npts=50, a=NULL, b=NULL, c=NULL)
fit4 <- BBMV::MH_MCMC_FPK(real$td$phy, real$td[['Tb']], bounds=c(0,49))
asr <- ACE_FPK(fit4)
mods <- lapply(c("OUrandomRoot", "OUfixedRoot", "BM"), function(x) phylolm(Tb~1, data=real$td$dat, phy=real$td$phy, measurement_error=TRUE, model=x))

height <- branching.times(real$td$phy)
height <- max(height)-height
height <- c(rep(166.2, 461), unname(height))
res <- asr
res <- lapply(res, function(x){y <- x; y[which(y[,2] <0.02),1] <- NA; y})
plot(0, max(height), type="n", xlab="Time", ylab="Tb", ylim=c(25,45), xlim=c(0, 200))
lapply(462:length(res), function(x) lines(height[x] + 50*res[[x]][,2],res[[x]][,1], col=col2))
points(rep(200, 461), real$td[['Tb']], pch="-", cex=2, col=col1)

tdcl <- paint_clades(real$td, 3)
noneuth <- filter(tdcl, clades %in% c(2,4))
euth <- filter(tdcl, clades==3)

ll_fitne<- BBMV::lnL_FPK(noneuth$phy, noneuth[['Tb']], Npts=50, a=NULL, b=NULL, c=NULL)
fitne <- BBMV::find.mle_FPK(ll_fitne)
asrne <- ACE_FPK(fitne)

ll_fite<- BBMV::lnL_FPK(euth$phy, euth[['Tb']], Npts=50, a=NULL, b=NULL, c=NULL)
fite <- BBMV::find.mle_FPK(ll_fite)
asre <- ACE_FPK(fite)

ll_fitsim <- lnL_FPK(sim$td$phy, sim$td[['Tb']], Npts=50, a=NULL, b=NULL, c=NULL)
fitsim <- find.mle_FPK(ll_fitsim)
asrsim <- ACE_FPK(fitsim)
