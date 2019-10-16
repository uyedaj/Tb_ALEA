## Figure 2 script
setwd("~/repos/Tb_ALEA/R/")
require(bayou)
require(OUwie)
require(treeplyr)
require(phytools)
#source("./R/utilityfunctions.R")
Tbsim <- read.csv("../data/trait_data_old")
#bt1 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_1.txt",header=TRUE, sep="\t")
#bt2 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_2.txt",header=TRUE, sep="\t")
#bt3 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_3.txt",header=TRUE, sep="\t")
#pairs(cbind(bt1$Mean.Scalar, bt2$Mean.Scalar, bt3$Mean.Scalar))
#tree1 <- read.nexus("../data/tree_old")
tree <- read.tree("../data/finaltree_BBMV")
tree <- reorder(tree, 'postorder')
#plot(tree$edge.length, bt1$Original.BL)
clarke <- read.csv("../data/Clarke.csv")
clarke[clarke=="*"] <- NA
clarke$Species <- gsub(" ", "_", clarke$Species)


td <- make.treedata(tree, clarke) %>% dplyr::select(., Tb) %>% filter(., !is.na(Tb))
colnames(td$dat) <- "Tb"

#taxalist <- lapply(2:nrow(bt1), function(x) strsplit(as.character(bt1$Taxa.List[x]), ",")[[1]])
#bt1 <- bt1[-1,-22]
#rtree <- tree
#rtree$edge.length <- bt1$Original.BL*bt1$Mean.Scalar
#plot(ladderize(rtree), cex=0.75)
#plot(log(bt2$Original.BL), log(bt2$Mean.Scalar))

pars <- list()
pars$alpha <- 0 # Phylogenetic half-life of 30 million years
pars$sig2 <- 0.02 # dunno if this is reasonable...
pars$k <- 1 # Number of shifts
pars$ntheta <- 2 # 2 optima, the root state and the Mammalian optimum
pars$theta <- 30 # Optima in degrees Cshifts <- list()
shifts <- list()
shifts$sb <- 895
shifts$loc <- 0
shifts$t2 <- 2
pars <- c(pars, shifts)

td$phy$edge.length[td$phy$edge.length==0] <- 0.001
## Prior
prior <- make.prior(td$phy, dists=list(dalpha="fixed", dsig2="dhalfcauchy", dsb="dsb", dk="ddunif", dtheta="dnorm", dloc="dunif"),
                    param=list(dalpha="fixed", dsig2=list(scale=0.1), dk=list(min=0, max=902), dtheta=list(mean=35, sd=10)),
                    fixed=list(alpha=0)) 
attributes(prior)$splitmergepars <- "sig2"

## Get some starting parameters, this is what your input for your likelihood function will be:
startpar <- pars
startpar$alpha <- 0
startpar$theta <- startpar$theta[1]
startpar$sig2 <- rep(startpar$sig2, startpar$ntheta)
#Test
prior(startpar)

## Fakish likelihood function, replace with your working version:
## Precalculate some things:
cache <- bayou:::.prepare.ou.univariate(td$phy, td[['Tb']])
custom.lik <- function(pars, cache, X, model="Custom"){
  phy <- cache$phy
  ## Bayou's method of mapping shift locations onto branches and segments (vectorized version of simmap)
  map <- bayou:::.pars2map(pars, cache)
  ## Transform branch lengths by sig2
  transfBL <- map$segs*pars$sig2[map$theta]
  ## Replace branch lengths on the tree
  phy$edge.length <- unname(tapply(transfBL, map$branch, sum))
  ## Get residuals from the root value
  X.c <- X - pars$theta
  ## Make likelihood function, can be made a lot faster by precalculating the cache object
  likfn <- geiger:::bm.lik(phy, X.c)
  ## Give fixed parameters sig2=1 and a SE of 0
  lnL <- likfn(c(1, 0))
  return(list(loglik=lnL, theta=pars$theta, resid=X.c))
}
custom.lik(startpar, cache, td[['Tb']])$loglik

monitorFn = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "sig2Root", "theta", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$sig2[1], pars$theta, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## Make bayou model object
model.auteur <- list(moves = list(alpha="fixed", sig2=".vectorMultiplier", 
                                  theta=".slidingWindowProposal", slide=".slide2",
                                  k=".splitmergebd"
),
## Relative proposal frequencies
control.weights = list(alpha=0, sig2=20, theta=4, slide=2, k=10),
## Tuning parameters, note that you need 1 for each reversible jump parameter
D = list(alpha=1, sig2=3, k=c(1), theta=6, slide=1),
parorder = c("alpha","theta","k", "ntheta", "sig2"),
rjpars = c("sig2"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = monitorFn,
lik.fn = custom.lik)

closeAllConnections()
mymcmc <- bayou.makeMCMC(td$phy, td[['Tb']], SE=0, model=model.auteur, prior=prior, startpar=startpar, outname="clarke_r001", new.dir="../output/bayou/", plot.freq=NULL, ticker.freq=1000, samp = 100, perform.checks = FALSE)
mymcmc$run(2000000)
#mymcmc$run(90000)

chain <- mymcmc$load()
chain <- set.burnin(chain, 0.3)
postburn <- floor(0.3*length(chain$gen)):length(chain$gen)
plot(chain)

## Function to pull out the transformed branch lengths for MCMC sample:
pullBL <- function(i, chain){
  pars <- pull.pars(i, chain, model=model.auteur)
  map <- bayou:::.pars2map(pars, cache)
  transfBL <- map$segs*pars$sig2[map$theta]
  cache$phy$edge.length <- unname(tapply(transfBL, map$branch, sum))/cache$phy$edge.length
  return(cache$phy)
}

############
allBLs <- sapply(postburn, function(x) pullBL(x, chain)$edge.length)
postmedBLs <- apply(allBLs, 1, mean)
postmedBLs <- postmedBLs/mean(sapply(chain$sig2, function(x) x[1]))

## Get all branches with a posterior probability above a certain cutoff:
L <- Lposterior(chain, cache$phy)
sb <- which(L$pp > 0.3)
pp <- list(k=length(sb), ntheta=length(sb)+1, sb=sb, loc=rep(0, length(sb)), t2=2:(length(sb)+1))
tr <- pars2simmap(pp, cache$phy)
tr$tree$edge.length <- postmedBLs * tr$tree$edge.length
#pdf("~/Downloads/auteurRun.pdf", width=10, height=10)

tipEdges <- which(cache$phy$edge[,2]  <= length(cache$phy$tip.label))
nH <- nodeHeights(tr$tree)
tipRates <- nH[tipEdges,2]
tipRates <- setNames(tipRates, tr$tree$tip.label[tr$tree$edge[tipEdges,2]])

td$dat$tipRates <- tipRates[cache$phy$tip.label]
.td <- filter(td, !is.na(Tb))
#.td <- treeply(.td, function(x) drop.tip(.td$phy, c("Tachyglossus_aculeatus", "Zaglossus_bruijni", "Ornithorhynchus_anatinus")))

#pdf(paste("../output/bayou/pathwiserates_Tb", j, ".pdf", sep="_"))
plot(chain)
#####
plotRegimes(tr$tree, cex=0.4, label.offset=0.01, show.tip.label=FALSE)
tiplabels(pch=21, cex=(scale(td[['Tb']])+2)/4, bg="red", adj=c(0.5, 0.5))
#####
plot(.td[['tipRates']], .td[['Tb']], pch=21, bg="gray80", col="gray80", xlab="Pathwise Rates (from bayou/auteur model)", ylab="Tb", ylim=c(15,45), xlim=c(0,1200))
lm1 <- phylolm(Tb~tipRates, data=.td$dat, phy=.td$phy, model="BM")
lm2 <- lm(Tb~tipRates, data=.td$dat)
abline(lm1, col="blue")
abline(lm2, col="red")
summary(lm1)
summary(lm2)
#dev.off()

empirical <- list(chain=chain, pgls=lm1, ols=lm2, td=.td, mcmc=mymcmc)
saveRDS(empirical, file="../output/empirical_bayou.rds")
##########Simulated Data#############
library(foreach)
library(doParallel)
library(phylolm)
library(extraDistr)
registerDoParallel(cores=20)
res <- foreach(j=6:10) %dopar% {
#source("./R/utilityfunctions.R")
Tbsim <- read.csv(paste("../data/trait_data_Tb", j, sep="_"))
#bt1 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_1.txt",header=TRUE, sep="\t")
#bt2 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_2.txt",header=TRUE, sep="\t")
#bt3 <- read.table("../data/var_rates_bayes_traits/Tb_real_BayesTraitvarrates_3.txt",header=TRUE, sep="\t")
#pairs(cbind(bt1$Mean.Scalar, bt2$Mean.Scalar, bt3$Mean.Scalar))
#tree1 <- read.nexus("../data/tree_old")
tree <- read.tree(paste("../data/finaltree_BBMV", j, sep="_"))
tree <- reorder(tree, "postorder")
#plot(tree$edge.length, bt1$Original.BL)

td <- make.treedata(tree, Tbsim) 
colnames(td$dat) <- "Tb"
td <- dplyr::select(td, Tb) %>% filter(., !is.na(Tb))

#taxalist <- lapply(2:nrow(bt1), function(x) strsplit(as.character(bt1$Taxa.List[x]), ",")[[1]])
#bt1 <- bt1[-1,-22]
#rtree <- tree
#rtree$edge.length <- bt1$Original.BL*bt1$Mean.Scalar
#plot(ladderize(rtree), cex=0.75)
#plot(log(bt2$Original.BL), log(bt2$Mean.Scalar))

pars <- list()
pars$alpha <- 0 # Phylogenetic half-life of 30 million years
pars$sig2 <- 0.02 # dunno if this is reasonable...
pars$k <- 1 # Number of shifts
pars$ntheta <- 2 # 2 optima, the root state and the Mammalian optimum
pars$theta <- 30 # Optima in degrees Cshifts <- list()
shifts <- list()
shifts$sb <- 895
shifts$loc <- 0
shifts$t2 <- 2
pars <- c(pars, shifts)

td$phy$edge.length[td$phy$edge.length==0] <- 0.001
## Prior
prior <- make.prior(td$phy, dists=list(dalpha="fixed", dsig2="dhalfcauchy", dsb="dsb", dk="ddunif", dtheta="dnorm", dloc="dunif"),
                    param=list(dalpha="fixed", dsig2=list(scale=0.1), dk=list(min=0, max=1128), dtheta=list(mean=35, sd=10)),
                    fixed=list(alpha=0)) 
attributes(prior)$splitmergepars <- "sig2"

## Get some starting parameters, this is what your input for your likelihood function will be:
startpar <- pars
startpar$alpha <- 0
startpar$theta <- startpar$theta[1]
startpar$sig2 <- rep(startpar$sig2, startpar$ntheta)
#Test
prior(startpar)

## Precalculate some things:
cache <- bayou:::.prepare.ou.univariate(td$phy, td[['Tb']])
custom.lik <- function(pars, cache, X, model="Custom"){
  phy <- cache$phy
  ## Bayou's method of mapping shift locations onto branches and segments (vectorized version of simmap)
  map <- bayou:::.pars2map(pars, cache)
  ## Transform branch lengths by sig2
  transfBL <- map$segs*pars$sig2[map$theta]
  ## Replace branch lengths on the tree
  phy$edge.length <- unname(tapply(transfBL, map$branch, sum))
  ## Get residuals from the root value
  X.c <- X - pars$theta
  ## Make likelihood function, can be made a lot faster by precalculating the cache object
  likfn <- geiger:::bm.lik(phy, X.c)
  ## Give fixed parameters sig2=1 and a SE of 0
  lnL <- likfn(c(1, 0))
  return(list(loglik=lnL, theta=pars$theta, resid=X.c))
}
custom.lik(startpar, cache, td[['Tb']])$loglik

monitorFn = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "sig2Root", "theta", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$sig2[1], pars$theta, pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

## Make bayou model object
model.auteur <- list(moves = list(alpha="fixed", sig2=".vectorMultiplier", 
                                  theta=".slidingWindowProposal", slide=".slide2",
                                  k=".splitmergebd"
),
## Relative proposal frequencies
control.weights = list(alpha=0, sig2=20, theta=4, slide=2, k=10),
## Tuning parameters, note that you need 1 for each reversible jump parameter
D = list(alpha=1, sig2=3, k=c(1), theta=6, slide=1),
parorder = c("alpha","theta","k", "ntheta", "sig2"),
rjpars = c("sig2"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = monitorFn,
lik.fn = custom.lik)

closeAllConnections()
mymcmc <- bayou.makeMCMC(td$phy, td[['Tb']], SE=0, model=model.auteur, prior=prior, startpar=startpar, outname=paste("BBMVsim_bayou_", j, sep="_"), new.dir="../output/bayou/", plot.freq=NULL, ticker.freq=1000, samp = 100, perform.checks = FALSE)
mymcmc$run(990000)
#mymcmc$run(90000)

chain <- mymcmc$load()
chain <- set.burnin(chain, 0.3)
postburn <- floor(0.3*length(chain$gen)):length(chain$gen)
saveRDS(chain, file=paste("BBMVsim_chain_",j,".rds", sep=""))
saveRDS(mymcmc, file=paste("BBMVsim_mcmc_",j,".rds", sep=""))


## Function to pull out the transformed branch lengths for MCMC sample:
pullBL <- function(i, chain){
  pars <- pull.pars(i, chain, model=model.auteur)
  map <- bayou:::.pars2map(pars, cache)
  transfBL <- map$segs*pars$sig2[map$theta]
  cache$phy$edge.length <- unname(tapply(transfBL, map$branch, sum))/cache$phy$edge.length
  return(cache$phy)
}

## Get averaged transformed edge lengths:
allBLs <- sapply(postburn, function(x) pullBL(x, chain)$edge.length)
postmedBLs <- apply(allBLs, 1, mean)
postmedBLs <- postmedBLs/median(sapply(chain$sig2, function(x) x[1]))

## Get all branches with a posterior probability above a certain cutoff:
L <- Lposterior(chain, td$phy)
sb <- which(L$pp > 0.3)
pp <- list(k=length(sb), ntheta=length(sb)+1, sb=sb, loc=rep(0, length(sb)), t2=2:(length(sb)+1))
tr <- pars2simmap(pp, td$phy)
tr$tree$edge.length <- postmedBLs * tr$tree$edge.length
#pdf("~/Downloads/auteurRun.pdf", width=10, height=10)

tipEdges <- which(td$phy$edge[,2]  <= length(td$phy$tip.label))
nH <- nodeHeights(tr$tree)
tipRates <- nH[tipEdges,2]
tipRates <- setNames(tipRates, tr$tree$tip.label[tr$tree$edge[tipEdges,2]])

td$dat$tipRates <- tipRates[td$phy$tip.label]
.td <- filter(td, !is.na(Tb))
#.td <- treeply(.td, function(x) drop.tip(.td$phy, c("Tachyglossus_aculeatus", "Zaglossus_bruijni", "Ornithorhynchus_anatinus")))

pdf(paste("../output/bayou/pathwiserates_Tb", j, ".pdf", sep="_"))
plot(chain)
#####
plotRegimes(tr$tree, cex=0.4, label.offset=0.01, show.tip.label=FALSE)
tiplabels(pch=21, cex=(scale(td[['Tb']])+2)/4, bg="red", adj=c(0.5, 0.5))
#####
plot(.td[['tipRates']], .td[['Tb']], pch=21, bg="gray80", col="gray80", xlab="Pathwise Rates (from bayou/auteur model)", ylab="Tb", ylim=c(15,45))
lm1 <- phylolm(Tb~tipRates, data=.td$dat, phy=.td$phy, model="BM")
lm2 <- lm(Tb~tipRates, data=.td$dat)
abline(lm1, col="blue")
abline(lm2, col="red")
summary(lm1)
summary(lm2)
dev.off()
}





