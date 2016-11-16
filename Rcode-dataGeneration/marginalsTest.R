library(pcalg)
library(corpcor)
library(RBGL)
library(igraph)

source('createMarginals.R')
source('sourceDir.R')
sourceDir('utilityFunctions')
sourceDir('codeFCI')

nDags <- 50;
nNodes <- 20;
numMargs <- 100;

#number of variables to be removed from the original dataset
nRem <- 2;
sel2 = NULL
possibleCombs <- choose(nNodes,nNodes-nRem)
for (i in 1:numMargs) {
  selected <- sample(nNodes,nNodes-nRem,replace=FALSE)
  sel2 = rbind(sel2,selected)
}
save(sel2, file = "sel2.RData")

folder <- "resultsFCImajrule"
marginalsFCI <- list()
for (i in 1:nDags) {
  print(i)
  marginalsFCI[[i]] <- createMarginals(datasets[[i]],nRem,numMargs,sel2)
}
save(marginalsFCI,file="resultsFCImajrule/marginalsFCI_rm2.RData")