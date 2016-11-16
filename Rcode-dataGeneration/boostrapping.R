#Script that creates 100 bootstrap datasets from the original ones and runs the FCI algorithm

load("datasets.RData")

nData <- length(datasets)
nBoots <- 100

bootDataset <- list()

for (i in 1:nData) {
  boot <- list()
  curData <- datasets[[i]]
  nSamples <- dim(curData)[1]
  for (j in 1:nBoots) {
    newData <- curData[sample(nSamples,nSamples,replace = T),]
    boot[[j]] <- newData
  }
  bootDataset[[i]] <- boot
}
save(bootDataset,file="bootDataset.RData")

bootFCI <- list()
bootFCIoriginal <- list()
bootFCImajrule <- list()

rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)

library(pcalg)
library(RBGL)
library(igraph)
library(corpcor)
source('sourceDir.R')
source('writeBoot.R')
sourceDir('utilityFunctions')
sourceDir('codeFCI')

for (i in 1:nData) {
  print(i)
  curData <- bootDataset[[i]]

  maj <- list()
  
  for (j in 1:length(curData)) {
    data <- curData[[j]]
    nSamples <- dim(data)[1]
    nNodes <- dim(data)[2]
    
    corMatrix<-cor(data)
    suffStat <- list(C = corMatrix, n = nSamples)
    
    maj[[j]] <- fci(suffStat,gaussCItest,alpha=0.05,p=nNodes, rules = rules, maj.rule = TRUE, aggressivelyPreventCycles = TRUE, doPdsep = FALSE)
    writeBoot(maj[[j]],"resultsFCImajrule",j)
    
  }
  bootFCImajrule[[i]] <- maj
}

save(bootFCImajrule,file="bootFCImajrule.RData")