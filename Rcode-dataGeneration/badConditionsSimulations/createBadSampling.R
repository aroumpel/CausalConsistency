library(pcalg)
library(RBGL)
library(igraph)
library(corpcor)
source('sourceDir.R')
sourceDir('utilityFunctions')
sourceDir('codeFCI')

load("datasets.RData")

nDags <- length(datasets)
perc <- 0.1
newData <- list()

for (i in 1:nDags) {
  d <- datasets[[i]]
  k <- sample(1:nDags,1)  #random dag from which we will draw data
  nSamples <- dim(d)[1]
  toReplace <- perc*nSamples
  d[1:toReplace,] <- datasets[[k]][1:toReplace,]
  newData[[i]] <- d
}
save(newData,file = "newData.RData")

rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)

for (i in 1:nDags) {
  data <- newData[[i]]
  nSamples <- dim(data)[1]
  nNodes <- dim(data)[2]
  corMatrix<-cor(data)
  
  # run FCI
  suffStat <- list(C = corMatrix, n = nSamples)
  
  print(i)
  pag <- fci(suffStat,gaussCItest,alpha=0.05,p=nNodes, rules = rules, maj.rule = TRUE, aggressivelyPreventCycles = TRUE)
  runFCIipml(pag,folder="resultsFCImajrule",doSepSet=TRUE,nNodes)
}