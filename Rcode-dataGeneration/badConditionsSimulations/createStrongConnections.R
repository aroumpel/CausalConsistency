library(pcalg)
library(RBGL)
library(igraph)
library(corpcor)
source('sourceDir.R')
sourceDir('utilityFunctions')
sourceDir('codeFCI')

load("dags.RData")

newDags <- dags

nDags <- length(newDags)
numStrong <- 2

for (i in 1:nDags) {
  #find how many edges there are, and change some of them randomly
  dag <- newDags[[i]]
  nEdges <- length(dag@edgeData@data)
  
  toChange <- sample(1:nEdges,numStrong)
  for (j in 1:numStrong) {
    dag@edgeData@data[[toChange[j]]]$weight <- 0.99
  }
  
  newDags[[i]] <- dag
}
save(newDags,file = "newDags.RData")

file.remove("dags.txt")
file.remove("resultsFCImajrule/pags.txt")
file.remove("resultsFCImajrule/mags.txt")
file.remove("resultsFCImajrule/sepSets.txt")

nSamples<-1000
nNodes <- 20
datasets <- list()
rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)

for (i in 1:nDags) {
  dag <- newDags[[i]]
  
  dag.amat <- matrix(0,nNodes,nNodes)
  
  for (j in 1:length(dag@edgeL)) {
    to <- dag@edgeL[[j]]$edges
    dag.amat[j,to] <- 1
  }
  
  write.table(dag.amat, file="dags.txt", append=TRUE, row.names=FALSE, col.names=FALSE)
  
  data <- rmvDAG(nSamples,dag)
  datasets[[i]] <- data
  
  corMatrix<-cor(data)
  
  # run FCI
  suffStat <- list(C = corMatrix, n = nSamples)
  print(i)
  pag <- fci(suffStat,gaussCItest,alpha=0.05,p=nNodes, rules = rules, maj.rule = TRUE, aggressivelyPreventCycles = TRUE)
  runFCIipml(pag,folder="resultsFCImajrule",doSepSet=TRUE,nNodes)
}
save(datasets,file = "datasets.RData")