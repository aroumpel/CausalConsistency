library(pcalg)
library(RBGL)
library(igraph)
library(corpcor)
source('sourceDir.R')
source('badConditionsSimulation/randomDAGparents.R')
sourceDir('utilityFunctions')
sourceDir('codeFCI')

nNodes <- 20
prob <- 0.1
nDags <- 50
numStrong <- 1
maxNeighbors <- 10
dags <- list()

#case when we work with one node, need to write moe code for multiple nodes
for (i in 1:nDags) {
  #pick a variable to connect with many others
  toChange <- sample(1:(nNodes-10),numStrong)

  dag <- randomDAGparents(nNodes, prob, toChange, maxNeighbors)
  
  dags[[i]] <- dag
}
save(dags,file = "dags.RData")

file.remove("dags.txt")
file.remove("resultsFCImajrule/pags.txt")
file.remove("resultsFCImajrule/mags.txt")
file.remove("resultsFCImajrule/sepSets.txt")

nSamples<-1000
datasets <- list()
rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)

for (i in 1:nDags) {
  dag <- dags[[i]]
  
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