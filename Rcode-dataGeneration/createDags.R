#Script that runs an example of creating DAGs, sampling data and running FCI.
#Package "pcalg" is used to run FCI algorithm.
#Utility "aggressivelly prevent cycles" has been added.

library(pcalg)
library(RBGL)
library(igraph)
library(corpcor)

source('sourceDir.R')
sourceDir('utilityFunctions')
sourceDir('codeFCI')

# first clear all the files
file.remove("dags.txt")

dir.create("resultsFCImajrule")
file.remove("resultsFCImajrule/pags.txt")
file.remove("resultsFCImajrule/mags.txt")
file.remove("resultsFCImajrule/sepSets.txt")

nDags <- 50;
nNodes<-20;
prob <- 0.1;
nSamples<-1000;
rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)

dags <- list()
datasets <- list()

for (i in 1:nDags) {
  print(i)
  dag <- randomDAG(nNodes,prob)
  #weights <- wgtMatrix(dag)
  dags[[i]] <- dag
  
  dag.amat <- matrix(0,nNodes,nNodes)
  
  for (j in 1:length(dag@edgeL)) {
    to <- dag@edgeL[[j]]$edges
    dag.amat[j,to] <- 1
  }
  
  write.table(dag.amat, file="dags.txt", append=TRUE, row.names=FALSE, col.names=FALSE)

  # create data
  data <- rmvDAG(nSamples,dag)
  datasets[[i]] <- data
  
  corMatrix<-cor(data)
  
  # run FCI
  suffStat <- list(C = corMatrix, n = nSamples)
  pag <- fci(suffStat,gaussCItest,alpha=0.05,p=nNodes, rules = rules, maj.rule = TRUE, aggressivelyPreventCycles = TRUE)
  runFCIipml(pag,folder="resultsFCImajrule",doSepSet=TRUE,nNodes)
}

save(dags,file = "dags.RData")
save(datasets,file = "datasets.RData")