createMarginals <- function (dataset, numRemoved, numMargs, sel) {
  nNodes <- dim(dataset)[2]
  nSamples <- dim(dataset)[1]
  marginals <- list()
  
  rules <- c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)
  
  for (i in 1:numMargs) {
    dataNew <- dataset
    selected <- sel[i,]
    dataNew <- dataNew[,selected]
    
    corMatrix<-cor(dataNew)
    suffStat <- list(C = corMatrix, n = nSamples)
    pagNew <- fci(suffStat,gaussCItest,alpha=0.05,p=nNodes-numRemoved, rules = rules, aggressivelyPreventCycles = TRUE)
    magNew <- pag2magAM(pagNew@amat,1)
    if(is.null(magNew)) magNew <- matrix(1,nNodes-numRemoved,nNodes-numRemoved)
    el <- list("pag"=pagNew,"mag"=magNew, "varsRm"=selected)
    marginals[[i]] <- el
  }
  
  return(marginals)
}