findDescendants <- function(graph) {
  dag <- matrix(0,dim(graph)[1],dim(graph)[2])
  dag[graph==3 & t(graph)==2] <- 1
  isDesc <- transitiveClosure(dag)
  
  return(isDesc)
}