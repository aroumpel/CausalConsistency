findAncestors <- function(graph) {
  dag <- matrix(0,dim(graph)[1],dim(graph)[2])
  dag[graph==2 & t(graph)==3] <- 1
  isAnc <- transitiveClosure(dag)
  
  return(isAnc)
}