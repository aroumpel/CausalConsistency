createsCycle <- function(aggressivelyPrevent ,graph, closure, from, to, valueFrom, valueTo) {
  cycle <- FALSE
  #res will be true if the addition of an orientation creates a cycle
  if (aggressivelyPrevent==TRUE) {
    tempGraph <- matrix(0,dim(graph)[1],dim(graph)[2])
    
    for (i in 1:length(from)) {
      graph[from[i],to[i]] <- valueTo[i]
      graph[to[i],from[i]] <- valueFrom[i]
    }
    
    closure <- updateTransitive(closure, from, to, valueFrom, valueTo)
    
    hasCycle <- any(closure & t(closure))
    hasAlmostCycle <- any(any(closure & graph==2 & t(graph)==2))
    
    if (hasCycle | hasAlmostCycle) {
      cycle <- TRUE
      print('cycle!!!')
    }
  }
  
  result <- list(cycle=cycle, closure=closure)
  
  return(result)
}