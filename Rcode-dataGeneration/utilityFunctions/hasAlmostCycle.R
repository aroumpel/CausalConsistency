hasAlmostCycle <- function(graph) {
  isDesc <- findDescendants(graph)
  res <- any(any((isDesc | t(isDesc)) & graph==2 & t(graph)==2))
  return(res)
}