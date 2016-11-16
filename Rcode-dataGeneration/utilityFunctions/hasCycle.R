hasCycle <- function(graph) {
  isDesc <- findDescendants(graph)
  res <- any(isDesc & t(isDesc))
  return(res)
}