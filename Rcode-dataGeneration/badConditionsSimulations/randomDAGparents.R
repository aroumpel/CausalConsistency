randomDAGparents <- function (n, prob, toChange, maxNeighbors, lB = 0.1, uB = 1, V = as.character(1:n))
{
  stopifnot(n >= 2, is.numeric(prob), length(prob) == 1,
            0 <= prob, prob <= 1,
            is.numeric(lB), is.numeric(uB), lB <= uB)
  edL <- vector("list", n)
  nmbEdges <- 0L
  for (i in seq_len(n - 2)) {
    if (i==toChange) {
      listSize <- maxNeighbors
      edgeList <- sample(seq(i+1, n), size = listSize)
    } else {
      listSize <- rbinom(1, n - i, prob)
      edgeList <- sample(seq(i+1, n), size = listSize)
    }
    nmbEdges <- nmbEdges + listSize
    weightList <- runif(length(edgeList), min = lB, max = uB)
    edL[[i]] <- list(edges = edgeList, weights = weightList)
  }
  ## i=n-1 separately
  ## (because of sample(7,1) is actually sample(1:7,1) and not 7)
  listSize <- rbinom(1, 1, prob)
  if (listSize > 0) {
    nmbEdges <- nmbEdges + 1
    edgeList <- n
    weightList <- runif(1, min = lB, max = uB)
  } else {
    edgeList <- integer(0)
    weightList <- numeric(0)
  }
  edL[[n-1]] <- list(edges = edgeList, weights = weightList)
  if (nmbEdges > 0) {
    edL[[n]] <- list(edges = integer(0), weights = numeric(0))
    names(edL) <- V
    new("graphNEL", nodes = V, edgeL = edL, edgemode = "directed")
  }
  else
    new("graphNEL", nodes = V, edgemode = "directed")
}