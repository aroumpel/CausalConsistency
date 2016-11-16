transitiveClosure <- function(amat) {
  require(relations)
  
  w <- which(amat != 0, arr.ind=T)
  R <- relation(graph = data.frame(A = w[,1], B = w[,2]))
  
  rowsInv <- unlist(R$.Data$domain$A)
  rowsToAppend <- setdiff(1:dim(amat)[1],rowsInv)
  colsInv <- unlist(R$.Data$domain$B)
  colsToAppend <- setdiff(1:dim(amat)[1],colsInv)
  
  t <- as.matrix(relation_incidence(R))
  
  for (i in rowsToAppend) {
    if (i<max(rowsInv)) {
      t <- rbind(t[1:i-1,],rep(0,dim(t)[2]),t[i:dim(t)[1],])
    } else {
      t <- rbind(t[1:nrow(t),],rep(0,dim(t)[2]))
    }
  }
  
  for(j in colsToAppend) {
    if (j<max(colsInv)) {
      t <- cbind(t[,1:j-1],rep(0,dim(t)[1]),t[,j:dim(t)[2]])
    } else {
      t <- cbind(t[,1:ncol(t)],rep(0,dim(t)[1]))
    }
  }
  
  r <- relation(incidence = t)
  
  closure <- transitive_closure(r)
  relation_incidence(closure)
  
  # the closure needs to be sorted
  rel <- as.matrix(relation_incidence(closure))
  sorted <- sort(as.numeric(rownames(rel)),index.return=TRUE)
  rel <- rel[sorted$ix,]
  rel <- rel[,sorted$ix]
  
  return(rel)
}