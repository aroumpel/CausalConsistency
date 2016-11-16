writeBoot <- function(pag,folder,numPag) {
  write.table(pag@amat,file=paste(folder,"/bootstrap/pags",i,".txt",sep=""), append=TRUE, row.names = FALSE, col.names=FALSE)
  for (iX in 1:nNodes) {
    sepMatrix <- matrix(0,nNodes,nNodes)
    currSet <- pag@sepset[[iX]]
    if ((is.numeric(currSet))) {  # true when we look for sepSet(i,i)
      write.table(sepMatrix,file=paste(folder,"/bootstrap/sepSets",i,".txt",sep=""),append=TRUE,row.names = FALSE,col.names = FALSE)
      next
    }
    for (iY in 1:(nNodes)) {
      sepSet <- currSet[[iY]]
      sepMatrix[iY,c(sepSet)]=1
    }
    write.table(sepMatrix,file=paste(folder,"/bootstrap/sepSets",i,".txt",sep=""),append=TRUE,row.names = FALSE,col.names = FALSE)
  }
}