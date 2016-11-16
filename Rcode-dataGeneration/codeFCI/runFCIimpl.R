runFCIipml <- function(pag,folder,doSepSet,nNodes) {
  mag <- matrix(0,nNodes,nNodes)
  mag <- pag2magAM(pag@amat,1, max.chordal = nNodes)
  if(is.null(mag)) mag <- matrix(1,nNodes,nNodes)
  
  write.table(mag,file=paste(folder,"/mags.txt",sep=""), append=TRUE, row.names = FALSE, col.names=FALSE)
  write.table(pag@amat,file=paste(folder,"/pags.txt",sep=""), append=TRUE, row.names = FALSE, col.names=FALSE)
  
  if (doSepSet==TRUE) {
  #write the separating sets (for every pag we get all pairs of variables(one line) and write the separating set (0 is its empty))
    for (iX in 1:nNodes) {
      sepMatrix <- matrix(0,nNodes,nNodes)
      currSet <- pag@sepset[[iX]]
      if ((is.numeric(currSet))) {  # true when we look for sepSet(i,i)
        write.table(sepMatrix,file=paste(folder,"/sepSets.txt",sep=""),append=TRUE,row.names = FALSE,col.names = FALSE)
        next
      }
      for (iY in 1:(nNodes)) {
        sepSet <- currSet[[iY]]
        sepMatrix[iY,c(sepSet)]=1
      }
      write.table(sepMatrix,file=paste(folder,"/sepSets.txt",sep=""),append=TRUE,row.names = FALSE,col.names = FALSE)
    }
  }
}