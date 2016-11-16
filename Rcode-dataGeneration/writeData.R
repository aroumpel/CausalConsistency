#Requires dataset to be loaded (for example marginalsFCI.RData)
#Also, folder must be empty...

nDags <- length(marginalsFCI)
nNodes <- 20;
nMargs <- 100;
numRemoved <- 2;
folder <- "resultsFCImajrule"
doSepSet <- "TRUE"

dir.create(paste0(folder,"/magsMarginals",numRemoved));
dir.create(paste0(folder,"/pagsMarginals",numRemoved));
dir.create(paste0(folder,"/sepSetsMarginals",numRemoved));

for (i in 1:nDags) {
  for (j in 1:length(marginalsFCI[[i]])) { #marginalsFCI[[i]][[1]]
      curMarg <- marginalsFCI[[i]][[j]] #marginalsFCI[[i]][[1]][[j]]
      
      pagNew <- curMarg$pag
      magNew <- curMarg$mag
      rm <- curMarg$varsRm
      
      write.table(magNew,file=paste0(folder,"/magsMarginals",numRemoved,"/mags",i,".txt"), append=TRUE, row.names = FALSE, col.names=FALSE)
      write.table(pagNew@amat,file=paste0(folder,"/pagsMarginals",numRemoved,"/pags",i,".txt"), append=TRUE, row.names = FALSE, col.names=FALSE)
      write.table(t(rm),file=paste0(folder,"/pagsMarginals",numRemoved,"/varsRemoved",i,".txt"), append=TRUE, row.names = FALSE, col.names=FALSE)
      
      
      if (doSepSet==TRUE) {
        #write the separating sets (for every pag we get all pairs of variables(one line) and write the separating set (0 is its empty))
        for (iX in 1:(nNodes-numRemoved)) {
          sepMatrix <- matrix(0,nNodes-numRemoved,nNodes-numRemoved)
          currSet <- pagNew@sepset[[iX]]
          if (is.numeric(currSet)) {  # true when we look for sepSet(i,i)
            write.table(sepMatrix,file=paste0(folder,"/sepSetsMarginals",numRemoved,"/sepSets",i,".txt"),append=TRUE,row.names = FALSE,col.names = FALSE)
            next
          }
          for (iY in 1:(nNodes-numRemoved)) {
            sepSet <- currSet[[iY]]
            sepMatrix[iY,c(sepSet)]=1
          }
          write.table(sepMatrix,file=paste0(folder,"/sepSetsMarginals",numRemoved,"/sepSets",i,".txt"),append=TRUE,row.names = FALSE,col.names = FALSE)
        }
      }
  }
}