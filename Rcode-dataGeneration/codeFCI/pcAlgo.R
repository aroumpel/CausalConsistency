pcAlgo <- function(dm = NA, C = NA, n = NA, alpha, corMethod = "standard",
                   verbose = FALSE, directed = FALSE,
                   G = NULL, datatype = 'continuous', NAdelete = TRUE,
                   m.max = Inf, u2pd = "rand", psepset = FALSE) {
  ## Purpose: Perform PC-Algorithm, i.e., estimate skeleton of DAG given data
  ## Output is an unoriented graph object
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - dm: Data matrix (rows: samples, cols: nodes)
  ## - C: correlation matrix (only for continuous)
  ## - n: sample size
  ## - alpha: Significance level of individual partial correlation tests
  ## - corMethod: "standard" or "Qn" for standard or robust correlation
  ##              estimation
  ## - G: the adjacency matrix of the graph from which the algorithm
  ##      should start (logical)
  ## - datatype: distinguish between discrete and continuous data
  ## - NAdelete: delete edge if pval=NA (for discrete data)
  ## - m.max: maximal size of conditioning set
  ## - u2pd: Function for converting udag to pdag
  ##   "rand": udag2pdagu
  ##   "relaxed": udag2pdagRelaxed
  ##   "retry": udag2pdagSpecial
  ## - psepset: Also check possible sep sets.
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date: 26 Jan 2006; Martin Maechler
  ## Modifications: Sarah Gerster, Diego Colombo
  
  .Deprecated(msg = "pcAlgo() is deprecated and only kept for backward compatibility.
              Please use skeleton, pc, or fci instead\n")
  cl <- match.call()
  
  if (any(is.na(dm))) {
    stopifnot(all(!is.na(C)),!is.na(n), (p <- ncol(C)) > 0)
  } else {
    n <- nrow(dm)
    p <- ncol(dm)
  }
  n <- as.integer(n)
  
  if (is.null(G)) {
    ## G := complete graph :
    G <- matrix(TRUE, p,p)
    diag(G) <- FALSE
  } else if (!(identical(dim(G),c(p,p))))
    stop("Dimensions of the dataset and G do not agree.")
  
  seq_p <- seq_len(p)
  sepset <- pl <- vector("list",p)
  for (i in seq_p) sepset[[i]] <- pl
  zMin <- matrix(Inf, p,p)
  n.edgetests <- numeric(1)# final length = max { ord}
  done <- FALSE
  ord <- 0
  
  if (datatype == 'continuous') {
    diag(zMin) <- 0
    if (any(is.na(C))) C <- mcor(dm, method = corMethod)
    cutoff <- qnorm(1 - alpha/2)
    while (!done && any(G) && ord <= m.max) {
      n.edgetests[ord+1] <- 0
      done <- TRUE
      ind <- which(G, arr.ind = TRUE)
      ## For comparison with C++ sort according to first row
      ind <- ind[order(ind[,1]), ]
      remEdges <- nrow(ind)
      if(verbose)
        cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
      for (i in 1:remEdges) {
        if(verbose && i%%100 == 0) cat("|i=",i,"|iMax=",remEdges,"\n")
        x <- ind[i,1]
        y <- ind[i,2]
        if (G[y,x]) {
          nbrsBool <- G[,x]
          nbrsBool[y] <- FALSE
          nbrs <- seq_p[nbrsBool]
          length_nbrs <- length(nbrs)
          if (length_nbrs >= ord) {
            if (length_nbrs > ord) done <- FALSE
            S <- seq(length = ord)
            repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
              n.edgetests[ord+1] <- n.edgetests[ord+1]+1
              z <- zStat(x,y, nbrs[S], C,n)
              if (verbose) cat(paste("x:",x,"y:",y,"S:"),nbrs[S],paste("z:",z,"\n"))
              if(abs(z) < zMin[x,y]) zMin[x,y] <- abs(z)
              if (abs(z) <= cutoff) {
                G[x,y] <- G[y,x] <- FALSE
                sepset[[x]][[y]] <- nbrs[S]
                break
              } else {
                nextSet <- getNextSet(length_nbrs, ord, S)
                if(nextSet$wasLast)
                  break
                S <- nextSet$nextSet
              }
            }
          }
        } ## end if(!done)
        
      } ## end for(i ..)
      ord <- ord+1
      ##    n.edgetests[ord] <- remEdges
    } ## while
    
    for (i in 1:(p-1)) {
      for (j in 2:p) {
        zMin[i,j] <- zMin[j,i] <- min(zMin[i,j],zMin[j,i])
      }
    }
  }
  else {
    ##
    ##
    ## DISCRETE DATA ######################################################
    ##
    if (datatype == 'discrete') {
      dm.df <- as.data.frame(dm)
      while (!done && any(G) && ord <= m.max) {
        n.edgetests[ord+1] <- 0
        done <- TRUE
        ind <- which(G, arr.ind = TRUE)
        ## For comparison with C++ sort according to first row
        ind <- ind[order(ind[,1]), ]
        remEdges <- nrow(ind)
        if(verbose)
          cat("Order=",ord,"; remaining edges:",remEdges,"\n", sep = '')
        for (i in 1:remEdges) {
          if(verbose) { if(i%%100 == 0) cat("|i=",i,"|iMax=",remEdges,"\n") }
          x <- ind[i,1]
          y <- ind[i,2]
          if (G[y,x]) {
            nbrsBool <- G[,x]
            nbrsBool[y] <- FALSE
            nbrs <- seq_p[nbrsBool]
            length_nbrs <- length(nbrs)
            if (length_nbrs >= ord) {
              if (length_nbrs > ord) done <- FALSE
              S <- seq(length = ord)
              repeat { ## condition w.r.to all  nbrs[S] of size 'ord'
                n.edgetests[ord+1] <- n.edgetests[ord+1]+1
                prob <- ci.test(x,y, nbrs[S], dm.df)
                if (verbose) cat("x=",x," y=",y," S=",nbrs[S],":",prob,"\n")
                if (is.na(prob)) prob <- if(NAdelete) 1 else 0
                if(prob >= alpha) { # independent
                  G[x,y] <- G[y,x] <- FALSE
                  sepset[[x]][[y]] <- nbrs[S]
                  break
                } else {
                  nextSet <- getNextSet(length_nbrs, ord, S)
                  if(nextSet$wasLast)
                    break
                  S <- nextSet$nextSet
                }
              }
            }
          } ## end if(!done)
          
        } ## end for(i ..)
        ord <- ord+1
        ##    n.edgetests[ord] <- remEdges
      } ## while
    } else
      stop("Datatype must be 'continuous' or 'discrete'.")
  }
  
  if (psepset) {
    amat <- G
    ind <- which(G, arr.ind = TRUE)
    storage.mode(amat) <- "integer" # (TRUE, FALSE) -->  (1, 0)
    ## Orient colliders
    for (i in seq_len(nrow(ind))) {
      x <- ind[i,1]
      y <- ind[i,2]
      allZ <- setdiff(which(amat[y,] == 1),x) ## x-y-z
      
      for (z in allZ) {
        if (amat[x,z] == 0 &&
            !((y %in% sepset[[x]][[z]]) ||
              (y %in% sepset[[z]][[x]]))) {
          if (verbose >= 2) {
            cat("\n",x,"*->",y,"<-*",z,"\n")
            cat("Sxz=",sepset[[z]][[x]],"and","Szx=",sepset[[x]][[z]],"\n")
          }
          
          ## x o-> y <-o z
          amat[x,y] <- amat[z,y] <- 2
          
        } ## for
      } ## if
    } ## for
    
    ## Compute poss. sepsets
    for (x in 1:p) {
      attr(x,'class') <- 'possibledsep'
      if (any(amat[x,] != 0)) {
        tf1 <- setdiff(reach(x,-1,-1,amat), x)
        for (y in seq_p[amat[x,] != 0]) {
          ## tf = possible_d_sep(amat,x,y)
          tf <- setdiff(tf1,y)
          ## test
          if (length(tf) > 0) {
            az <- abs(zStat(x,y,tf,C,n))
            if (az < zMin[x,y]) zMin[x,y] <- az
            if (az <= cutoff) {
              ## delete x-y
              amat[x, y] <- amat[y, x] <- 0
              ## save pos d-sepset in sepset
              sepset[[x]][[y]] <- tf
            }
            if (verbose >= 2)
              cat("Possible-D-Sep of", x, "and", y, "is", tf, " - |z| = ",az,"\n")
          }
        }
      }
    }
    G[amat == 0] <- FALSE
    G[amat == 1] <- TRUE
  } ## end if(psepset)
  
  if(verbose) { cat("Final graph adjacency matrix:\n"); print(symnum(G)) }
  
  ## transform matrix to graph object (if not deprecated anyway: FIX to use correct node names!)
  Gobject <- if (sum(G) == 0) {
    new("graphNEL", nodes = as.character(seq_p))
  } else {
    colnames(G) <- rownames(G) <- as.character(seq_p)
    as(G,"graphNEL")
  }
  
  res <- new("pcAlgo", graph = Gobject,
             call = cl, n = n, max.ord = as.integer(ord-1),
             n.edgetests = n.edgetests, sepset = sepset,
             zMin = zMin)
  if (directed)
    switch (u2pd,
            "rand"    = udag2pdag       (res),
            "retry"   = udag2pdagSpecial(res)$pcObj,
            "relaxed" = udag2pdagRelaxed(res))
  else
    res
} ## {pcAlgo} __ deprecated __