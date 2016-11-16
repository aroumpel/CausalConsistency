udag2pag <- function (pag, sepset, rules = rep(TRUE, 10), unfVect = NULL, 
    verbose = FALSE, orientCollider = TRUE, aggressivelyPreventCycles) 
{
    ###### CHANGES BY ROUMPELAKI ANNA ######
    #Added the aggressivelyPreventCycles option
    #We start with an empty transitive closure of the graph and incrementally update it
  
    stopifnot(is.logical(rules), length(rules) == 10)
    if (any(pag != 0)) {
        closure <- matrix(0,dim(pag)[1],dim(pag)[1])
      
        p <- as.numeric(dim(pag)[1])
        if (orientCollider) {
            ind <- which(pag == 1, arr.ind = TRUE)
            for (i in seq_len(nrow(ind))) {
                x <- ind[i, 1]
                y <- ind[i, 2]
                allZ <- setdiff(which(pag[y, ] != 0), x)
                for (z in allZ) {
                  if (pag[x, z] == 0 && !((y %in% sepset[[x]][[z]]) || (y %in% sepset[[z]][[x]]))) {
                    if (length(unfVect) == 0) {
                      from <- c(pag[y,x],pag[y,z])
                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, c(x,z), c(y,y), from, c(2,2))
                      cycleCreated <- result$cycle
                      
                      if (!cycleCreated) {
                        closure <- result$closure
                        if (verbose) {
                          cat("\n", x, "*->", y, "<-*", z, "\n")
                          cat("Sxz=", sepset[[z]][[x]], "and", "Szx=", sepset[[x]][[z]], "\n")
                        }
                        pag[x, y] <- pag[z, y] <- 2
                      }
                    }
                    else {
                      if (!any(unfVect == triple2numb(p, x, y, z), na.rm = TRUE) && !any(unfVect == triple2numb(p, z, y, x), na.rm = TRUE)) {
                        from <- c(pag[y,x],pag[y,z])
                        result <- createsCycle(aggressivelyPreventCycles, pag, closure, c(x,z), c(y,y), from, c(2,2))
                        cycleCreated <- result$cycle
                        
                        if (!cycleCreated) {
                          closure <- result$closure
                          if (verbose) {
                            cat("\n", x, "*->", y, "<-*", z, "\n")
                            cat("Sxz=", sepset[[z]][[x]], "and", "Szx=", sepset[[x]][[z]], "\n")
                          }
                          pag[x, y] <- pag[z, y] <- 2
                        }
                      }
                    }
                  }
                }
            }
        }
        old_pag1 <- matrix(0, p, p)
        while (any(old_pag1 != pag)) {
            old_pag1 <- pag
            if (rules[1]) {
                ind <- which((pag == 2 & t(pag) != 0), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                  a <- ind[i, 1]
                  b <- ind[i, 2]
                  indC <- which((pag[b, ] != 0 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
                  indC <- setdiff(indC, a)
                  if (length(indC) > 0) {
                    if (length(unfVect) == 0) {
                      # inject aggressivelly prevent cycles
                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, b, indC, 3, 2)
                      cycleCreated <- result$cycle
                      
                      if (!cycleCreated) {
                        closure <- result$closure
                        
                        pag[b, indC] <- 2
                        pag[indC, b] <- 3
                        if (verbose) 
                          cat("\nRule 1", "\nOrient:", a, "*->", b, "o-*", indC, "as:", b, "->", indC, "\n")
                      }
                    }
                    else {
                      for (c in indC) {
                        if (!any(unfVect == triple2numb(p, a, 
                          b, c), na.rm = TRUE) && !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                          # inject aggressivelly prevent cycles
                          result <- createsCycle(aggressivelyPreventCycles, pag, closure, b, c, 3, 2)
                          cycleCreated <- result$cycle
                          
                          if (!cycleCreated) {
                            closure <- result$closure
                            pag[b, c] <- 2
                            pag[c, b] <- 3
                            if (verbose) 
                              cat("\nRule 1", "\nConservatively orient:", a, "*->", b, "o-*", c, "as:", b, "->", c, "\n")
                          }
                        }
                      }
                    }
                  }
                }
            }
            if (rules[2]) {
                ind <- which((pag == 1 & t(pag) != 0), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                  a <- ind[i, 1]
                  c <- ind[i, 2]
                  indB <- which((pag[a, ] == 2 & pag[, a] == 3 & pag[c, ] != 0 & pag[, c] == 2) |
                                  (pag[a, ] == 2 & pag[, a] != 0 & pag[c, ] == 3 & pag[, c] == 2))
                  if (length(indB) > 0) {
                    result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, pag[c,a], 2)
                    cycleCreated <- result$cycle
                    
                    if (!cycleCreated) {
                      closure <- result$closure
                      
                      pag[a, c] <- 2
                      if (verbose) {
                        cat("\nRule 2", "\n")
                        cat("Orient:", a, "->", indB, "*->", c, "or", a, "*->", indB, "->", c, "with", a, "*-o", c, "as:", a, "*->", c, "\n")
                      }
                    }
                  }
                }
            }
            if (rules[3]) {
                ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                  b <- ind[i, 1]
                  d <- ind[i, 2]
                  indAC <- which((pag[b, ] != 0 & pag[, b] == 2) & (pag[, d] == 1 & pag[d, ] != 0))
                  if (length(indAC) >= 2) {
                    if (length(unfVect) == 0) {
                      counter <- 0
                      while ((counter < (length(indAC) - 1)) && (pag[d, b] != 2)) {
                        counter <- counter + 1
                        ii <- counter
                        while (ii < length(indAC) && pag[d, b] != 2) {
                          ii <- ii + 1
                          if (pag[indAC[counter], indAC[ii]] == 0 && pag[indAC[ii], indAC[counter]] == 0) {
                            result <- createsCycle(aggressivelyPreventCycles, pag, closure, d, b, pag[b,d], 2)
                            cycleCreated <- result$cycle
                            
                            if (!cycleCreated) {
                              closure <- result$closure
                              
                              if (verbose) {
                                cat("\nRule 3", "\n")
                                cat("Orient:", d, "*->", b, "\n")
                              }
                              pag[d, b] <- 2
                            }
                          }
                        }
                      }
                    }
                    else {
                      comb.indAC <- combn(indAC, 2)
                      for (j in 1:dim(comb.indAC)[2]) {
                        a <- comb.indAC[1, j]
                        c <- comb.indAC[2, j]
                        if (pag[a, c] == 0 && pag[c, a] == 0 && c != a) {
                          if (!any(unfVect == triple2numb(p, a, d, c), na.rm = TRUE) && !any(unfVect == triple2numb(p, c, d, a), na.rm = TRUE)) {
                            result <- createsCycle(aggressivelyPreventCycles, pag, closure, d, b, pag[b,d], 2)
                            cycleCreated <- result$cycle
                            
                            if (!cycleCreated) {
                              closure <- result$closure
                              
                              pag[d, b] <- 2
                              if (verbose) {
                                cat("\nRule 3", "\n")
                                cat("Conservatively orient:", d, "*->", b, "\n")
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
            if (rules[4]) {
                ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                  b <- ind[1, 1]
                  c <- ind[1, 2]
                  ind <- ind[-1, , drop = FALSE]
                  indA <- which((pag[b, ] == 2 & pag[, b] != 0) & (pag[c, ] == 3 & pag[, c] == 2))
                  while (length(indA) > 0 && pag[c, b] == 1) {
                    a <- indA[1]
                    indA <- indA[-1]
                    Done <- FALSE
                    while (!Done && pag[a, b] != 0 && pag[a, c] != 0 && pag[b, c] != 0) {
                      md.path <- minDiscrPath(pag, a, b, c, verbose = verbose)
                      if ((N.md <- length(md.path)) == 1) {
                        Done <- TRUE
                      }
                      else {
                        if ((b %in% sepset[[md.path[1]]][[md.path[N.md]]]) || 
                          (b %in% sepset[[md.path[N.md]]][[md.path[1]]])) {
                          result <- createsCycle(aggressivelyPreventCycles, pag, closure, b, c, 3, 2)
                          cycleCreated <- result$cycle
                          
                          if (!cycleCreated) {
                            closure <- result$closure
                            
                            if (verbose) 
                              cat("\nRule 4", "\nThere is a discriminating path between", md.path[1], "and", c, "for", b, 
                                ",and", b, "is in Sepset of", c, "and", md.path[1], ". Orient:", b, "->", c, "\n")
                              pag[b, c] <- 2
                              pag[c, b] <- 3
                          }
                        }
                        else {
                          result <- createsCycle(aggressivelyPreventCycles, pag, closure, c(b,a), c(c,b), c(2,2), c(2,2))
                          cycleCreated <- result$cycle
                          
                          if (!cycleCreated) {
                            closure <- result$closure
                            
                            if (verbose) 
                              cat("\nRule 4", "\nThere is a discriminating path between", md.path[1], "and", c, "for", b, 
                                ",and", b, "is not in Sepset of", c, "and", md.path[1], ". Orient:", a, "<->", b, "<->", c, "\n")
                              pag[a, b] <- pag[b, c] <- pag[c, b] <- 2
                          }
                        }
                        Done <- TRUE
                      }
                    }
                  }
                }
            }
            if (rules[5]) {
                ind <- which((pag == 1 & t(pag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                  a <- ind[1, 1]
                  b <- ind[1, 2]
                  ind <- ind[-1, , drop = FALSE]
                  indC <- which((pag[a, ] == 1 & pag[, a] == 1) & (pag[b, ] == 0 & pag[, b] == 0))
                  indC <- setdiff(indC, b)
                  indD <- which((pag[b, ] == 1 & pag[, b] == 1) & (pag[a, ] == 0 & pag[, a] == 0))
                  indD <- setdiff(indD, a)
                  if (length(indC) > 0 && length(indD) > 0) {
                    counterC <- 0
                    while ((counterC < length(indC)) && pag[a, b] == 1) {
                      counterC <- counterC + 1
                      c <- indC[counterC]
                      counterD <- 0
                      while ((counterD < length(indD)) && pag[a, b] == 1) {
                        counterD <- counterD + 1
                        d <- indD[counterD]
                        if (pag[c, d] == 1 && pag[d, c] == 1) {
                          if (length(unfVect) == 0) {
                            from <- c(a,a,c,d)
                            to <- c(b,c,d,b)
                            values <- rep(3,4)
                            
                            result <- createsCycle(aggressivelyPreventCycles, pag, closure, from, to, values, values)
                            cycleCreated <- result$cycle
                            
                            if (!cycleCreated) {
                              closure <- result$closure
                              
                              pag[a, b] <- pag[b, a] <- 3
                              pag[a, c] <- pag[c, a] <- 3
                              pag[c, d] <- pag[d, c] <- 3
                              pag[d, b] <- pag[b, d] <- 3
                              if (verbose) 
                                cat("\nRule 5", "\nThere exists an uncovered circle path between", a, "and", b, ". Orient:", a, 
                                  "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                            }
                          }
                          else {
                            path2check <- c(a, c, d, b)
                            if (faith.check(path2check, unfVect, p)) {
                              from <- c(a,b,a,c,c,d,d,b)
                              to <- c(b,a,c,a,d,c,b,d)
                              values <- rep(3,8)
                              
                              result <- createsCycle(aggressivelyPreventCycles, pag, closure, from, to, values, values)
                              cycleCreated <- result$cycle
                              
                              if (!cycleCreated) {
                                closure <- result$closure
                                
                                pag[a, b] <- pag[b, a] <- 3
                                pag[a, c] <- pag[c, a] <- 3
                                pag[c, d] <- pag[d, c] <- 3
                                pag[d, b] <- pag[b, d] <- 3
                                if (verbose) 
                                  cat("\nRule 5", "\nThere exists a faithful uncovered circle path between", a, "and", b, ". Conservatively orient:", 
                                    a, "-", b, "and", a, "-", c, "-", d, "-", b, "\n")
                              }
                            }
                          }
                        }
                        else {
                          ucp <- minUncovCircPath(p, pag = pag, path = c(a, c, d, b), unfVect = unfVect, verbose = verbose)
                          if (length(ucp) > 1) {
                            n <- length(ucp)
                            
                            result <- createsCycle(aggressivelyPreventCycles, pag, closure, ucp[1], ucp[n], 3, 3)
                            cycleCreated <- result$cycle
                            
                            if (!cycleCreated) {
                              closure <- result$closure
                              
                              pag[ucp[1], ucp[n]] <- pag[ucp[n], ucp[1]] <- 3
                            }
                            for (j in 1:(length(ucp) - 1)) {
                              result <- createsCycle(aggressivelyPreventCycles, pag, closure, c(ucp[j],ucp[j+1]), c(ucp[j+1],ucp[j]), c(3,3), c(3,3))
                              cycleCreated <- result$cycle
                              
                              if (!cycleCreated) {
                                closure <- result$closure
                                
                                pag[ucp[j], ucp[j + 1]] <- pag[ucp[j + 1], ucp[j]] <- 3
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
            if (rules[6]) {
                ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                  b <- ind[i, 1]
                  c <- ind[i, 2]
                  if (any(pag[b, ] == 3 & pag[, b] == 3)) {
                    result <- createsCycle(aggressivelyPreventCycles, pag, closure, c, b, pag[b,c], 3)
                    cycleCreated <- result$cycle
                    
                    if (!cycleCreated) {
                      closure <- result <- closure
                      
                      pag[c, b] <- 3
                      if (verbose) 
                        cat("\nRule 6", "\nOrient:", b, "o-*", c, "as", b, "-*", c, "\n")
                    }
                  }
                }
            }
            if (rules[7]) {
                ind <- which((pag != 0 & t(pag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                  b <- ind[i, 1]
                  c <- ind[i, 2]
                  indA <- which((pag[b, ] == 3 & pag[, b] == 1) & (pag[c, ] == 0 & pag[, c] == 0))
                  indA <- setdiff(indA, c)
                  if (length(indA) > 0) {
                    if (length(unfVect) == 0) {
                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, c, b, pag[b,c], 3)
                      cycleCreated <- result$cycle
                      
                      if (!cycleCreated) {
                        closure <- result$closure
                        
                        pag[c, b] <- 3
                        if (verbose) 
                          cat("\nRule 7", "\nOrient:", indA, "-o", b, "o-*", c, "as", b, "-*", c, "\n")
                      }
                    }
                    else for (a in indA) if (!any(unfVect == triple2numb(p, a, b, c), na.rm = TRUE) && 
                      !any(unfVect == triple2numb(p, c, b, a), na.rm = TRUE)) {
                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, c, b, pag[b,c], 3)
                      cycleCreated <- result$cycle
                      
                      if (!cycleCreated) {
                        closure <- result$closure
                        
                        pag[c, b] <- 3
                        if (verbose) 
                          cat("\nRule 7", "\nConservatively orient:", a, "-o", b, "o-*", c, "as", b, "-*", c, "\n")
                      }
                    }
                  }
                }
            }
            if (rules[8]) {
                ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
                for (i in seq_len(nrow(ind))) {
                  a <- ind[i, 1]
                  c <- ind[i, 2]
                  indB <- which(pag[, a] == 3 & (pag[a, ] == 2 | pag[a, ] == 1) & pag[c, ] == 3 & pag[, c] == 2)
                  if (length(indB) > 0) {
                    result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, 3, pag[a,c])
                    cycleCreated <- result$cycle
                    
                    if (!cycleCreated) {
                      closure <- result$closure
                      
                      pag[c, a] <- 3
                      if (verbose) 
                        cat("\nRule 8", "\nOrient:", a, "->", indB, "->", c, "or", a, "-o", indB, "->", c, 
                          "with", a, "o->", c, "as", a, "->", c, "\n")
                    }
                  }
                }
            }
            if (rules[9]) {
                ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                  a <- ind[1, 1]
                  c <- ind[1, 2]
                  ind <- ind[-1, , drop = FALSE]
                  indB <- which((pag[a, ] == 2 | pag[a, ] == 1) & (pag[, a] == 1 | pag[, a] == 3) & (pag[c, ] == 0 & pag[, c] == 0))
                  indB <- setdiff(indB, c)
                  while ((length(indB) > 0) && (pag[c, a] == 1)) {
                    b <- indB[1]
                    indB <- indB[-1]
                    upd <- minUncovPdPath(p, pag, a, b, c, unfVect = unfVect, verbose = verbose)
                    if (length(upd) > 1) {
                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, 3, pag[a,c])
                      cycleCreated <- result$cycle
                      
                      if (!cycleCreated) {
                        closure <- result$closure
                        
                        pag[c, a] <- 3
                        if (verbose) 
                          cat("\nRule 9", "\nThere exists an uncovered potentially directed path between", 
                            a, "and", c, ". Orient:", a, " ->", c, "\n")
                      }
                    }
                  }
                }
            }
            if (rules[10]) {
                ind <- which((pag == 2 & t(pag) == 1), arr.ind = TRUE)
                while (length(ind) > 0) {
                  a <- ind[1, 1]
                  c <- ind[1, 2]
                  ind <- ind[-1, , drop = FALSE]
                  indB <- which((pag[c, ] == 3 & pag[, c] == 2))
                  if (length(indB) >= 2) {
                    counterB <- 0
                    while (counterB < length(indB) && (pag[c, a] == 1)) {
                      counterB <- counterB + 1
                      b <- indB[counterB]
                      indD <- setdiff(indB, b)
                      counterD <- 0
                      while ((counterD < length(indD)) && (pag[c, a] == 1)) {
                        counterD <- counterD + 1
                        d <- indD[counterD]
                        if ((pag[a, b] == 1 || pag[a, b] == 2) && 
                          (pag[b, a] == 1 || pag[b, a] == 3) && 
                          (pag[a, d] == 1 || pag[a, d] == 2) && 
                          (pag[d, a] == 1 || pag[d, a] == 3) && 
                          pag[d, b] == 0 && pag[b, d] == 0) {
                          if (length(unfVect) == 0) {
                            result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, 3, pag[a,c])
                            cycleCreated <- result$cycle
                            
                            if (!cycleCreated) {
                              closure <- result$closure
                              
                              pag[c, a] <- 3
                              if (verbose) 
                                cat("\nRule 10 [easy]", "\nOrient:", a, "->", c, "\n")
                            }
                          }
                          else if (!any(unfVect == triple2numb(p, b, a, d), na.rm = TRUE) && 
                                   !any(unfVect == triple2numb(p, d, a, b), na.rm = TRUE)) {
                            result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, 3, pag[a,c])
                            cycleCreated <- result$cycle
                            
                            if (!cycleCreated) {
                              closure <- result$closure
                              
                              pag[c, a] <- 3
                              if (verbose) 
                                cat("\nRule 10 [easy]", "\nConservatively orient:", a, "->", c, "\n")
                            }
                          }
                        }
                        else {
                          indX <- which((pag[a, ] == 1 | pag[a, ] == 2) & (pag[, a] == 1 | pag[, a] == 3), arr.ind = TRUE)
                          indX <- setdiff(indX, c)
                          if (length(indX >= 2)) {
                            counterX1 <- 0
                            while (counterX1 < length(indX) && pag[c, a] == 1) {
                              counterX1 <- counterX1 + 1
                              first.pos <- indA[counterX1]
                              indX2 <- setdiff(indX, first.pos)
                              counterX2 <- 0
                              while (counterX2 < length(indX2) && pag[c, a] == 1) {
                                counterX2 <- counterX2 + 1
                                sec.pos <- indX2[counterX2]
                                t1 <- minUncovPdPath(p, pag, a, first.pos, b, unfVect = unfVect, verbose = verbose)
                                if (length(t1) > 1) {
                                  t2 <- minUncovPdPath(p, pag, a, sec.pos, d, unfVect = unfVect, verbose = verbose)
                                  if (length(t2) > 1 && first.pos != sec.pos && pag[first.pos, sec.pos] == 0) {
                                    if (length(unfVect) == 0) {
                                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, 3, pag[a,c])
                                      cycleCreated <- result$cycle
                                      
                                      if (!cycleCreated) {
                                        closure <- result$closure
                                        
                                        pag[c, a] <- 3
                                        if (verbose) 
                                          cat("\nRule 10", "\nOrient:", a, "->", c, "\n")
                                      }
                                    }
                                    else if (!any(unfVect == triple2numb(p, first.pos, a, sec.pos), na.rm = TRUE) && 
                                             !any(unfVect == triple2numb(p, sec.pos, a, first.pos), na.rm = TRUE)) {
                                      result <- createsCycle(aggressivelyPreventCycles, pag, closure, a, c, 3, pag[a,c])
                                      cycleCreated <- result$cycle
                                      
                                      if (!cycleCreated) {
                                        closure <- result$closure
                                        
                                        pag[c, a] <- 3
                                        if (verbose) 
                                          cat("\nRule 10", "\nConservatively orient:", a, "->", c, "\n")
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
            }
        }
    }
  pag
}