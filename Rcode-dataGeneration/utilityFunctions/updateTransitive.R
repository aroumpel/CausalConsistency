updateTransitive <- function(closure, from, to, valueFrom, valueTo) {
  #to update the transitive closure when we add i->j, we need five steps
  #1. identify the ancestors of i, A
  #2. identify the descendants of j, D
  #3. for every a in A, set closure(a,j)=1
  #4. for every d in D, set closure(i,d)=1
  #5. for every a in A, d in D, set closure(a,d)=1
  # of course set closure(i,j)=1
  
  a1 <- which(valueFrom==3);   a2 <- which(valueFrom==2);
  b1 <- which(valueTo==2);     b2 <- which(valueTo==3);
  
  ind1 <- intersect(a1,b1);   ind2 <- intersect(a2,b2);
  
  closure[from[ind1],to[ind1]] <- 1
  closure[to[ind2],from[ind2]] <- 1
  
  # find the ancestors of "from" and the descendants of "to"
  for (i in 1:length(from)) {
    A <- which(closure[,from[i]]==1)
    D <- which(closure[to[i],]==1)
    
    closure[A,to[i]] <- 1
    closure[from[i],D] <- 1
    closure[A,D] <- 1
  }
  
  return(closure)
}