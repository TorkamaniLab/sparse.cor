rrBind <- function(x){
  ## recursively binds rows of sparse matrices
  ## this uses less memory than do.call(cBind, ...)
  
  require("Matrix")
  
  if(length(x) == 1) return(x[[1]])
  
  n = length(x)
  k = round(n/2)
  x1 = rrBind(x[1:k])
  x2 = rrBind(x[(k+1):n])
  
  rBind(x1, x2)
}

rcBind <- function(x){
  ## recursively binds columns of sparse matrices
  ## this uses less memory than do.call(cBind, ...)
  
  require("Matrix")
  
  if(length(x) == 1) return(x[[1]])
  
  n = length(x)
  k = round(n/2)
  x1 = rcBind(x[1:k])
  x2 = rcBind(x[(k+1):n])
  
  cBind(x1, x2)
}

get.vkd <- function(x, lambda){
  n = nrow(x)  
  w = rep(1/n, n)
  alpha = -1
  
  xs = wt.scale(x, w, center=TRUE, scale=TRUE) # standardize data matrix
  
  # bias correction factor
  h1 = 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)
  p = ncol(xs)
  
  svdxs = fast.svd(xs)
  m = length(svdxs$d)  # rank of xs
  
  UTWU = t(svdxs$u) %*% sweep(svdxs$u, 1, w, "*") #  t(U) %*% diag(w) %*% U
  C = sweep(sweep(UTWU, 1, svdxs$d, "*"), 2, svdxs$d, "*") # D %*% UTWU %*% D
  C = (1-lambda) * h1 * C
  
  C = (C + t(C))/2  # symmetrize for numerical reasons (mpower() checks symmetry)
  
  # note: C is of size m x m, and diagonal if w=1/n
  
  f = diag(m) - mpower(C/lambda + diag(m), alpha)
  v = svdxs$v
  k = f %*% t(v)
  d = (1 - sapply(1:ncol(k), function(i) sum(v[i,] * k[,i]))) * (lambda^alpha)
  
  return(list(v=v, k=k, d=d))
}


sparse.cor <- function(x, threshold=0.9, n.pieces=10, shrink=F){
  if(shrink) sparse.cor.shrink(x, threshold, n.pieces) else
         sparse.cor.noshrink(x, threshold, n.pieces)
}

sparse.cor.noshrink <- function(x, threshold, n.pieces){
  
  ## returns a sparse correlation matrix 
  ## all correlations with absolute value less than
  ## threshold are set to zero
  
  require("plyr")
  require("Matrix")
  
  cor.slice <- function(i,j){
    ind.i = (ind[i]+1):ind[i+1]
    ind.j = (ind[j]+1):ind[j+1]
    
    s = cor(x[,ind.i], x[,ind.j])
    s[abs(s) < threshold] = 0
    Matrix(s)
  }
  
  ind = round(seq(0, ncol(x), length.out=n.pieces+1))
  m = llply(1:n.pieces, 
            function(i) llply(1:n.pieces, 
                              function(j) cor.slice(i,j)),
            .progress="text")
  
  rrBind(llply(m, rcBind))
}

sparse.cor.shrink <- function(x, threshold, n.pieces){
  
  ## returns a sparse correlation matrix with shrinkage estimation
  ## see Schaefer and Strimmer (2005) for theory
  ## all correlations with absolute value less than
  ## threshold are set to zero
  
  require("plyr")
  require("Matrix")
  require("corpcor")
  
  cor.slice <- function(i,j){
    ind.i = (ind[i]+1):ind[i+1]
    ind.j = (ind[j]+1):ind[j+1]
    
    s = cor(x[,ind.i], x[,ind.j])
    a = matrix(0, nrow=nrow(s), ncol=ncol(s))
    if(i == j) diag(a) = 1
    
    s.hat = lambda*a + (1 - lambda)*s
    s.hat[abs(s.hat) < threshold] = 0
    Matrix(s.hat)
  }

  ind = round(seq(0, ncol(x), length.out=n.pieces+1))
  lambda = estimate.lambda(x)
  m = llply(1:n.pieces, 
            function(i) llply(1:n.pieces, 
                              function(j) cor.slice(i,j)),
            .progress="text")
  
  rrBind(llply(m, rcBind))
}

sparse.pcor <- function(x, threshold, n.pieces=10){
  require("corpcor")
  require("plyr")
  require("Matrix")
  
  sparse.pcor.slice <- function(i,j){
    ind.i = (ind[i]+1):ind[i+1]
    ind.j = (ind[j]+1):ind[j+1]
    
    z = (v[ind.i, ] %*% k[,ind.j])
    if(i == j) z = z - diag(nrow(z))
    z = z*(lambda^alpha)
    
    d.i = sqrt(abs(d[ind.i]))
    d.j = sqrt(abs(d[ind.j]))
    z = z / (d.i %o% d.j)
    
    if(i == j) diag(z) = 1
    
    z[abs(z) < threshold] = 0
    Matrix(z)
  }
  
  lambda = estimate.lambda(x)
  alpha = -1

  vkd = get.vkd(x, lambda)
  v = vkd$v
  k = vkd$k
  d = vkd$d
  
  ind = round(seq(0, ncol(x), length.out=n.pieces+1))
  m = llply(1:n.pieces, 
            function(i) llply(1:n.pieces, 
                              function(j) sparse.pcor.slice(i,j)),
            .progress="text")
  
  
  rrBind(llply(m, rcBind))
}

sample.cor <- function(x, n=1, shrink=F){
  if(shrink) sample.cor.shrink(x, n) else
         sample.cor.noshrink(x, n)
}

sample.cor.noshrink <- function(x, n=1){
  i = sample(1:ncol(x), n, replace=T)
  j = sample(1:ncol(x), n, replace=T)
  z = scale(x)
  k = nrow(z)
  colSums(z[ , i, drop=F] * z[, j, drop=F]) / (k - 1)
}

sample.cor.shrink <- function(x, n=1){
  require("corpcor")
  
  lambda = estimate.lambda(x)
  
  i = sample(1:ncol(x), n, replace=T)
  j = sample(1:ncol(x), n, replace=T)
  z = scale(x)
  k = nrow(z)
  
  s = colSums(z[,i,drop=F] * z[,j,drop=F]) / (k - 1)
  
  s.hat = (1 - lambda) * s
  s.hat[i == j] = 1
  s.hat
}

sample.pcor <- function(x, n=1){
  require("corpcor")
  
  lambda = estimate.lambda(x)
  alpha = -1
  
  vkd = get.vkd(x, lambda)
  v = vkd$v
  k = vkd$k
  d = vkd$d
  
  i = sample(1:ncol(x), n, replace=T)
  j = sample(1:ncol(x), n, replace=T)
  
  d.i = sqrt(abs(d[i]))
  d.j = sqrt(abs(d[j]))
  
  z = rowSums(v[i,] * t(k[,j]))
  z = z*(lambda^alpha)
  z = z / (d.i * d.j)
  z[i == j ] = 1
  
  z
}