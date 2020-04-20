options(warn = -1, digits = 4)

#' Compute the riemann sum
#' @param domain is the a vector consisting of time points where the integrand is evaluated.
#' @param integrand is a matrix. Each column correspond to one integrand. 
#' length(integrand) == length(domain) if integrand is a vector
#' nrow(integrand) == length(domain) if integrand is a matrix

riemannSum = function (domain, integrand) {  
  if (is.vector(integrand))
    result = sum(integrand) * diff(range(domain)) / (length(domain) - 1)
  if (is.matrix(integrand))
    result = colSums(integrand) * diff(range(domain)) / (length(domain) - 1)
  return(result)
}

#' Expand 1-dim functional data onto a B-spline basis.
#' @param x is a matrix. Each row corresponds to a curve.
#' @param nOrder = 4 by default indicates the cubic B-spline.

expand.bspline = function(x, nOrder = 4){
  if (!("fda" %in% rownames(installed.packages()))) install.packages("fda")
  if (!all(c(
    exists("create.bspline.basis", mode="function"),
    exists("int2Lfd", mode="function"),
    exists("fdPar", mode="function"),
    exists("Data2fd", mode="function"),
    exists("inprod", mode="function")
  ))) library("fda")
  
  timePts = 1:ncol(x)
  K = min(nOrder+length(timePts)-2, nrow(x)-1) # number of basis functions
  splineBasis = create.bspline.basis(rangeval=c(1, length(timePts)), nbasis=K, norder=nOrder, names="bspl")
  D2Lfd = int2Lfd(m=2)
  
  # tune the value of lambda for smoothing x
  log10lambda = 0:6
  gcvsave = NULL
  for (i in log10lambda) {
    lambda   = 10^i
    D2fdPar = fdPar(splineBasis, D2Lfd, lambda)
    smooth = smooth.basis(argvals=timePts, y=t(x), D2fdPar, dfscale=1.2)
    gcvsave = c(gcvsave, sum(smooth$gcv))
  }
  
  lambda = 10^log10lambda[which.min(gcvsave)]
  D2fdPar = fdPar(splineBasis, D2Lfd, lambda)
  fd.train = Data2fd(y=t(x), argvals=timePts, D2fdPar, nderiv=2, lambda)
  
  return(list(coef=t(fd.train$coef), basis=splineBasis, D2fdPar=D2fdPar, lambda=lambda))
} 

#' search for the minimizer of f
#' @param ninitial is the number of initial values.

manymodes = function(f, ninitial){
  
  initial1 = -1/(ninitial+1)*(1:ninitial)
  initial2 = 1e8/(ninitial+1)*(1:ninitial)
  
  for (i in 1:length(initial1)) {
    result.tmp = optim(initial1[i], f, method=c("L-BFGS-B"), lower=-1+1e-6, upper=-1e-6)
    if (i==1) result = result.tmp
    if (result.tmp$value<result$value) result = result.tmp
  }
  for (i in 1:length(initial2)) {
    result.tmp = optim(initial2[i], f, method=c("L-BFGS-B"), lower=1e-5, upper=Inf)
    if (result.tmp$value<result$value) result = result.tmp
  }
  
  x.star = result$par
  return(x.star)
}

#' compute A^.5 for symmetric matrix A

half = function(A){
  svdA = svd(A)
  E.half = diag(sqrt(svdA$d))
  U = svdA$u
  V = svdA$v
  return(U %*% E.half %*% t(V))
}

#' compute A^-.5 for symmetric matrix A

halfinv = function(A){
  svdA = svd(A)
  tmp = svdA$d[abs(svdA$d) > 1e-6]
  E.halfinv = diag(c(1/sqrt(tmp), rep(0, times=length(svdA$d)-length(tmp))))
  U = svdA$u
  V = svdA$v
  return(U %*% E.halfinv %*% t(V))
}

#' implement SIMPLS algorithm in de Jong (1993, Chemometr. Intell. Lab. vol. 18, pp. 251-263).
#' interim function for @function pFPLS
#' @input n*p matrix x, n*m matrix y and number of components picked

simpls = function(y, x, ncomp){

  Rmat = NULL
  Tmat = NULL
  Pmat = NULL
  Qmat = NULL
  Umat = NULL
  Vmat = NULL
  
  if (ncol(y) > 1) Y0 = sweep(y, 2, colMeans(y))
  else Y0 = y - mean(y)
  X = x
  
  S = crossprod(X, Y0)
  
  for (a in 1:ncomp){
    q = svd(crossprod(S), nu=1, nv=1)$u # Y block factor weights
    r = S %*% q # X block factor weights
    t = X %*% r # X block factor scores
    t = t - mean(t) # center scores
    normt = sum(t^2)^.5
    t = t/normt # normalize scores
    r = r/normt # adapt weights a~~ordingIy
    p = crossprod(X, t) # X block factor loadings
    q = crossprod(Y0, t) # Y block factor loadings
    u = Y0 %*% q # Y block factor scores
    v = p # initialize orthogonal loadings
    if (a>1){
      v = v - Vmat %*% crossprod(Vmat, p) # make v perpendicular of previous loadings
      u = u - Tmat %*% crossprod(Tmat, u) # make u perpendicular of previous t values
    }
    v = v / sum(v^2)^.5 # normalize orthogonal loadings
    S = S - v %*% crossprod(v, S) # deflate S with respect to current loadings
    
    Rmat = cbind(Rmat, r)
    Tmat = cbind(Tmat, t)
    Pmat = cbind(Pmat, p)
    Qmat = cbind(Qmat, q)
    Umat = cbind(Umat, u)
    Vmat = cbind(Vmat, v)
  }
  
  B = tcrossprod(Rmat, Qmat) # PLS regression coefficients
  
  return(list(
    R = Rmat,
    coef = B
  ))
}

#' interim function for @function fpls

fpcr.setup <- function(y, xfuncs = NULL, nbasis = NULL, basismat = NULL, penmat = NULL, argvals = NULL,
                       covt = NULL, mean.signal.term = FALSE, spline.order = NULL, fdobj = NULL, pen.order = 2) {
  if (!is.null(fdobj)) {
    # cat("performing scalar on functional regression\n")
    basis <- fdobj$basis
    if (!is.null(nbasis)){
      if (length(nbasis) > 1){
        warning("nbasis =", nbasis[-1], "will be ignored. Only the first nbasis, nbasis =",
                nbasis[1],  "will be considerred")
      }
      if (nbasis[1] != basis$nbasis){
        warning(paste("'nbasis =", nbasis, "'overridden since the supplied 'fdobj'
                      has basis dimension", basis$nbasis))
      }
      }
    if (!is.null(spline.order)) {
      if(spline.order != norder(basis)){
        warning(paste("'spline.order =", spline.order, "'overridden since the supplied
                      'fdobj' has a basis of order", norder(basis)))
      }
      }
    if (is.null(argvals)){
      argvals <- seq(basis$rangeval[1], basis$rangeval[2], length.out = 401)
    }
    xfuncs <- t(eval.fd(argvals, fdobj))
    }
  
  if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:3)
    stop("xfuncs must either be a 3D or 2D array")
  
  dim.sig <- length(dim(xfuncs)) - 1
  if (dim.sig == 1) {
    if (is.null(fdobj)){
      if (!is.null(nbasis) && length(nbasis) > 1) {
        warning("nbasis = ", nbasis[-1], " will be ignored. Only the first nbasis, nbasis = ",
                nbasis[1], " will be considerred")
      }
      if (is.null(argvals)) argvals <- seq(0, 1, length.out=ncol(xfuncs))
      if (is.null(nbasis)) nbasis <- 40
      if (is.null(spline.order)) spline.order <- 4
      basis <- create.bspline.basis(rangeval = c(min(argvals), max(argvals)),
                                    nbasis = nbasis[1], norder = spline.order)
    }
    if (is.null(basismat)) basismat <- eval.basis(argvals, basis)
    if (is.null(penmat)) penmat <- getbasispenalty(basis, pen.order)
  }
  else {
    d1 <- dim(xfuncs)[2]
    d2 <- dim(xfuncs)[3]
    if (is.null(argvals)) argvals <- list(seq(0, 1, length.out=d1), seq(0, 1, length.out=d2))
    if (is.null(nbasis)) {
      nbasis <- c(15, 15)
    } else if (length(nbasis) == 1){
      warning("nbasis = ", nbasis[1], " will be applied to both direction")
      nbasis <- c(nbasis, nbasis)
    } else if (length(nbasis) > 2){
      warning("only nbasis = ", nbasis[1:2], " will be considerred")
      nbasis <- nbasis[1:2]
    }
    if (is.null(spline.order)) spline.order <- 4
    basis1 <- create.bspline.basis(rangeval = c(min(argvals[[1]]), max(argvals[[1]])),
                                   nbasis = nbasis[1], norder = spline.order)
    basis2 <- create.bspline.basis(rangeval = c(min(argvals[[2]]), max(argvals[[2]])),
                                   nbasis = nbasis[2], norder = spline.order)
    if (is.null(basismat)){
      basismat <- eval.basis(argvals[[2]], basis2) %x% eval.basis(argvals[[1]], basis1)
    }
    if (is.null(penmat)) penmat <- osplinepen2d(basis1, basis2)
    dim(xfuncs) <- c(dim(xfuncs)[1], d1 * d2)
  }
  X0 <- if (mean.signal.term) cbind(1, apply(xfuncs, 1, mean))
  else matrix(1, length(y), 1)
  if (!is.null(covt)) X0 <- cbind(X0, as.matrix(covt))
  sigs.decor <- if (ncol(X0) == 1) scale(xfuncs, center = TRUE, scale = FALSE)
  else lm(xfuncs ~ X0 - 1)$resid
  SB <- sigs.decor %*% basismat
  return(list(X0 = X0, SB = SB, penmat = penmat, basismat = basismat, xfuncs = xfuncs,
              nbasis = nbasis, argvals = argvals, dim.sig = dim.sig))
}

#' interim function for @function FPLS_R

fpls = function (y, xfuncs = NULL, fdobj = NULL, ncomp = NULL, pve = 0.99, 
                 nbasis = NULL, basismat = NULL, penmat = NULL, argvals = NULL, 
                 covt = NULL, mean.signal.term = FALSE, spline.order = NULL, 
                 family = "gaussian", method = "REML", sp = NULL, pen.order = 2, 
                 cv1 = FALSE, nfold = 5, store.cv = FALSE, store.gam = TRUE) {
  
  if (!("MASS" %in% rownames(installed.packages()))) install.packages("MASS")
  
  n <- length(y)
  do.cv <- FALSE
  if (!nfold %in% 1:n) 
    stop("Argument 'nfold' is invalid: must be an integer between 1 and ", 
         n, ".")
  if (!family %in% c("gaussian", "binomial")) 
    stop("Only 'gaussian' and 'binomial' models are implemented in the current version.")
  if (is.null(fdobj)) {
    if (!is.array(xfuncs) || !length(dim(xfuncs)) %in% 2:3) 
      stop("xfuncs must either be a 2D or 3D array")
    dim.sig <- length(dim(xfuncs)) - 1
    if (is.null(nbasis)) 
      nbasis <- ifelse(dim.sig == 1, 40, 15)
    if (dim.sig == 1) {
      if (!is.list(nbasis) && is.vector(nbasis)) 
        nbs <- matrix(nbasis, ncol = 1)
      else stop("for 1D predictors, 'nbasis' must be a vector")
      siglength <- ncol(xfuncs)
    }
    else {
      if (is.list(nbasis) && length(nbasis) == 2) 
        nbs <- cbind(rep(nbasis[[1]], length(nbasis[[2]])), 
                     rep(nbasis[[2]], each = length(nbasis[[1]])))
      else if (!is.list(nbasis) && is.vector(nbasis)) 
        nbs <- matrix(rep(nbasis, 2), ncol = 2)
      else stop("for 2D predictors, 'nbasis' must either be a vector or a list of length 2")
      d1 <- dim(xfuncs)[2L]
      d2 <- dim(xfuncs)[3L]
      siglength <- d1 * d2
    }
    if (cv1 || nrow(nbs) > 2 || (!is.null(ncomp) && length(ncomp) > 
                                 1)) 
      do.cv <- TRUE
  }
  else if (cv1 || (!is.null(ncomp) && length(ncomp) > 1)) 
    do.cv <- TRUE
  if (!do.cv) {
    store.cv <- FALSE
    nfold <- 1
  }
  groups <- split(sample(1:n), rep(1:nfold, length = n))
  n.unpen.cols <- 1 + mean.signal.term + ifelse(is.null(covt), 
                                                0, ncol(as.matrix(covt)))
  cv.table <- array(0, dim = c(ifelse(is.null(fdobj), nrow(nbs), 
                                      1), ifelse(is.null(ncomp), 1, length(ncomp))))
  for (inb in 1:dim(cv.table)[1]) {
    st <- fpcr.setup(y = y, xfuncs = xfuncs, fdobj = fdobj, 
                     nbasis = if (is.null(fdobj)) 
                       nbs[inb, ]
                     else NULL, basismat = basismat, penmat = penmat, 
                     argvals = argvals, covt = covt, mean.signal.term = mean.signal.term, 
                     spline.order = spline.order, pen.order = pen.order)
    argvals <- st$argvals

    for (ifold in 1:nfold) {
      if (do.cv) {
        idxTest <- groups[[ifold]]
        idxTrain <- (1:n)[-idxTest]
      }
      else idxTest <- idxTrain <- 1:n
      X0.tr <- st$X0[idxTrain, ]
      SB.tr <- st$SB[idxTrain, ]
      svdSB <- svd(SB.tr)      
      if (is.null(ncomp)) 
        ncomp <- min(which(cumsum(svdSB$d) > pve * sum(svdSB$d)))
      
      R <- simpls(y = y, x = SB.tr, ncomp = min(length(y)-1, ncol(SB.tr)-1, max(ncomp)))$R   
      
      for (incomp in 1:length(ncomp)) {
        R.ncomp <- R[, 1:ncomp[incomp]]
        X <- cbind(X0.tr, SB.tr %*% R.ncomp)
        S <- list(matrix(0, ncol(X), ncol(X)))
        S[[1]][-(1:n.unpen.cols), -(1:n.unpen.cols)] <- crossprod(R.ncomp, 
                                                                  st$penmat %*% R.ncomp)
        obje <- mgcv::gam(y[idxTrain] ~ X - 1, paraPen = list(X = S), 
                    family = get(family), method = method, sp = sp)
        BV <- st$basismat %*% R.ncomp
        fhat <- BV %*% obje$coef[-(1:n.unpen.cols)]
        undecor.coef <- obje$coef[1:n.unpen.cols] - MASS::ginv(X0.tr) %*% 
          st$xfuncs[idxTrain, ] %*% fhat
        yhat <- st$X0[idxTest, ] %*% undecor.coef + st$xfuncs[idxTest, 
                                                              ] %*% fhat
        if (family == "gaussian") 
          cv.table[inb, incomp] <- cv.table[inb, incomp] + 
          mean((yhat - y[idxTest])^2)
        else if (family == "binomial") {
          phat <- exp(yhat)/(1 + exp(yhat))
          phat <- replace(phat, exp(yhat) == Inf, 1)
          cv.table[inb, incomp] <- cv.table[inb, incomp] + 
            mean((phat > mean(y[idxTrain])) != y[idxTest])
        }
      }
    }
  }
  if (do.cv) {
    idxmin <- which(cv.table == min(cv.table[cv.table != 
                                               0], na.rm = TRUE), arr.ind = TRUE)
    if (nrow(idxmin) > 1) 
      idxmin <- idxmin[1, ]
    if (is.list(nbasis)) {
      dim(cv.table) <- c(length(nbasis[[1]]), length(nbasis[[2]]), 
                         length(ncomp))
      dimnames(cv.table) <- list(paste("nbasis1=", nbasis[[1]]), 
                                 paste("nbasis2=", nbasis[[2]]), paste("ncomp", 
                                                                       ncomp))
    }
    else dimnames(cv.table) <- list(paste("nbasis=", nbasis), 
                                    paste("ncomp", ncomp))
    if (dim.sig == 1) 
      nbasis <- nbs[idxmin[1], ]
    else {
      nbasis <- list()
      nbasis[[1]] <- nbs[idxmin[1], 1]
      nbasis[[2]] <- nbs[idxmin[1], 2]
    }
    obje <- fpls(y = y, xfuncs = xfuncs, fdobj = fdobj, ncomp = ncomp[idxmin[2]], 
                 nbasis = nbasis, basismat = basismat, penmat = penmat, 
                 argvals = argvals, covt = covt, mean.signal.term = mean.signal.term, 
                 spline.order = spline.order, family = family, method = method, 
                 sp = sp, pen.order = pen.order, cv1 = FALSE, store.gam = store.gam)
    ncomp <- ncomp[idxmin[2]]
    if (store.cv) 
      obje$cv.table <- cv.table
    else obje$cv.table <- min(cv.table[cv.table != 0])
  }
  else {
    se <- sqrt(rowSums((BV %*% obje$Vp[-(1:n.unpen.cols), 
                                       -(1:n.unpen.cols)]) * BV))
    if (!store.gam) {
      yhat <- obje$fitted.values
      obje <- list()
      obje$fitted.values <- yhat
    }
    obje$se <- se
    obje$argvals <- argvals
    obje$undecor.coef <- as.matrix(undecor.coef, nrow = 1)
    colnames(obje$undecor.coef) <- if (is.null(dimnames(covt)) || 
                                       is.null(dimnames(covt)[[2]])) 
      paste("X", 0:(n.unpen.cols - 1), sep = "")
    else c("Intercept", dimnames(covt)[[2]])
    if (st$dim.sig == 2) 
      dim(fhat) <- c(d1, d2)
    obje$fhat <- fhat
    obje$nbasis <- nbasis
    obje$ncomp <- ncomp
  }
  class(obje) = "fpcr"
  return(obje)
}

#' implement FCR with 1-dim functional x and 1-dim y.
#' @param y.old is the training response vector of length N
#' @param x.old is the training predictor matrix of N*m
#' @param y.new is the training response vector
#' @param x.new is the training predictor matrix
#' @param hyperPars is of two columns: 1st column for alpha and 2nd column for upper bound of p

FCR = function(y.old, x.old, y.new = NULL, x.new = NULL, hyperPars){
  
  Alpha = sort(hyperPars[, 1])
  pMax = hyperPars[, 2][order(hyperPars[, 1])]
  
  x.train = sweep(x.old, 2, colMeans(x.old)) # centering x
  y.train = y.old - mean(y.old) # centering y
  N = length(y.train)
  
  expansion = expand.bspline(x.train)
  splineBasis = expansion$basis
  basismatrix = eval.basis(1:ncol(x.old), expansion$basis)
  C = expansion$coef
  W = inprod(splineBasis, splineBasis)
  
  if (!is.null(x.new)){
    x.test = sweep(x.new, 2, colMeans(x.old))
    fd.test = Data2fd(y=t(x.test), argvals=1:ncol(x.new), expansion$D2fdPar, nderiv=2, lambda=expansion$lambda)
    C.star = t(fd.test$coef)
  }

  betaCoef = array(NA, dim=c(length(Alpha), max(pMax), ncol(W)))
  Beta = array(NA, dim=c(length(Alpha), max(pMax), ncol(x.old)))
  GCV = matrix(NA, nrow=length(Alpha), ncol=max(pMax))
  CV = matrix(NA, nrow=length(Alpha), ncol=max(pMax))
  
  for(j in 1:length(Alpha)){
    
    alpha = Alpha[j]
    gam = alpha / (1-alpha)
    
    b = list()
    G = list()
    H = matrix(NA, nrow = N, ncol = pMax[j]) # initiate H
    wCoef = matrix(NA, nrow = ncol(W), ncol = pMax[j])
    
    Whalf = half(W)
    Whalfinv = halfinv(W)
    CWhalf = C %*% half(W)
    PCWhalf = CWhalf
    r = qr(CWhalf)$rank
    V = list()
    
    for(i in 1:pMax[j]){
      
      if (i==1){
        P = diag(N)
      }
      else{
        P = P %*% (diag(N) - H[, i-1] %*% crossprod(H[, i-1])^{-1} %*% t(H[, i-1]))
      }
      
      PCWhalf = P %*% PCWhalf
      svdPCWhalf = svd(PCWhalf)
      E = diag(svdPCWhalf$d[1:(r-i+1)])
      V[[i]] = svdPCWhalf$v[, 1:(r-i+1)]
      G[[i]] = svdPCWhalf$u[, 1:(r-i+1)] %*% E
      Gy = crossprod(G[[i]], y.train)
      GG = E^2
      zeta = max(GG)
      
      L = function(delta) {
        return(solve(GG + delta^-1 * zeta * diag(r-i+1)))
      }

      neglnQ = function(delta) {
        Ldelta = L(delta)
        result = -as.numeric(2 * log(abs(t(Gy) %*% Ldelta %*% Gy)) -
                               gam * log(abs(t(Gy) %*% Ldelta %*% Ldelta %*% Gy)) +
                               (gam-1) * log(abs(t(Gy) %*% Ldelta %*% GG %*% Ldelta %*% Gy)))
        if (is.infinite(result)) return(1e10)
        else return(result)
      }

      deltaOpt = manymodes(neglnQ, ninitial = ninitial)
      LdeltaOpt = L(deltaOpt)
      b[[i]] = LdeltaOpt %*% Gy / as.numeric(crossprod(LdeltaOpt %*% Gy))^.5
      
      if (length(b)<i) break
      H[, i] = CWhalf %*% V[[i]] %*% b[[i]]
      wCoef[, i] = Whalfinv %*% V[[i]] %*% b[[i]]
      
      betaCoef[j ,i, ] = as.matrix(wCoef[, 1:i]) %*% lm.fit(x = as.matrix(H[, 1:i]), y = y.train)$coefficients
      residual = y.train - C %*% W %*% betaCoef[j,i,]
      GCV[j,i] = sum(residual^2)/(N-i-1)^2
      Beta[j,i,] = basismatrix %*% betaCoef[j,i,]
      
      if (!is.null(x.new)) {
        pred = mean(y.old) + C.star %*% W %*% betaCoef[j,i,]
        if (!is.null(y.new)){
          CV[j,i] = sum((y.new-pred)^2)
        }
      }
    }
  }
  
  if (nrow(GCV)==1) {
    criterion.value = matrix(GCV[, 1:max(pMax)], nrow=1)
  }else {
    criterion.value = GCV[, 1:max(pMax)]
  }
  loc.opt = which(criterion.value == min(criterion.value, na.rm=TRUE), arr.ind=TRUE)
  
  sse = NULL
  pars.opt = NULL
  beta.opt = matrix(NA, nrow=nrow(loc.opt), ncol=ncol(x.old))
  
  for (i in 1:nrow(loc.opt)){ # in case there are multi optimizers
    pars.opt = rbind(pars.opt, c(Alpha[loc.opt[i,1]], loc.opt[i,2]))
    beta.opt[i,] = Beta[loc.opt[i,1], loc.opt[i,2], ]
  }
  colnames(pars.opt) = c("alpha", "ncomp")
  
  if (!is.null(x.new) & !is.null(y.new)){
    for (i in 1:nrow(loc.opt)){ # in case there are multi optimizers
      if (i==1) {
        sse = CV[loc.opt[i,1], loc.opt[i,2]]
      }else {
        sse = c(sse, CV[loc.opt[i,1], loc.opt[i,2]])
      } 
    }
  }
  
  return(list(C=C, W=W, CV=CV, GCV=GCV, D2fdPar=expansion$D2fdPar, pars.opt=pars.opt, sse=sse, beta=beta.opt))
}

#' implement pFPLS in Aguilera et al. (2016, Chemometr. Intell. Lab. vol. 154, pp. 80-92).

pFPLS = function(y.old, x.old, y.new=NULL, x.new=NULL, Theta.pFPLS, pMax, pMin){

  N = length(y.old)
  nfold = 10 # 10-fold 
  test = split(1:N, 1:nfold) 
  CV = matrix(NA, nrow=length(Theta.pFPLS), ncol=pMax-pMin+1)
  
  for (k in 1:nfold) {
    
    x.train = sweep(x.old[-test[[k]],], 2, colMeans(x.old[-test[[k]],]))
    y.train = as.matrix(y.old[-test[[k]]] - mean(y.old[-test[[k]]]))
    x.test = sweep(x.old[test[[k]],], 2, colMeans(x.old[-test[[k]],]))
    y.test = as.matrix(y.old[test[[k]]])
    
    expansion = expand.bspline(x.train)
    splineBasis = expansion$basis
    C = expansion$coef
    W = inprod(splineBasis, splineBasis)
    G = inprod(splineBasis, splineBasis, Lfdobj1=2, Lfdobj2=2)
    
    fd.test = Data2fd(y=t(x.test), argvals=1:ncol(x.test), expansion$D2fdPar, nderiv=2, lambda=expansion$lambda)
    C.star = t(fd.test$coef)
    
    CV.tmp = matrix(NA, nrow=length(Theta.pFPLS), ncol=pMax-pMin+1)
  
    for (j in 1:length(Theta.pFPLS)) {
      theta = Theta.pFPLS[j]
      Wtheta = C %*% W %*% halfinv(W + theta*G)
      Wtheta.star = C.star %*% W %*% halfinv(W + theta*G)
    
      for (i in pMin:pMax){
        simpls.coef = simpls(y=y.train, x=Wtheta, ncomp=i)$coef
        pred = mean(y.old[-test[[k]]]) + Wtheta.star %*% simpls.coef
        CV.tmp[j, i-pMin+1] = sum((y.test - pred)^2)
      }
    }
    
    if (k==1) CV = CV.tmp
    else CV = CV + CV.tmp
  }
  
  pred = NULL
  sse = NULL
  pars.opt = NULL
  beta.opt = NULL
  
  criterion.value = CV
  loc.opt = which(criterion.value == min(criterion.value, na.rm=TRUE), arr.ind=TRUE)
  
  x.train = sweep(x.old, 2, colMeans(x.old))
  y.train = as.matrix(y.old - mean(y.old))
  N = length(y.train)
  
  expansion = expand.bspline(x.train)
  splineBasis = expansion$basis
  C = expansion$coef
  W = inprod(splineBasis, splineBasis)
  G = inprod(splineBasis, splineBasis, Lfdobj1=2, Lfdobj2=2)
  basismatrix = eval.basis(1:ncol(x.old), expansion$basis)
  
  if (is.null(x.new)) {
    for (i in 1:nrow(loc.opt)){ # in case there are multi optimizers
      theta = Theta.pFPLS[loc.opt[i,1]]
      p = loc.opt[i,2]+pMin-1
      pars.opt = rbind(pars.opt, c(theta, p))
      
      Wtheta = C %*% W %*% halfinv(W + theta*G)
      simpls.coef = simpls(y=y.train, x=Wtheta, ncomp=p)$coef
      beta.opt = rbind(beta.opt, matrix(basismatrix %*% halfinv(W + theta*G) %*% simpls.coef, nrow=1))
    }
  } else {
    
    x.test = sweep(x.new, 2, colMeans(x.old))
    fd.test = Data2fd(y=t(x.test), argvals=1:ncol(x.test), expansion$D2fdPar, nderiv=2, lambda=expansion$lambda)
    C.star = t(fd.test$coef)
    
    for (i in 1:nrow(loc.opt)){ # in case there are multi optimizers
      
      theta = Theta.pFPLS[loc.opt[i,1]]
      p = loc.opt[i,2]+pMin-1
      pars.opt = rbind(pars.opt, c(theta, p))
      
      Wtheta = C %*% W %*% halfinv(W + theta*G)
      Wtheta.star = C.star %*% W %*% halfinv(W + theta*G)
      simpls.coef = simpls(y=y.train, x=Wtheta, ncomp=p)$coef
      beta.opt = rbind(beta.opt, matrix(basismatrix %*% halfinv(W + theta*G) %*% simpls.coef, nrow=1))
      pred.tmp = mean(y.old) + Wtheta.star %*% simpls.coef
      pred = cbind(pred, pred.tmp)
      if (!is.null(y.new)) sse = cbind(sse, sum((pred.tmp-y.new)^2))
    }
    
    colnames(pars.opt) = c("theta", "ncomp")
  }
  
  return(list(C=C, W=W, G=G, 
              pMax=pMax, pars.opt=pars.opt,
              CV=CV, D2fdPar=expansion$D2fdPar,
              pred=pred, sse=sse, beta=beta.opt))
}

#' implement FPCR_R-REML in Reiss & Ogden (2007, JASA, vol. 102, pp. 984-996)

FPCR_R = function(y.old, x.old, y.new=NULL, x.new=NULL){
  
  if (!("refund" %in% rownames(installed.packages()))) install.packages("refund")
  if (!all(c(
    exists("fpcr", mode="function")
  ))) library("refund")
  
  x.train = sweep(x.old, 2, colMeans(x.old)) # centering x
  y.train = y.old - mean(y.old) # centering y
  expansion = expand.bspline(x.train)
  
  model = fpcr(y=y.train, xfuncs=x.train, nbasis=ncol(expansion$coef))
  CoefFun = model$fhat
  
  if (is.null(x.new)){
    pred = NULL
    sse = NULL
  }else{
    x.test = sweep(x.new, 2, colMeans(x.old))
    pred = mean(y.old) + x.test %*% CoefFun
    if (is.null(y.new)){
      sse = NULL
    }
    sse = sum((pred-y.new)^2)
  }
  return(list(pred=pred, sse=sse, beta=CoefFun))
}

#' implement FPLS_R-REML in Reiss & Ogden (2007, JASA, vol. 102, pp. 984-996)

FPLS_R = function(y.old, x.old, y.new=NULL, x.new=NULL){
  
  if (!("refund" %in% rownames(installed.packages()))) install.packages("refund")
  library("refund")
  
  x.train = sweep(x.old, 2, colMeans(x.old)) # centering x
  y.train = as.matrix(y.old - mean(y.old)) # centering y
  expansion = expand.bspline(x.train)
  
  model = fpls(y=y.train, xfuncs=x.train, nbasis=ncol(expansion$coef))
  CoefFun = model$fhat
  
  if (is.null(x.new)){
    pred = NULL
    sse = NULL
  }else{
    x.test = sweep(x.new, 2, colMeans(x.old))
    pred = mean(y.old) + x.test %*% CoefFun
    if (is.null(y.new)){
      sse = NULL
    }else{
      sse = sum((pred-y.new)^2)
    }
  }

  return(list(pred=pred, sse=sse, beta=CoefFun))
}

#' implement sFPCR 

sFPCR = function(y.old, x.old, y.new=NULL, x.new=NULL, Alpha, Theta, pMax=5, pMin=1){

  x.train = sweep(x.old, 2, colMeans(x.old)) # centering x
  y.train = y.old - mean(y.old) # centering y
  N = length(y.train)
  parsGrid = expand.grid(Alpha, Theta) #parameter grid of alpha & theta

  expansion = expand.bspline(x.train)
  splineBasis = expansion$basis
  basismatrix = eval.basis(1:ncol(x.train), expansion$basis)
  S = expansion$coef
  W = inprod(splineBasis, splineBasis)
  D = inprod(splineBasis, splineBasis, Lfdobj1=2, Lfdobj2=2)
  M = W %*% t(S) %*% y.train

  Beta = list()
  GCV = matrix(NA, nrow=dim(parsGrid)[1], ncol=pMax) # initiation

  for(j in 1:dim(parsGrid)[1]){

    theta = parsGrid[j,1]
    lambda = parsGrid[j,2]

    G = W + lambda*D
    U = theta/N * crossprod(S %*% W) + (1-theta)/N^2 * tcrossprod(M)

    G.halfinv = halfinv(G)

    Delta = svd(crossprod(G.halfinv, U %*% G.halfinv), nu=pMax, nv=pMax)$u
    Beta[[j]] = G.halfinv %*% Delta

    for (i in 1:pMax){
      model = lm.fit(x=S %*% W %*% Beta[[j]][,1:i], y=y.train)
      residual = model$residuals
      GCV[j,i] = N*sum(residual^2)/(N-i-1)^2
    }
  }

  pred = NULL
  sse = NULL
  pars.opt = NULL
  beta.opt = NULL
  
  if (nrow(GCV)==1) criterion.value = matrix(GCV[, pMin:pMax], nrow=1)
  else criterion.value = GCV[, pMin:pMax]
  loc.opt = which(criterion.value == min(criterion.value, na.rm=TRUE), arr.ind=TRUE)

  if (is.null(x.new)) {
    for (i in 1:nrow(loc.opt)){ # in case there are multi optimizers
      pars.opt = rbind(pars.opt, c(parsGrid[loc.opt[i,1],], loc.opt[i,2]+pMin-1))
      Beta.tmp = Beta[[loc.opt[i,1]]][, 1:(loc.opt[i,2]+pMin-1)]
      model.tmp = lm.fit(x=S %*% W %*% Beta.tmp, y=y.train)
      beta.opt = rbind(beta.opt, matrix(basismatrix %*% Beta.tmp %*% (model.tmp$coefficients), nrow=1))
    }
  } else {
    x.test = sweep(x.new, 2, colMeans(x.old))
    fd.test = Data2fd(y=t(x.test), argvals=1:ncol(x.new),
                      expansion$D2fdPar, nderiv=2, lambda=expansion$lambda)
    S.star = t(fd.test$coef)

    for (i in 1:nrow(loc.opt)){ # in case there are multi optimizers

      pars.opt = rbind(pars.opt, c(parsGrid[loc.opt[i,1], ], loc.opt[i,2]+pMin-1))
      Beta.tmp = Beta[[loc.opt[i,1]]][, 1:(loc.opt[i,2]+pMin-1)]
      model.tmp = lm.fit(x=S %*% W %*% Beta.tmp, y=y.train)
      beta.opt = rbind(beta.opt, matrix(basismatrix %*% Beta.tmp %*% (model.tmp$coefficients), nrow=1))
      pred.tmp = mean(y.old) + S.star %*% W %*% Beta.tmp %*% (model.tmp$coefficients)
      pred = cbind(pred, pred.tmp)
      if (!is.null(y.new)) sse = cbind(sse, sum((pred.tmp-y.new)^2))
    }

    colnames(pars.opt) = c("theta", "lambda", "ncomp")
  }

  return(list(S=S, W=W, G=G, Beta=Beta,
              parsGrid=parsGrid, pMax=pMax,
              GCV=GCV, D2fdPar=expansion$D2fdPar,
              pred=pred, pars.opt=pars.opt, sse=sse, beta=beta.opt))

}

#' integrate ReMSPEs resulted from different methods

integrate.remspe = function(ReMSPE.FCR.random 
                            ,ReMSPE.FCR.grid
                            ,ReMSPE.sFPCR
                            ,ReMSPE.pFPLS 
                            ,ReMSPE.FPLSR
                            ,ReMSPE.FPCRR
                            ,ReMSPE.FPCR
                            ){
  nrow.ReMSPE = min(
    length(ReMSPE.FCR.random)
    ,length(ReMSPE.FCR.grid)
    ,length(ReMSPE.sFPCR)
    ,length(ReMSPE.pFPLS)
    ,length(ReMSPE.FPLSR)
    ,length(ReMSPE.FPCRR)
    ,length(ReMSPE.FPCR)
  )
  if (nrow.ReMSPE>0){
    ReMSPE = cbind(
      ReMSPE.FCR.random[1:nrow.ReMSPE]
      ,ReMSPE.FCR.grid[1:nrow.ReMSPE]
      ,ReMSPE.sFPCR[1:nrow.ReMSPE]
      ,ReMSPE.pFPLS[1:nrow.ReMSPE]
      ,ReMSPE.FPLSR[1:nrow.ReMSPE]
      ,ReMSPE.FPCRR[1:nrow.ReMSPE]
      ,ReMSPE.FPCR[1:nrow.ReMSPE]
    )
    colnames(ReMSPE) = c(
      'fun CR random'
      ,'fun CR grid'
      ,'sup FPCA'
      ,'pFPLS'
      ,'FPLS_R-REML'
      ,'FPCR_R-REML'
      ,'sm fun PCA'
    )
  }
  return(ReMSPE)
}

#' integrate RMSE curves resulted from different methods

integrate.rmse = function(betaRMSE.FCR.random
                          ,betaRMSE.FCR.grid
                          # ,betaRMSE.sFPCR
                          ,betaRMSE.pFPLS
                          ,betaRMSE.FPLSR
                          ,betaRMSE.FPCRR
                          ,betaRMSE.FPCR
                          ){
  nrow.RMSE = min(
    length(betaRMSE.FCR.random)
    ,length(betaRMSE.FCR.grid)
    # ,length(betaRMSE.sFPCR)
    ,length(betaRMSE.pFPLS)
    ,length(betaRMSE.FPLSR)
    ,length(betaRMSE.FPCRR)
    ,length(betaRMSE.FPCR)
  )
  if (nrow.RMSE>0){
    RMSE = cbind(
      betaRMSE.FCR.random[1:nrow.RMSE]
      ,betaRMSE.FCR.grid[1:nrow.RMSE]
      # ,betaRMSE.sFPCR[1:nrow.RMSE]
      ,betaRMSE.pFPLS[1:nrow.RMSE]
      ,betaRMSE.FPLSR[1:nrow.RMSE]
      ,betaRMSE.FPCRR[1:nrow.RMSE]
      ,betaRMSE.FPCR[1:nrow.RMSE]
    )
    colnames(RMSE) = c(
      'fun CR random'
      ,'fun CR grid'
      # ,'sup FPCA'
      ,'pFPLS'
      ,'FPLS_R-REML'
      ,'FPCR_R-REML'
      ,'sm fun PCA'
    )
  }
  return(RMSE)
}

#' boxplot of ReMSPE

boxplot.remspe = function(nrowplot, ReMSPE){
  
  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  library(ggplot2)
  library(reshape2)
  
  melton = melt(data.frame(ReMSPE[1:nrowplot,], Replication=1:nrow(ReMSPE[1:nrowplot,])), id.vars = "Replication")
  bplot = ggplot(melton, aes(x = variable, y = value, colour = variable)) + 
    geom_boxplot(width=.5, outlier.shape=NA) +
    coord_cartesian(ylim = c(0, switch(case.no, .05, .15, .2, .5))) +
    labs(x = '', y="ReMSPE") + 
    theme_bw() +
    theme(legend.position = "none", panel.border = element_blank(), plot.title = element_text(hjust = 0.5)) 
  
  plot(bplot)
  ggsave(file = paste('figure/real_',
                      switch(case.no, 'medfly_', 'fat_', 'moisture_', 'protein_'), 
                      'seed=', seed, '_',
                      'loop=', RR, '_',
                      'ninitial=', ninitial, '_600dpi.pdf',
                      sep=''), 
         dpi=600)
}

# curve plot of RMSE

curveplot.rmse = function(RMSE){
  
  if (!("ggplot2" %in% rownames(installed.packages()))) install.packages("ggplot2")
  library(ggplot2)
  library(reshape2)
  
  melton = melt(data.frame(RMSE, Day = denseGrid), id.vars="Day", variable.name="Method", value.name="RMSE")
  curveplot = ggplot(melton, aes(x=Day, y=RMSE, color=Method, linetype=Method)) + 
    geom_line() +
    # coord_cartesian(ylim = c(0,  .1)) +
    labs(x=expression(t), y="RMSE") +
    theme_bw() +
    theme(legend.position=c(0.8, 0.7), panel.border = element_blank(), plot.title = element_text(hjust = 0.5))
  
  plot(curveplot)
  ggsave(file = paste('figure/simul_',
                      'scene=', scene, '_',
                      switch(distri, 'normal_', 'exp_'),
                      'seed=', seed, '_',
                      'loop=', RR, '_',
                      'ninitial=', ninitial, '_',
                      'snr=', snr, '.pdf',
                      sep=''), 
         dpi=600)
}