#' Put this r file into the folder containing functions.R
#' In the same folder, create two subfolders named 'figure' & 'Rimage', respectively.

if (!("rstudioapi" %in% rownames(installed.packages()))) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#### global setting
source("functions.R")
seed = 1
set.seed(seed)
RR = 500 # number of replication
N = 100 # number of curves
ninitial = 1 # number of initial values used in optimization
snr = 4 # signal-to-noise ratio used in simulation
scene = 2 # 1 for the scene with beta = 1st eigenfunction; 2 for the scene with beta = 3rd eigenfunction
distri = 2 # 1 for N(0,1); 2 for exp(1)

### candidate pools for parameters to be tuned
pUpper = 2 # max number of components
pLower = 1 # min number of components
Theta.pFPLS = c(0, 1, 1e1, 1e2, 1e3, 1e4, 1e5)
Alpha.sFPCR = seq(from=0, to=1, len=11)
Theta.sFPCR = c(0, 1, 1e1, 1e2, 1e3, 1e4, 1e5)

#### data generation

## simulation scene I: beta = 1st eigenfunction 

if (scene == 1) {
  if (!("fda" %in% rownames(installed.packages()))) 
    install.packages("fda")
  library(fda)
  data(CanadianWeather)
  x = t(CanadianWeather$dailyAv[,,3])
  denseGrid = 1:ncol(x)
  expansion.x = expand.bspline(x)
  W = inprod(expansion.x$basis, expansion.x$basis)
  basismatrix = eval.basis(denseGrid, expansion.x$basis)
  
  x.fd = smooth.basisPar(argvals=denseGrid,
                         y=t(x), fdobj=expansion.x$basis,
                         Lfdobj=int2Lfd(2), lambda=expansion.x$lambda)
  x.pca.obj = pca.fd(x.fd$fd, nharm = 6)
  x.eigenfun = t(basismatrix %*% x.pca.obj$harmonics$coefs[,1:3])
  beta.true = x.eigenfun[1,]
  
  mu.x = apply(x, 2, mean)
  mu.x.fd = Data2fd(y=mu.x, argvals=denseGrid, expansion.x$D2fdPar, nderiv=2, lambda=expansion.x$lambda)
  mu.y = t(mu.x.fd$coefs) %*% W %*% matrix(x.pca.obj$harmonics$coefs[,1], ncol=1)
  
  x.train = list()
  y.train = list()
  for (R in 1:RR){
    x.train[[R]] = matrix(NA, nrow=N, ncol=ncol(x))
    y.train[[R]] = numeric(N)
    for (n in 1:N){
      score.tmp = rnorm(n=3, mean=0, sd=sqrt(x.pca.obj$values[1:3]))
      x.train[[R]][n,] = matrix(score.tmp, nrow=1) %*% x.eigenfun + mu.x
      y.train[[R]][n] = score.tmp[1] + mu.y
    }
    y.train[[R]] = y.train[[R]] + ((x.pca.obj$values[1])/snr)^.5 * switch(distri,
      rnorm(n = N),
      rexp(n = N) - 1
    )
  }
}

## simulation scene II: beta = 3rd eigenfunction

if (scene == 2) {
  if (!("fda" %in% rownames(installed.packages()))) 
    install.packages("fda")
  library(fda)
  data(CanadianWeather)
  x = t(CanadianWeather$dailyAv[,,3])
  denseGrid = 1:ncol(x)
  expansion.x = expand.bspline(x)
  W = inprod(expansion.x$basis, expansion.x$basis)
  basismatrix = eval.basis(denseGrid, expansion.x$basis)
  
  x.fd = smooth.basisPar(argvals=denseGrid,
                         y=t(x), fdobj=expansion.x$basis,
                         Lfdobj=int2Lfd(2), lambda=expansion.x$lambda)
  x.pca.obj = pca.fd(x.fd$fd, nharm = 6)
  cumsum(x.pca.obj$varprop[1:3])
  x.eigenfun = t(basismatrix %*% x.pca.obj$harmonics$coefs[,1:3])
  beta.true = x.eigenfun[3,]
  
  mu.x = apply(x, 2, mean)
  mu.x.fd = Data2fd(y=mu.x, argvals=denseGrid, expansion.x$D2fdPar, nderiv=2, lambda=expansion.x$lambda)
  mu.y = t(mu.x.fd$coefs) %*% W %*% matrix(x.pca.obj$harmonics$coefs[,3], ncol=1)
  
  x.train = list()
  y.train = list()
  for (R in 1:RR){
    x.train[[R]] = matrix(NA, nrow=N, ncol=ncol(x))
    y.train[[R]] = numeric(N)
    for (n in 1:N){
      score.tmp = rnorm(n=3, mean=0, sd=sqrt(x.pca.obj$values[1:3]))
      x.train[[R]][n,] = matrix(score.tmp, nrow=1) %*% x.eigenfun + mu.x
      y.train[[R]][n] = score.tmp[3] + mu.y
    }
    y.train[[R]] = y.train[[R]] +  + ((x.pca.obj$values[3])/snr)^.5 * switch(distri,
      rnorm(n = N, mean = 0, sd = 1),
      rexp(n = N) - 1
    )
  }
}

# Dangerous! erase existing dat of squared error of etimated beta
betaSE.FCR.random = NULL
betaSE.FCR.grid = NULL
betaSE.pFPLS = NULL
betaSE.FPCRR = NULL
betaSE.FPLSR = NULL
betaSE.sFPCR = NULL
betaSE.FPCR = NULL

check.FCR.random = ifelse(is.null(betaSE.FCR.random), 0, nrow(betaSE.FCR.random))
check.FCR.grid = ifelse(is.null(betaSE.FCR.grid), 0, nrow(betaSE.FCR.grid))
check.pFPLS = ifelse(is.null(betaSE.pFPLS), 0, nrow(betaSE.pFPLS))
check.FPCRR = ifelse(is.null(betaSE.FPCRR), 0, nrow(betaSE.FPCRR))
check.FPLSR = ifelse(is.null(betaSE.FPLSR), 0, nrow(betaSE.FPLSR))
check.sFPCR = ifelse(is.null(betaSE.sFPCR), 0, nrow(betaSE.sFPCR))
check.FPCR = ifelse(is.null(betaSE.FPCR), 0, nrow(betaSE.FPCR))

################ Replica & Time ##################

ptm = NULL
ptm[1] = proc.time()[3]

for (R in 1:RR){
  if (R > check.FCR.random){
    hyperPars.FCR.random = runif(10, min = 0, max = 1)
    hyperPars.FCR.random = cbind(hyperPars.FCR.random, sample.int(pUpper, size=length(hyperPars.FCR.random), replace=T))
    FCR.random.result = FCR(y.old=y.train[[R]], x.old=x.train[[R]], hyperPars=hyperPars.FCR.random)
    betaSE.FCR.random = rbind(betaSE.FCR.random, apply(sweep(FCR.random.result$beta, 2, beta.true)^2, 2, mean))
    if (R==RR) betaRMSE.FCR.random = sqrt(apply(betaSE.FCR.random, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[2] = proc.time()[3]

for (R in 1:RR){
  if (R > check.sFPCR){
    sFPCR.result = sFPCR(y.old=y.train[[R]], x.old=x.train[[R]], Alpha=Alpha.sFPCR, Theta=Theta.sFPCR, pMax=pUpper, pMin=pLower)
    betaSE.sFPCR = rbind(betaSE.sFPCR, apply(sweep(sFPCR.result$beta, 2, beta.true)^2, 2, mean))
    if (R==RR) betaRMSE.sFPCR = sqrt(apply(betaSE.sFPCR, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[3] = proc.time()[3]

for (R in 1:RR){
  if (R > check.pFPLS){
    pFPLS.result = pFPLS(y.old=y.train[[R]], x.old=x.train[[R]], Theta.pFPLS=Theta.pFPLS, pMax=pUpper, pMin=pLower)
    betaSE.pFPLS = rbind(betaSE.pFPLS, apply(sweep(pFPLS.result$beta, 2, beta.true)^2, 2, mean))
    if (R==RR) betaRMSE.pFPLS = sqrt(apply(betaSE.pFPLS, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[4] = proc.time()[3]

for (R in 1:RR){
  if (R > check.FPLSR){
    FPLSR.result = FPLS_R(y.old=y.train[[R]], x.old=x.train[[R]])
    betaSE.FPLSR = rbind(betaSE.FPLSR, sweep(matrix(FPLSR.result$beta, nrow=1), 2, beta.true)^2)
    if (R==RR) betaRMSE.FPLSR = sqrt(apply(betaSE.FPLSR, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[5] = proc.time()[3]

for (R in 1:RR){
  if (R > check.FPCRR){
    FPCRR.result = FPCR_R(y.old=y.train[[R]], x.old=x.train[[R]])
    betaSE.FPCRR = rbind(betaSE.FPCRR, sweep(matrix(FPCRR.result$beta, nrow=1), 2, beta.true)^2)
    if (R==RR) betaRMSE.FPCRR = sqrt(apply(betaSE.FPCRR, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[6] = proc.time()[3]

for (R in 1:RR){
  if (R > check.FPCR){
    FPCR.result = sFPCR(y.old=y.train[[R]], x.old=x.train[[R]], Alpha=1, Theta=Theta.sFPCR, pMax=pUpper, pMin=pLower)
    betaSE.FPCR = rbind(betaSE.FPCR, apply(sweep(FPCR.result$beta, 2, beta.true)^2, 2, mean))
    if (R==RR) betaRMSE.FPCR = sqrt(apply(betaSE.FPCR, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[7] = proc.time()[3]

for (R in 1:RR){
  if (R > check.FCR.grid){
    hyperPars.FCR.grid = c(seq(from=0, to=.9, len=10), .999)
    hyperPars.FCR.grid = cbind(hyperPars.FCR.grid, pUpper)
    FCR.grid.result = FCR(y.old=y.train[[R]], x.old=x.train[[R]], hyperPars=hyperPars.FCR.grid)
    betaSE.FCR.grid = rbind(betaSE.FCR.grid, apply(sweep(FCR.grid.result$beta, 2, beta.true)^2, 2, mean))
    if (R==RR) betaRMSE.FCR.grid = sqrt(apply(betaSE.FCR.grid, 2, mean))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[8] = proc.time()[3]
diff(ptm)

file = paste('Rimage/simul_',
             'scene=', scene, '_',
             switch(distri, 'normal_', 'exp_'),
             'seed=', seed, '_',
             'loop=', RR, '_',
             'ninitial=', ninitial, '_',
             'snr=', snr, '.RData',
             sep='')

save.image(file)

source("functions.R")
RMSE = integrate.rmse(betaRMSE.FCR.random
                      ,betaRMSE.FCR.grid
                      # ,betaRMSE.sFPCR
                      ,betaRMSE.pFPLS
                      ,betaRMSE.FPLSR
                      ,betaRMSE.FPCRR
                      ,betaRMSE.FPCR
                      )
source("functions.R")
curveplot.rmse(RMSE)
