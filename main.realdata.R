#' Put this r file into the folder containing functions.R
#' In the same folder, create two subfolders named 'figure' & 'Rimage', respectively.

if (!("rstudioapi" %in% rownames(installed.packages()))) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#### global setting
source("functions.R")
seed = 1
set.seed(seed)
RR = 500 # number of random splits
ninitial = 1 # number of initial values used in optimization
case.no = 4 # 1 for medfly data, others for tecator data (2 for meat vs fat, 3 for meat vs moisture, 3 for meat vs protein)

#### data

if (case.no == 1) { # medfly data
  load("data/medfly.rda")
  x = t(medfly$eggcount)
  y = as.matrix(medfly$lifetime)
}

if (case.no != 1) { # meat data
  meat = t(matrix(t(as.matrix(read.table("data/meat.txt"))), ncol=240))
  x = meat[, 1:100]
  
  if (case.no == 2) # meat vs fat
    y = as.matrix(meat[, 124])
  if (case.no == 3) # meat vs moisture
    y = as.matrix(meat[, 123])
  if (case.no == 4) # meat vs protein
    y = as.matrix(meat[, 125])
}  

# candidate pools for parameters to be tuned
N = nrow(x)
Theta.pFPLS = c(0, 1, 1e1, 1e2, 1e3, 1e4, 1e5)
Alpha.sFPCR = seq(from=0, to=1, len=11)
Theta.sFPCR = c(0, 1, 1e1, 1e2, 1e3, 1e4, 1e5)
nSubInt = 1e4L # number of subintervals in root-finding
pUpper = 5 # max number of components
pLower = 1 # min number of components

# Dangerous! clear existing ReMSPE
ReMSPE.FCR.random = NULL
ReMSPE.FCR.grid = NULL
ReMSPE.pFPLS = NULL
ReMSPE.FPCRR = NULL
ReMSPE.FPLSR = NULL
ReMSPE.sFPCR = NULL
ReMSPE.FPCR = NULL

# Check current progress
check.FCR.random = ifelse(is.null(ReMSPE.FCR.random), 0, length(ReMSPE.FCR.random))
check.FCR.grid = ifelse(is.null(ReMSPE.FCR.grid), 0, length(ReMSPE.FCR.grid))
check.pFPLS = ifelse(is.null(ReMSPE.pFPLS), 0, length(ReMSPE.pFPLS))
check.FPCRR = ifelse(is.null(ReMSPE.FPCRR), 0, length(ReMSPE.FPCRR))
check.FPLSR = ifelse(is.null(ReMSPE.FPLSR), 0, length(ReMSPE.FPLSR))
check.sFPCR = ifelse(is.null(ReMSPE.sFPCR), 0, length(ReMSPE.sFPCR))
check.FPCR = ifelse(is.null(ReMSPE.FPCR), 0, length(ReMSPE.FPCR))


################## Replica & Time ####################

index = list()
for (R in 1:RR){
  index[[R]] = sample(1:N, round(N * .9))
}

ptm = NULL
ptm[1] = proc.time()[3]

for (R in 1:RR){

  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]

  if (R > check.FCR.random){
    hyperPars.FCR.random = runif(10, min = 0, max = 1)
    hyperPars.FCR.random = cbind(hyperPars.FCR.random, sample.int(pUpper, size=length(hyperPars.FCR.random), replace=T))
    FCR.random.result = FCR(y.old, x.old, y.new, x.new, hyperPars=hyperPars.FCR.random)
    ReMSPE.FCR.random = c(ReMSPE.FCR.random, mean(FCR.random.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[2] = proc.time()[3]

for (R in 1:RR){
  
  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]
  
  if (R > check.sFPCR){
    sFPCR.result = sFPCR(y.old, x.old, y.new, x.new, Alpha.sFPCR, Theta.sFPCR, pUpper, pLower)
    ReMSPE.sFPCR = c(ReMSPE.sFPCR, mean(sFPCR.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[3] = proc.time()[3]

for (R in 1:RR){

  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]

  if (R > check.pFPLS){
    pFPLS.result = pFPLS(y.old, x.old, y.new, x.new, Theta.pFPLS, pUpper, pLower)
    ReMSPE.pFPLS = c(ReMSPE.pFPLS, mean(pFPLS.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[4] = proc.time()[3]

for (R in 1:RR){
  
  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]
  
  if (R > check.FPLSR){
    FPLSR.result = FPLS_R(y.old, x.old, y.new, x.new)
    ReMSPE.FPLSR = c(ReMSPE.FPLSR, mean(FPLSR.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[5] = proc.time()[3]

for (R in 1:RR){

  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]

  if (R > check.FPCRR){
    FPCRR.result = FPCR_R(y.old, x.old, y.new, x.new)
    ReMSPE.FPCRR = c(ReMSPE.FPCRR, mean(FPCRR.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[6] = proc.time()[3]

for (R in 1:RR){

  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]

  if (R > check.FPCR){
    FPCR.result = sFPCR(y.old, x.old, y.new, x.new, Alpha=1, Theta.sFPCR, pUpper, pLower)
    ReMSPE.FPCR = c(ReMSPE.FPCR, mean(FPCR.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[7] = proc.time()[3]

for (R in 1:RR){
  
  samp.idx = index[[R]]
  x.old = x[samp.idx, ]
  y.old = y[samp.idx]
  x.new = x[-samp.idx, ]
  y.new = y[-samp.idx]
  
  if (R > check.FCR.grid){
    hyperPars.FCR.grid = c(seq(from = 0, to = .9, length.out = 10), .999)
    hyperPars.FCR.grid = cbind(hyperPars.FCR.grid, pUpper)
    FCR.grid.result = FCR(y.old, x.old, y.new, x.new, hyperPars=hyperPars.FCR.grid)
    ReMSPE.FCR.grid = c(ReMSPE.FCR.grid, mean(FCR.grid.result$sse)/sum((y.new-mean(y.old))^2))
  }
  if (R %% 40 == 0) cat(R,'\n')
  else cat(R)
}

ptm[8] = proc.time()[3]
diff(ptm) # computing time consumed 

file = paste('Rimage/real_',
             switch(case.no, 'medfly_', 'fat_', 'moisture_', 'protein_'), 
             'seed=', seed, '_',
             'loop=', RR, '_',
             'ninitial=', ninitial, '.RData', 
             sep='')
save.image(file)

#' boxplots

ReMSPE = integrate.remspe(ReMSPE.FCR.random 
                          ,ReMSPE.FCR.grid
                          ,ReMSPE.sFPCR
                          ,ReMSPE.pFPLS 
                          ,ReMSPE.FPLSR
                          ,ReMSPE.FPCRR
                          ,ReMSPE.FPCR
                          )
boxplot.remspe(nrow(ReMSPE), ReMSPE)
