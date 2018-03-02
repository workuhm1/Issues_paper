
source("dgenerate_binary.r")
source("dualModel.r")
source("dualCFA.r")
source("dualMIMIC.r")

## ---------------------
## STEP-1: Simulate data
## ---------------------

## sample size
A = c(50, 100,300,3000)

## no. variables per factor (B)
B = c(6,10,16)

## correlation (C)
C = c(0.0,0.4,0.8)

## communality (D)
D_lambda = list(c(0.316,0.447), c(0.316,0.632), c(0.632,0.775))   

## marginal probability (pi)
E = list(c(0.05,0.15), c(0.40,0.50))

## prototype
# all_seed <- round(runif(50,1000,9999))   
all_seed <- round(runif(50,10000,19999))   

## run the data generating function, rep=100
for(rep in all_seed){
  ## set the data generating seed
  
  set.seed(rep)  
  
  dgenerateBinary(N=A, P=B, rho_lat=C, lambda=D_lambda, prob=E, seed=rep)
}


## -------------------------
## using Parallel computing 
## -------------------------

# library(doMC)
library(doParallel)
library(foreach)

# Find out how many cores are available (if you don't already know)
detectCores()

# Create cluster with desired number of cores
# registerDoMC(2)
cl <- makeCluster(2)
registerDoParallel(cl)

## START
# all_files <- dir()
all_files <- dir("./dataBinary/")
# files_RData <- all_files[grep(".RData", dir())]
files_RData <- all_files[grep(".RData", all_files)]

# start time
strt <- Sys.time()

# ptime <- system.time({
# r <- foreach(file = files_RData, .combine = cbind)%dopar% {
# ls <- foreach(file = files_RData) %dopar% {
ls <- foreach(file = files_RData) %do% {
  # for(file in files_RData){
  # r <- foreach(icount(file = files_RData))%dopar% {
  #   cat(file)
  
  file_out <- gsub(pattern=".RData", replacement=".txt", x=file)
  #     load(file)
  
  file2 <- paste("./dataBinary/", file, sep="")
  load(file2)
  
  dualModel(dataY=YBin, dataXX = XX, file=file_out, seed=seed)
  
}
# ptime
print(Sys.time() - strt)
stopCluster(cl)

