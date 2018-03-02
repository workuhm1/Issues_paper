

## ------------------------------------
## Bias, Type-I & Power Analysis of SEM
## ------------------------------------

all.files <- dir("./MIMIC/gamma/")
# all.files <- dir("./corr_ThetaCFA_ThetaMIMIC/")

files_out <- all.files[grep(".txt", all.files)]
bias_out <- data.frame(seed=integer(), a=integer(), b=integer(), c=numeric(), d=integer(), 
                       paramHeader=character(), param=character(),
                       est=numeric(), pval=numeric(), gammaTRUE=numeric())

# FIres <- data.frame(seed=integer(), a=integer(), b=integer(), c=numeric(), d=integer(),  
#                  rho1=numeric(), rho2=numeric(), 
#                  err1=character(), err2=character(), flipped=character())

# lambdaRes <- data.frame(seed=integer(), a=integer(), b=integer(), c=numeric(), d=integer(),  
#                      param=character(), est=numeric(), low2.5=numeric(), up2.5=numeric(), 
#                      True=numeric())

# corrCFAMIMIC <- data.frame(seed=integer(), 
#                            a=integer(), b=integer(), c=numeric(), d=integer(),  
#                            rho1_CFAMIMIC=numeric(), rho2_CFAMIMIC=numeric())


for(file in files_out){
  
  file2 <- paste("./MIMIC/gamma/", file, sep="")
  # file2 <- paste("./corr_ThetaCFA_ThetaMIMIC/", file, sep="")
  
  tmp <- read.table(file2, header=TRUE)
  bias_out <- rbind(bias_out, tmp)
  # lambdaRes <- rbind(lambdaRes, tmp)
  # FIres <- rbind(FIres, tmp)
  # corrCFAMIMIC <- rbind(corrCFAMIMIC, tmp)
}

 
## export the final data
write.table(bias_out, "./bias_out_interval.txt", sep=";", row.names=FALSE, col.names=TRUE)
# write.table(lambdaRes, "./MIMIC_lambda_interval.txt", sep=";", row.names=FALSE, col.names=TRUE)

## read the data
# bias_out_bin = read.table("bias_out_binary.txt", header=TRUE, sep=";")
# bias_out_interv = read.table("bias_out_interval.txt", header=TRUE, sep=";")

## merge the data
# bias_out = merge(bias_out_bin, bias_out_interv, all=TRUE)
# library(car)
# bias_out$e = Recode(bias_out$e, "NA=3; 1=1; 2=2 ")
# ## write
# write.table(bias_out, "./bias_out_combined.txt", sep=";", row.names=FALSE, col.names=TRUE)

bias_out = read.table("bias_out_combined.txt", header=TRUE, sep=";")

## ----------------
## Type-I Analysis 
## ----------------

## extract the beta's
betaGEN1 = bias_out[(bias_out$paramHeader=="F1.ON")&(bias_out$param=="GEN"), ]
betaGEN2 = bias_out[(bias_out$paramHeader=="F2.ON")&(bias_out$param=="GEN"), ]

## divide the data by communality
betaGEN1_commLow_N50 = betaGEN1[(betaGEN1$a==50)&(betaGEN1$d==1), ]
betaGEN1_commLow_N100 = betaGEN1[(betaGEN1$a==100)&(betaGEN1$d==1), ]
betaGEN1_commLow_N300 = betaGEN1[(betaGEN1$a==300)&(betaGEN1$d==1), ]
betaGEN1_commLow_N3000 = betaGEN1[(betaGEN1$a==3000)&(betaGEN1$d==1), ]

betaGEN1_commWide_N50 = betaGEN1[(betaGEN1$a==50)&(betaGEN1$d==2), ]
betaGEN1_commWide_N100 = betaGEN1[(betaGEN1$a==100)&(betaGEN1$d==2), ]
betaGEN1_commWide_N300 = betaGEN1[(betaGEN1$a==300)&(betaGEN1$d==2), ]
betaGEN1_commWide_N3000 = betaGEN1[(betaGEN1$a==3000)&(betaGEN1$d==2), ]

betaGEN1_commHigh_N50 = betaGEN1[(betaGEN1$a==50)&(betaGEN1$d==3), ]
betaGEN1_commHigh_N100 = betaGEN1[(betaGEN1$a==100)&(betaGEN1$d==3), ]
betaGEN1_commHigh_N300 = betaGEN1[(betaGEN1$a==300)&(betaGEN1$d==3), ]
betaGEN1_commHigh_N3000 = betaGEN1[(betaGEN1$a==3000)&(betaGEN1$d==3), ]

## --- GENDER - 2 ---
## divide the data by communality
betaGEN2_commLow_N50 = betaGEN2[(betaGEN2$a==50)&(betaGEN2$d==1), ]
betaGEN2_commLow_N100 = betaGEN2[(betaGEN2$a==100)&(betaGEN2$d==1), ]
betaGEN2_commLow_N300 = betaGEN2[(betaGEN2$a==300)&(betaGEN2$d==1), ]
betaGEN2_commLow_N3000 = betaGEN2[(betaGEN2$a==3000)&(betaGEN2$d==1), ]

betaGEN2_commWide_N50 = betaGEN2[(betaGEN2$a==50)&(betaGEN2$d==2), ]
betaGEN2_commWide_N100 = betaGEN2[(betaGEN2$a==100)&(betaGEN2$d==2), ]
betaGEN2_commWide_N300 = betaGEN2[(betaGEN2$a==300)&(betaGEN2$d==2), ]
betaGEN2_commWide_N3000 = betaGEN2[(betaGEN2$a==3000)&(betaGEN2$d==2), ]

betaGEN2_commHigh_N50 = betaGEN2[(betaGEN2$a==50)&(betaGEN2$d==3), ]
betaGEN2_commHigh_N100 = betaGEN2[(betaGEN2$a==100)&(betaGEN2$d==3), ]
betaGEN2_commHigh_N300 = betaGEN2[(betaGEN2$a==300)&(betaGEN2$d==3), ]
betaGEN2_commHigh_N3000 = betaGEN2[(betaGEN2$a==3000)&(betaGEN2$d==3), ]

## --- Age ---
## extract the beta's
betaAGE1 = bias_out[(bias_out$paramHeader=="F1.ON")&(bias_out$param=="AGE"), ]
betaAGE2 = bias_out[(bias_out$paramHeader=="F2.ON")&(bias_out$param=="AGE"), ]

## divide the data by communality
betaAGE1_commLow_N50 = betaAGE1[(betaAGE1$a==50)&(betaAGE1$d==1), ]
betaAGE1_commLow_N100 = betaAGE1[(betaAGE1$a==100)&(betaAGE1$d==1), ]
betaAGE1_commLow_N300 = betaAGE1[(betaAGE1$a==300)&(betaAGE1$d==1), ]
betaAGE1_commLow_N3000 = betaAGE1[(betaAGE1$a==3000)&(betaAGE1$d==1), ]

betaAGE1_commWide_N50 = betaAGE1[(betaAGE1$a==50)&(betaAGE1$d==2), ]
betaAGE1_commWide_N100 = betaAGE1[(betaAGE1$a==100)&(betaAGE1$d==2), ]
betaAGE1_commWide_N300 = betaAGE1[(betaAGE1$a==300)&(betaAGE1$d==2), ]
betaAGE1_commWide_N3000 = betaAGE1[(betaAGE1$a==3000)&(betaAGE1$d==2), ]

betaAGE1_commHigh_N50 = betaAGE1[(betaAGE1$a==50)&(betaAGE1$d==3), ]
betaAGE1_commHigh_N100 = betaAGE1[(betaAGE1$a==100)&(betaAGE1$d==3), ]
betaAGE1_commHigh_N300 = betaAGE1[(betaAGE1$a==300)&(betaAGE1$d==3), ]
betaAGE1_commHigh_N3000 = betaAGE1[(betaAGE1$a==3000)&(betaAGE1$d==3), ]


## --- factor 2 vs age
betaAGE2_commLow_N50 = betaAGE2[(betaAGE2$a==50)&(betaAGE2$d==1), ]
betaAGE2_commLow_N100 = betaAGE2[(betaAGE2$a==100)&(betaAGE2$d==1), ]
betaAGE2_commLow_N300 = betaAGE2[(betaAGE2$a==300)&(betaAGE2$d==1), ]
betaAGE2_commLow_N3000 = betaAGE2[(betaAGE2$a==3000)&(betaAGE2$d==1), ]

betaAGE2_commWide_N50 = betaAGE2[(betaAGE2$a==50)&(betaAGE2$d==2), ]
betaAGE2_commWide_N100 = betaAGE2[(betaAGE2$a==100)&(betaAGE2$d==2), ]
betaAGE2_commWide_N300 = betaAGE2[(betaAGE2$a==300)&(betaAGE2$d==2), ]
betaAGE2_commWide_N3000 = betaAGE2[(betaAGE2$a==3000)&(betaAGE2$d==2), ]

betaAGE2_commHigh_N50 = betaAGE2[(betaAGE2$a==50)&(betaAGE2$d==3), ]
betaAGE2_commHigh_N100 = betaAGE2[(betaAGE2$a==100)&(betaAGE2$d==3), ]
betaAGE2_commHigh_N300 = betaAGE2[(betaAGE2$a==300)&(betaAGE2$d==3), ]
betaAGE2_commHigh_N3000 = betaAGE2[(betaAGE2$a==3000)&(betaAGE2$d==3), ]



## --- calculate type-1 error proportion

## 1. Number of indicators
library(binom)
agresti_coull_CI = function(dat){
  
  tab = table(dat)
  if(is.na(tab[2])){
    x = 0
    n = tab[1] + x
  }else{
    x = tab[2]
    n = tab[1] + tab[2]
  }
  
  ac_CI = binom.confint(x=x, n=n, 0.95, methods="ac")
  return(ac_CI)
}

## 1a. communality=LOW; N=50
typeI_CommLowN50_P6 = agresti_coull_CI(betaGEN1_commLow_N50[betaGEN1_commLow_N50$b==6, 9] < 0.05)
typeI_CommLowN50_P10 = agresti_coull_CI(betaGEN1_commLow_N50[betaGEN1_commLow_N50$b==10, 9] < 0.05)
typeI_CommLowN50_P16 = agresti_coull_CI(betaGEN1_commLow_N50[betaGEN1_commLow_N50$b==16, 9] < 0.05)

## 1b. communality=LOW; N=100
typeI_CommLowN100_P6 = agresti_coull_CI(betaGEN1_commLow_N100[betaGEN1_commLow_N100$b==6, 9] < 0.05)
typeI_CommLowN100_P10 = agresti_coull_CI(betaGEN1_commLow_N100[betaGEN1_commLow_N100$b==10, 9] < 0.05)
typeI_CommLowN100_P16 = agresti_coull_CI(betaGEN1_commLow_N100[betaGEN1_commLow_N100$b==16, 9] < 0.05)

## 1c. communality=LOW; N=300
typeI_CommLowN300_P6 = agresti_coull_CI(betaGEN1_commLow_N300[betaGEN1_commLow_N300$b==6, 9] < 0.05)
typeI_CommLowN300_P10 = agresti_coull_CI(betaGEN1_commLow_N300[betaGEN1_commLow_N300$b==10, 9] < 0.05)
typeI_CommLowN300_P16 = agresti_coull_CI(betaGEN1_commLow_N300[betaGEN1_commLow_N300$b==16, 9] < 0.05)

## 1d. communality=LOW; N=3000
typeI_CommLowN3000_P6 = agresti_coull_CI(betaGEN1_commLow_N3000[betaGEN1_commLow_N3000$b==6, 9] < 0.05)
typeI_CommLowN3000_P10 = agresti_coull_CI(betaGEN1_commLow_N3000[betaGEN1_commLow_N3000$b==10, 9] < 0.05)
typeI_CommLowN3000_P16 = agresti_coull_CI(betaGEN1_commLow_N3000[betaGEN1_commLow_N3000$b==16, 9] < 0.05)

## 1e. communality=Wide; N=50
typeI_CommWideN50_P6 = agresti_coull_CI(betaGEN1_commWide_N50[betaGEN1_commWide_N50$b==6, 9] < 0.05)
typeI_CommWideN50_P10 = agresti_coull_CI(betaGEN1_commWide_N50[betaGEN1_commWide_N50$b==10, 9] < 0.05)
typeI_CommWideN50_P16 = agresti_coull_CI(betaGEN1_commWide_N50[betaGEN1_commWide_N50$b==16, 9] < 0.05)

## 1f. communality=Wide; N=100
typeI_CommWideN100_P6 = agresti_coull_CI(betaGEN1_commWide_N100[betaGEN1_commWide_N100$b==6, 9] < 0.05)
typeI_CommWideN100_P10 = agresti_coull_CI(betaGEN1_commWide_N100[betaGEN1_commWide_N100$b==10, 9] < 0.05)
typeI_CommWideN100_P16 = agresti_coull_CI(betaGEN1_commWide_N100[betaGEN1_commWide_N100$b==16, 9] < 0.05)

## 1g. communality=Wide; N=300
typeI_CommWideN300_P6 = agresti_coull_CI(betaGEN1_commWide_N300[betaGEN1_commWide_N300$b==6, 9] < 0.05)
typeI_CommWideN300_P10 = agresti_coull_CI(betaGEN1_commWide_N300[betaGEN1_commWide_N300$b==10, 9] < 0.05)
typeI_CommWideN300_P16 = agresti_coull_CI(betaGEN1_commWide_N300[betaGEN1_commWide_N300$b==16, 9] < 0.05)

## 1h. communality=Wide; N=3000
typeI_CommWideN3000_P6 = agresti_coull_CI(betaGEN1_commWide_N3000[betaGEN1_commWide_N3000$b==6, 9] < 0.05)
typeI_CommWideN3000_P10 = agresti_coull_CI(betaGEN1_commWide_N3000[betaGEN1_commWide_N3000$b==10, 9] < 0.05)
typeI_CommWideN3000_P16 = agresti_coull_CI(betaGEN1_commWide_N3000[betaGEN1_commWide_N3000$b==16, 9] < 0.05)

## 1h. communality=High; N=50
typeI_CommHighN50_P6 = agresti_coull_CI(betaGEN1_commHigh_N50[betaGEN1_commHigh_N50$b==6, 9] < 0.05)
typeI_CommHighN50_P10 = agresti_coull_CI(betaGEN1_commHigh_N50[betaGEN1_commHigh_N50$b==10, 9] < 0.05)
typeI_CommHighN50_P16 = agresti_coull_CI(betaGEN1_commHigh_N50[betaGEN1_commHigh_N50$b==16, 9] < 0.05)

## 1i. communality=High; N=100
typeI_CommHighN100_P6 = agresti_coull_CI(betaGEN1_commHigh_N100[betaGEN1_commHigh_N100$b==6, 9] < 0.05)
typeI_CommHighN100_P10 = agresti_coull_CI(betaGEN1_commHigh_N100[betaGEN1_commHigh_N100$b==10, 9] < 0.05)
typeI_CommHighN100_P16 = agresti_coull_CI(betaGEN1_commHigh_N100[betaGEN1_commHigh_N100$b==16, 9] < 0.05)

## 1j. communality=High; N=300
typeI_CommHighN300_P6 = agresti_coull_CI(betaGEN1_commHigh_N300[betaGEN1_commHigh_N300$b==6, 9] < 0.05)
typeI_CommHighN300_P10 = agresti_coull_CI(betaGEN1_commHigh_N300[betaGEN1_commHigh_N300$b==10, 9] < 0.05)
typeI_CommHighN300_P16 = agresti_coull_CI(betaGEN1_commHigh_N300[betaGEN1_commHigh_N300$b==16, 9] < 0.05)

## 1k. communality=High; N=3000
typeI_CommHighN3000_P6 = agresti_coull_CI(betaGEN1_commHigh_N3000[betaGEN1_commHigh_N3000$b==6, 9] < 0.05)
typeI_CommHighN3000_P10 = agresti_coull_CI(betaGEN1_commHigh_N3000[betaGEN1_commHigh_N3000$b==10, 9] < 0.05)
typeI_CommHighN3000_P16 = agresti_coull_CI(betaGEN1_commHigh_N3000[betaGEN1_commHigh_N3000$b==16, 9] < 0.05)

## --- GENDER - 2 ----

## 1.1. Number of indicators

## 1.1a. communality=LOW; N=50
typeI_CommLowN50_P6 = agresti_coull_CI(betaGEN2_commLow_N50[betaGEN2_commLow_N50$b==6, 9] < 0.05)
typeI_CommLowN50_P10 = agresti_coull_CI(betaGEN2_commLow_N50[betaGEN2_commLow_N50$b==10, 9] < 0.05)
typeI_CommLowN50_P16 = agresti_coull_CI(betaGEN2_commLow_N50[betaGEN2_commLow_N50$b==16, 9] < 0.05)

## 1.1b. communality=LOW; N=100
typeI_CommLowN100_P6 = agresti_coull_CI(betaGEN2_commLow_N100[betaGEN2_commLow_N100$b==6, 9] < 0.05)
typeI_CommLowN100_P10 = agresti_coull_CI(betaGEN2_commLow_N100[betaGEN2_commLow_N100$b==10, 9] < 0.05)
typeI_CommLowN100_P16 = agresti_coull_CI(betaGEN2_commLow_N100[betaGEN2_commLow_N100$b==16, 9] < 0.05)

## 1.1c. communality=LOW; N=300
typeI_CommLowN300_P6 = agresti_coull_CI(betaGEN2_commLow_N300[betaGEN2_commLow_N300$b==6, 9] < 0.05)
typeI_CommLowN300_P10 = agresti_coull_CI(betaGEN2_commLow_N300[betaGEN2_commLow_N300$b==10, 9] < 0.05)
typeI_CommLowN300_P16 = agresti_coull_CI(betaGEN2_commLow_N300[betaGEN2_commLow_N300$b==16, 9] < 0.05)

## 1.1d. communality=LOW; N=3000
typeI_CommLowN3000_P6 = agresti_coull_CI(betaGEN2_commLow_N3000[betaGEN2_commLow_N3000$b==6, 9] < 0.05)
typeI_CommLowN3000_P10 = agresti_coull_CI(betaGEN2_commLow_N3000[betaGEN2_commLow_N3000$b==10, 9] < 0.05)
typeI_CommLowN3000_P16 = agresti_coull_CI(betaGEN2_commLow_N3000[betaGEN2_commLow_N3000$b==16, 9] < 0.05)

## 1.1e. communality=Wide; N=50
typeI_CommWideN50_P6 = agresti_coull_CI(betaGEN2_commWide_N50[betaGEN2_commWide_N50$b==6, 9] < 0.05)
typeI_CommWideN50_P10 = agresti_coull_CI(betaGEN2_commWide_N50[betaGEN2_commWide_N50$b==10, 9] < 0.05)
typeI_CommWideN50_P16 = agresti_coull_CI(betaGEN2_commWide_N50[betaGEN2_commWide_N50$b==16, 9] < 0.05)

## 1.1f. communality=Wide; N=100
typeI_CommWideN100_P6 = agresti_coull_CI(betaGEN2_commWide_N100[betaGEN2_commWide_N100$b==6, 9] < 0.05)
typeI_CommWideN100_P10 = agresti_coull_CI(betaGEN2_commWide_N100[betaGEN2_commWide_N100$b==10, 9] < 0.05)
typeI_CommWideN100_P16 = agresti_coull_CI(betaGEN2_commWide_N100[betaGEN2_commWide_N100$b==16, 9] < 0.05)

## 1.1g. communality=Wide; N=300
typeI_CommWideN300_P6 = agresti_coull_CI(betaGEN2_commWide_N300[betaGEN2_commWide_N300$b==6, 9] < 0.05)
typeI_CommWideN300_P10 = agresti_coull_CI(betaGEN2_commWide_N300[betaGEN2_commWide_N300$b==10, 9] < 0.05)
typeI_CommWideN300_P16 = agresti_coull_CI(betaGEN2_commWide_N300[betaGEN2_commWide_N300$b==16, 9] < 0.05)

## 1.1h. communality=Wide; N=3000
typeI_CommWideN3000_P6 = agresti_coull_CI(betaGEN2_commWide_N3000[betaGEN2_commWide_N3000$b==6, 9] < 0.05)
typeI_CommWideN3000_P10 = agresti_coull_CI(betaGEN2_commWide_N3000[betaGEN2_commWide_N3000$b==10, 9] < 0.05)
typeI_CommWideN3000_P16 = agresti_coull_CI(betaGEN2_commWide_N3000[betaGEN2_commWide_N3000$b==16, 9] < 0.05)

## 1.1h. communality=High; N=50
typeI_CommHighN50_P6 = agresti_coull_CI(betaGEN2_commHigh_N50[betaGEN2_commHigh_N50$b==6, 9] < 0.05)
typeI_CommHighN50_P10 = agresti_coull_CI(betaGEN2_commHigh_N50[betaGEN2_commHigh_N50$b==10, 9] < 0.05)
typeI_CommHighN50_P16 = agresti_coull_CI(betaGEN2_commHigh_N50[betaGEN2_commHigh_N50$b==16, 9] < 0.05)

## 1.1i. communality=High; N=100
typeI_CommHighN100_P6 = agresti_coull_CI(betaGEN2_commHigh_N100[betaGEN2_commHigh_N100$b==6, 9] < 0.05)
typeI_CommHighN100_P10 = agresti_coull_CI(betaGEN2_commHigh_N100[betaGEN2_commHigh_N100$b==10, 9] < 0.05)
typeI_CommHighN100_P16 = agresti_coull_CI(betaGEN2_commHigh_N100[betaGEN2_commHigh_N100$b==16, 9] < 0.05)

## 1.1j. communality=High; N=300
typeI_CommHighN300_P6 = agresti_coull_CI(betaGEN2_commHigh_N300[betaGEN2_commHigh_N300$b==6, 9] < 0.05)
typeI_CommHighN300_P10 = agresti_coull_CI(betaGEN2_commHigh_N300[betaGEN2_commHigh_N300$b==10, 9] < 0.05)
typeI_CommHighN300_P16 = agresti_coull_CI(betaGEN2_commHigh_N300[betaGEN2_commHigh_N300$b==16, 9] < 0.05)

## 1.1k. communality=High; N=3000
typeI_CommHighN3000_P6 = agresti_coull_CI(betaGEN2_commHigh_N3000[betaGEN2_commHigh_N3000$b==6, 9] < 0.05)
typeI_CommHighN3000_P10 = agresti_coull_CI(betaGEN2_commHigh_N3000[betaGEN2_commHigh_N3000$b==10, 9] < 0.05)
typeI_CommHighN3000_P16 = agresti_coull_CI(betaGEN2_commHigh_N3000[betaGEN2_commHigh_N3000$b==16, 9] < 0.05)


## ------------
## AGE: Type-I 
## ------------

## 2a. communality=LOW; N=50
typeI_CommLowN50_P6 = agresti_coull_CI(betaAGE1_commLow_N50[betaAGE1_commLow_N50$b==6, 9] < 0.05)
typeI_CommLowN50_P10 = agresti_coull_CI(betaAGE1_commLow_N50[betaAGE1_commLow_N50$b==10, 9] < 0.05)
typeI_CommLowN50_P16 = agresti_coull_CI(betaAGE1_commLow_N50[betaAGE1_commLow_N50$b==16, 9] < 0.05)

## 2b. communality=LOW; N=100
typeI_CommLowN100_P6 = agresti_coull_CI(betaAGE1_commLow_N100[betaAGE1_commLow_N100$b==6, 9] < 0.05)
typeI_CommLowN100_P10 = agresti_coull_CI(betaAGE1_commLow_N100[betaAGE1_commLow_N100$b==10, 9] < 0.05)
typeI_CommLowN100_P16 = agresti_coull_CI(betaAGE1_commLow_N100[betaAGE1_commLow_N100$b==16, 9] < 0.05)

## 2c. communality=LOW; N=300
typeI_CommLowN300_P6 = agresti_coull_CI(betaAGE1_commLow_N300[betaAGE1_commLow_N300$b==6, 9] < 0.05)
typeI_CommLowN300_P10 = agresti_coull_CI(betaAGE1_commLow_N300[betaAGE1_commLow_N300$b==10, 9] < 0.05)
typeI_CommLowN300_P16 = agresti_coull_CI(betaAGE1_commLow_N300[betaAGE1_commLow_N300$b==16, 9] < 0.05)

## 2d. communality=LOW; N=3000
typeI_CommLowN3000_P6 = agresti_coull_CI(betaAGE1_commLow_N3000[betaAGE1_commLow_N3000$b==6, 9] < 0.05)
typeI_CommLowN3000_P10 = agresti_coull_CI(betaAGE1_commLow_N3000[betaAGE1_commLow_N3000$b==10, 9] < 0.05)
typeI_CommLowN3000_P16 = agresti_coull_CI(betaAGE1_commLow_N3000[betaAGE1_commLow_N3000$b==16, 9] < 0.05)

## 2e. communality=WIDE; N=50
typeI_CommWideN50_P6 = agresti_coull_CI(betaAGE1_commWide_N50[betaAGE1_commWide_N50$b==6, 9] < 0.05)
typeI_CommWideN50_P10 = agresti_coull_CI(betaAGE1_commWide_N50[betaAGE1_commWide_N50$b==10, 9] < 0.05)
typeI_CommWideN50_P16 = agresti_coull_CI(betaAGE1_commWide_N50[betaAGE1_commWide_N50$b==16, 9] < 0.05)

## 2f. communality=WIDE; N=100
typeI_CommWideN100_P6 = agresti_coull_CI(betaAGE1_commWide_N100[betaAGE1_commWide_N100$b==6, 9] < 0.05)
typeI_CommWideN100_P10 = agresti_coull_CI(betaAGE1_commWide_N100[betaAGE1_commWide_N100$b==10, 9] < 0.05)
typeI_CommWideN100_P16 = agresti_coull_CI(betaAGE1_commWide_N100[betaAGE1_commWide_N100$b==16, 9] < 0.05)

## 2g. communality=WIDE; N=300
typeI_CommWideN300_P6 = agresti_coull_CI(betaAGE1_commWide_N300[betaAGE1_commWide_N300$b==6, 9] < 0.05)
typeI_CommWideN300_P10 = agresti_coull_CI(betaAGE1_commWide_N300[betaAGE1_commWide_N300$b==10, 9] < 0.05)
typeI_CommWideN300_P16 = agresti_coull_CI(betaAGE1_commWide_N300[betaAGE1_commWide_N300$b==16, 9] < 0.05)

## 2h. communality=WIDE; N=3000
typeI_CommWideN3000_P6 = agresti_coull_CI(betaAGE1_commWide_N3000[betaAGE1_commWide_N3000$b==6, 9] < 0.05)
typeI_CommWideN3000_P10 = agresti_coull_CI(betaAGE1_commWide_N3000[betaAGE1_commWide_N3000$b==10, 9] < 0.05)
typeI_CommWideN3000_P16 = agresti_coull_CI(betaAGE1_commWide_N3000[betaAGE1_commWide_N3000$b==16, 9] < 0.05)

## 2i. communality=HIGH; N=50
typeI_CommHighN50_P6 = agresti_coull_CI(betaAGE1_commHigh_N50[betaAGE1_commHigh_N50$b==6, 9] < 0.05)
typeI_CommHighN50_P10 = agresti_coull_CI(betaAGE1_commHigh_N50[betaAGE1_commHigh_N50$b==10, 9] < 0.05)
typeI_CommHighN50_P16 = agresti_coull_CI(betaAGE1_commHigh_N50[betaAGE1_commHigh_N50$b==16, 9] < 0.05)

## 2i. communality=HIGH; N=100
typeI_CommHighN100_P6 = agresti_coull_CI(betaAGE1_commHigh_N100[betaAGE1_commHigh_N100$b==6, 9] < 0.05)
typeI_CommHighN100_P10 = agresti_coull_CI(betaAGE1_commHigh_N100[betaAGE1_commHigh_N100$b==10, 9] < 0.05)
typeI_CommHighN100_P16 = agresti_coull_CI(betaAGE1_commHigh_N100[betaAGE1_commHigh_N100$b==16, 9] < 0.05)

## 2j. communality=HIGH; N=300
typeI_CommHighN300_P6 = agresti_coull_CI(betaAGE1_commHigh_N300[betaAGE1_commHigh_N300$b==6, 9] < 0.05)
typeI_CommHighN300_P10 = agresti_coull_CI(betaAGE1_commHigh_N300[betaAGE1_commHigh_N300$b==10, 9] < 0.05)
typeI_CommHighN300_P16 = agresti_coull_CI(betaAGE1_commHigh_N300[betaAGE1_commHigh_N300$b==16, 9] < 0.05)

## 2k. communality=HIGH; N=3000
typeI_CommHighN3000_P6 = agresti_coull_CI(betaAGE1_commHigh_N3000[betaAGE1_commHigh_N3000$b==6, 9] < 0.05)
typeI_CommHighN3000_P10 = agresti_coull_CI(betaAGE1_commHigh_N3000[betaAGE1_commHigh_N3000$b==10, 9] < 0.05)
typeI_CommHighN3000_P16 = agresti_coull_CI(betaAGE1_commHigh_N3000[betaAGE1_commHigh_N3000$b==16, 9] < 0.05)


## ---------------
## Age Vs factor 2
## ---------------

## 2.1a. communality=LOW; N=50
typeI_CommLowN50_P6 = agresti_coull_CI(betaAGE2_commLow_N50[betaAGE2_commLow_N50$b==6, 9] < 0.05)
typeI_CommLowN50_P10 = agresti_coull_CI(betaAGE2_commLow_N50[betaAGE2_commLow_N50$b==10, 9] < 0.05)
typeI_CommLowN50_P16 = agresti_coull_CI(betaAGE2_commLow_N50[betaAGE2_commLow_N50$b==16, 9] < 0.05)

## 2.1b. communality=LOW; N=100
typeI_CommLowN100_P6 = agresti_coull_CI(betaAGE2_commLow_N100[betaAGE2_commLow_N100$b==6, 9] < 0.05)
typeI_CommLowN100_P10 = agresti_coull_CI(betaAGE2_commLow_N100[betaAGE2_commLow_N100$b==10, 9] < 0.05)
typeI_CommLowN100_P16 = agresti_coull_CI(betaAGE2_commLow_N100[betaAGE2_commLow_N100$b==16, 9] < 0.05)

## 2.1c. communality=LOW; N=300
typeI_CommLowN300_P6 = agresti_coull_CI(betaAGE2_commLow_N300[betaAGE2_commLow_N300$b==6, 9] < 0.05)
typeI_CommLowN300_P10 = agresti_coull_CI(betaAGE2_commLow_N300[betaAGE2_commLow_N300$b==10, 9] < 0.05)
typeI_CommLowN300_P16 = agresti_coull_CI(betaAGE2_commLow_N300[betaAGE2_commLow_N300$b==16, 9] < 0.05)

## 2.1d. communality=LOW; N=3000
typeI_CommLowN3000_P6 = agresti_coull_CI(betaAGE2_commLow_N3000[betaAGE2_commLow_N3000$b==6, 9] < 0.05)
typeI_CommLowN3000_P10 = agresti_coull_CI(betaAGE2_commLow_N3000[betaAGE2_commLow_N3000$b==10, 9] < 0.05)
typeI_CommLowN3000_P16 = agresti_coull_CI(betaAGE2_commLow_N3000[betaAGE2_commLow_N3000$b==16, 9] < 0.05)

## 2.1e. communality=WIDE; N=50
typeI_CommWideN50_P6 = agresti_coull_CI(betaAGE2_commWide_N50[betaAGE2_commWide_N50$b==6, 9] < 0.05)
typeI_CommWideN50_P10 = agresti_coull_CI(betaAGE2_commWide_N50[betaAGE2_commWide_N50$b==10, 9] < 0.05)
typeI_CommWideN50_P16 = agresti_coull_CI(betaAGE2_commWide_N50[betaAGE2_commWide_N50$b==16, 9] < 0.05)

## 2.1f. communality=WIDE; N=100
typeI_CommWideN100_P6 = agresti_coull_CI(betaAGE2_commWide_N100[betaAGE2_commWide_N100$b==6, 9] < 0.05)
typeI_CommWideN100_P10 = agresti_coull_CI(betaAGE2_commWide_N100[betaAGE2_commWide_N100$b==10, 9] < 0.05)
typeI_CommWideN100_P16 = agresti_coull_CI(betaAGE2_commWide_N100[betaAGE2_commWide_N100$b==16, 9] < 0.05)

## 2.1g. communality=WIDE; N=300
typeI_CommWideN300_P6 = agresti_coull_CI(betaAGE2_commWide_N300[betaAGE2_commWide_N300$b==6, 9] < 0.05)
typeI_CommWideN300_P10 = agresti_coull_CI(betaAGE2_commWide_N300[betaAGE2_commWide_N300$b==10, 9] < 0.05)
typeI_CommWideN300_P16 = agresti_coull_CI(betaAGE2_commWide_N300[betaAGE2_commWide_N300$b==16, 9] < 0.05)

## 2.1h. communality=WIDE; N=3000
typeI_CommWideN3000_P6 = agresti_coull_CI(betaAGE2_commWide_N3000[betaAGE2_commWide_N3000$b==6, 9] < 0.05)
typeI_CommWideN3000_P10 = agresti_coull_CI(betaAGE2_commWide_N3000[betaAGE2_commWide_N3000$b==10, 9] < 0.05)
typeI_CommWideN3000_P16 = agresti_coull_CI(betaAGE2_commWide_N3000[betaAGE2_commWide_N3000$b==16, 9] < 0.05)

## 2.1i. communality=HIGH; N=50
typeI_CommHighN50_P6 = agresti_coull_CI(betaAGE2_commHigh_N50[betaAGE2_commHigh_N50$b==6, 9] < 0.05)
typeI_CommHighN50_P10 = agresti_coull_CI(betaAGE2_commHigh_N50[betaAGE2_commHigh_N50$b==10, 9] < 0.05)
typeI_CommHighN50_P16 = agresti_coull_CI(betaAGE2_commHigh_N50[betaAGE2_commHigh_N50$b==16, 9] < 0.05)

## 2.1i. communality=HIGH; N=100
typeI_CommHighN100_P6 = agresti_coull_CI(betaAGE2_commHigh_N100[betaAGE2_commHigh_N100$b==6, 9] < 0.05)
typeI_CommHighN100_P10 = agresti_coull_CI(betaAGE2_commHigh_N100[betaAGE2_commHigh_N100$b==10, 9] < 0.05)
typeI_CommHighN100_P16 = agresti_coull_CI(betaAGE2_commHigh_N100[betaAGE2_commHigh_N100$b==16, 9] < 0.05)

## 2.1j. communality=HIGH; N=300
typeI_CommHighN300_P6 = agresti_coull_CI(betaAGE2_commHigh_N300[betaAGE2_commHigh_N300$b==6, 9] < 0.05)
typeI_CommHighN300_P10 = agresti_coull_CI(betaAGE2_commHigh_N300[betaAGE2_commHigh_N300$b==10, 9] < 0.05)
typeI_CommHighN300_P16 = agresti_coull_CI(betaAGE2_commHigh_N300[betaAGE2_commHigh_N300$b==16, 9] < 0.05)

## 2.1.k. communality=HIGH; N=3000
typeI_CommHighN3000_P6 = agresti_coull_CI(betaAGE2_commHigh_N3000[betaAGE2_commHigh_N3000$b==6, 9] < 0.05)
typeI_CommHighN3000_P10 = agresti_coull_CI(betaAGE2_commHigh_N3000[betaAGE2_commHigh_N3000$b==10, 9] < 0.05)
typeI_CommHighN3000_P16 = agresti_coull_CI(betaAGE2_commHigh_N3000[betaAGE2_commHigh_N3000$b==16, 9] < 0.05)


## -------------------------------
## Power Analysis, consider Neurot
## -------------------------------

## extract the beta's
betaN1 = bias_out[(bias_out$paramHeader=="F1.ON")&(bias_out$param=="N"), ]
betaN2 = bias_out[(bias_out$paramHeader=="F2.ON")&(bias_out$param=="N"), ]

## divide the data by communality
## first factor
betaN1_commLow_N50 = betaN1[(betaN1$a==50)&(betaN1$d==1), ]
betaN1_commLow_N100 = betaN1[(betaN1$a==100)&(betaN1$d==1), ]
betaN1_commLow_N300 = betaN1[(betaN1$a==300)&(betaN1$d==1), ]
betaN1_commLow_N3000 = betaN1[(betaN1$a==3000)&(betaN1$d==1), ]

betaN1_commWide_N50 = betaN1[(betaN1$a==50)&(betaN1$d==2), ]
betaN1_commWide_N100 = betaN1[(betaN1$a==100)&(betaN1$d==2), ]
betaN1_commWide_N300 = betaN1[(betaN1$a==300)&(betaN1$d==2), ]
betaN1_commWide_N3000 = betaN1[(betaN1$a==3000)&(betaN1$d==2), ]

betaN1_commHigh_N50 = betaN1[(betaN1$a==50)&(betaN1$d==3), ]
betaN1_commHigh_N100 = betaN1[(betaN1$a==100)&(betaN1$d==3), ]
betaN1_commHigh_N300 = betaN1[(betaN1$a==300)&(betaN1$d==3), ]
betaN1_commHigh_N3000 = betaN1[(betaN1$a==3000)&(betaN1$d==3), ]


## second factor
betaN2_commLow_N50 = betaN2[(betaN2$a==50)&(betaN2$d==1), ]
betaN2_commLow_N100 = betaN2[(betaN2$a==100)&(betaN2$d==1), ]
betaN2_commLow_N300 = betaN2[(betaN2$a==300)&(betaN2$d==1), ]
betaN2_commLow_N3000 = betaN2[(betaN2$a==3000)&(betaN2$d==1), ]

betaN2_commWide_N50 = betaN2[(betaN2$a==50)&(betaN2$d==2), ]
betaN2_commWide_N100 = betaN2[(betaN2$a==100)&(betaN2$d==2), ]
betaN2_commWide_N300 = betaN2[(betaN2$a==300)&(betaN2$d==2), ]
betaN2_commWide_N3000 = betaN2[(betaN2$a==3000)&(betaN2$d==2), ]

betaN2_commHigh_N50 = betaN2[(betaN2$a==50)&(betaN2$d==3), ]
betaN2_commHigh_N100 = betaN2[(betaN2$a==100)&(betaN2$d==3), ]
betaN2_commHigh_N300 = betaN2[(betaN2$a==300)&(betaN2$d==3), ]
betaN2_commHigh_N3000 = betaN2[(betaN2$a==3000)&(betaN2$d==3), ]


## --- calculate power ---

## 3. Number of indicators

## 3a. communality=LOW; N=50
power_CommLowN50_P6 = round(mean(betaN1_commLow_N50[betaN1_commLow_N50$b==6, 9] < 0.05), 3)
power_CommLowN50_P10 = round(mean(betaN1_commLow_N50[betaN1_commLow_N50$b==10, 9] < 0.05), 3)
power_CommLowN50_P16 = round(mean(betaN1_commLow_N50[betaN1_commLow_N50$b==16, 9] < 0.05), 3)

## 3b. communality=LOW; N=100
power_CommLowN100_P6 = round(mean(betaN1_commLow_N100[betaN1_commLow_N100$b==6, 9] < 0.05), 3)
power_CommLowN100_P10 = round(mean(betaN1_commLow_N100[betaN1_commLow_N100$b==10, 9] < 0.05), 3)
power_CommLowN100_P16 = round(mean(betaN1_commLow_N100[betaN1_commLow_N100$b==16, 9] < 0.05), 3)

## 3c. communality=LOW; N=300
power_CommLowN300_P6 = round(mean(betaN1_commLow_N300[betaN1_commLow_N300$b==6, 9] < 0.05), 3)
power_CommLowN300_P10 = round(mean(betaN1_commLow_N300[betaN1_commLow_N300$b==10, 9] < 0.05), 3)
power_CommLowN300_P16 = round(mean(betaN1_commLow_N300[betaN1_commLow_N300$b==16, 9] < 0.05), 3)

## 3d. communality=LOW; N=3000
power_CommLowN3000_P6 = round(mean(betaN1_commLow_N3000[betaN1_commLow_N3000$b==6, 9] < 0.05), 3)
power_CommLowN3000_P10 = round(mean(betaN1_commLow_N3000[betaN1_commLow_N3000$b==10, 9] < 0.05), 3)
power_CommLowN3000_P16 = round(mean(betaN1_commLow_N3000[betaN1_commLow_N3000$b==16, 9] < 0.05), 3)

## 3e. communality=WIDE; N=50
power_CommWideN50_P6 = round(mean(betaN1_commWide_N50[betaN1_commWide_N50$b==6, 9] < 0.05), 3)
power_CommWideN50_P10 = round(mean(betaN1_commWide_N50[betaN1_commWide_N50$b==10, 9] < 0.05), 3)
power_CommWideN50_P16 = round(mean(betaN1_commWide_N50[betaN1_commWide_N50$b==16, 9] < 0.05), 3)

## 3f. communality=WIDE; N=100
power_CommWideN100_P6 = round(mean(betaN1_commWide_N100[betaN1_commWide_N100$b==6, 9] < 0.05), 3)
power_CommWideN100_P10 = round(mean(betaN1_commWide_N100[betaN1_commWide_N100$b==10, 9] < 0.05), 3)
power_CommWideN100_P16 = round(mean(betaN1_commWide_N100[betaN1_commWide_N100$b==16, 9] < 0.05), 3)

## 3f. communality=WIDE; N=300
power_CommWideN300_P6 = round(mean(betaN1_commWide_N300[betaN1_commWide_N300$b==6, 9] < 0.05), 3)
power_CommWideN300_P10 = round(mean(betaN1_commWide_N300[betaN1_commWide_N300$b==10, 9] < 0.05), 3)
power_CommWideN300_P16 = round(mean(betaN1_commWide_N300[betaN1_commWide_N300$b==16, 9] < 0.05), 3)

## 3g. communality=WIDE; N=3000
power_CommWideN3000_P6 = round(mean(betaN1_commWide_N3000[betaN1_commWide_N3000$b==6, 9] < 0.05), 3)
power_CommWideN3000_P10 = round(mean(betaN1_commWide_N3000[betaN1_commWide_N3000$b==10, 9] < 0.05), 3)
power_CommWideN3000_P16 = round(mean(betaN1_commWide_N3000[betaN1_commWide_N3000$b==16, 9] < 0.05), 3)

## 3h. communality=HIGH; N=50
power_CommHighN50_P6 = round(mean(betaN1_commHigh_N50[betaN1_commHigh_N50$b==6, 9] < 0.05), 3)
power_CommHighN50_P10 = round(mean(betaN1_commHigh_N50[betaN1_commHigh_N50$b==10, 9] < 0.05), 3)
power_CommHighN50_P16 = round(mean(betaN1_commHigh_N50[betaN1_commHigh_N50$b==16, 9] < 0.05), 3)

## 3i. communality=HIGH; N=100
power_CommHighN100_P6 = round(mean(betaN1_commHigh_N100[betaN1_commHigh_N100$b==6, 9] < 0.05), 3)
power_CommHighN100_P10 = round(mean(betaN1_commHigh_N100[betaN1_commHigh_N100$b==10, 9] < 0.05), 3)
power_CommHighN100_P16 = round(mean(betaN1_commHigh_N100[betaN1_commHigh_N100$b==16, 9] < 0.05), 3)

## 3j. communality=HIGH; N=300
power_CommHighN300_P6 = round(mean(betaN1_commHigh_N300[betaN1_commHigh_N300$b==6, 9] < 0.05), 3)
power_CommHighN300_P10 = round(mean(betaN1_commHigh_N300[betaN1_commHigh_N300$b==10, 9] < 0.05), 3)
power_CommHighN300_P16 = round(mean(betaN1_commHigh_N300[betaN1_commHigh_N300$b==16, 9] < 0.05), 3)

## 3k. communality=HIGH; N=3000
power_CommHighN3000_P6 = round(mean(betaN1_commHigh_N3000[betaN1_commHigh_N3000$b==6, 9] < 0.05), 3)
power_CommHighN3000_P10 = round(mean(betaN1_commHigh_N3000[betaN1_commHigh_N3000$b==10, 9] < 0.05), 3)
power_CommHighN3000_P16 = round(mean(betaN1_commHigh_N3000[betaN1_commHigh_N3000$b==16, 9] < 0.05), 3)


## --- second factor vs neurot -----

## 3.1. Number of indicators

## 3.1a. communality=LOW; N=50
power_CommLowN50_P6 = round(mean(betaN2_commLow_N50[betaN2_commLow_N50$b==6, 9] < 0.05), 3)
power_CommLowN50_P10 = round(mean(betaN2_commLow_N50[betaN2_commLow_N50$b==10, 9] < 0.05), 3)
power_CommLowN50_P16 = round(mean(betaN2_commLow_N50[betaN2_commLow_N50$b==16, 9] < 0.05), 3)

## 3.1b. communality=LOW; N=100
power_CommLowN100_P6 = round(mean(betaN2_commLow_N100[betaN2_commLow_N100$b==6, 9] < 0.05), 3)
power_CommLowN100_P10 = round(mean(betaN2_commLow_N100[betaN2_commLow_N100$b==10, 9] < 0.05), 3)
power_CommLowN100_P16 = round(mean(betaN2_commLow_N100[betaN2_commLow_N100$b==16, 9] < 0.05), 3)

## 3.1c. communality=LOW; N=300
power_CommLowN300_P6 = round(mean(betaN2_commLow_N300[betaN2_commLow_N300$b==6, 9] < 0.05), 3)
power_CommLowN300_P10 = round(mean(betaN2_commLow_N300[betaN2_commLow_N300$b==10, 9] < 0.05), 3)
power_CommLowN300_P16 = round(mean(betaN2_commLow_N300[betaN2_commLow_N300$b==16, 9] < 0.05), 3)

## 3.1d. communality=LOW; N=3000
power_CommLowN3000_P6 = round(mean(betaN2_commLow_N3000[betaN2_commLow_N3000$b==6, 9] < 0.05), 3)
power_CommLowN3000_P10 = round(mean(betaN2_commLow_N3000[betaN2_commLow_N3000$b==10, 9] < 0.05), 3)
power_CommLowN3000_P16 = round(mean(betaN2_commLow_N3000[betaN2_commLow_N3000$b==16, 9] < 0.05), 3)

## 3.1e. communality=WIDE; N=50
power_CommWideN50_P6 = round(mean(betaN2_commWide_N50[betaN2_commWide_N50$b==6, 9] < 0.05), 3)
power_CommWideN50_P10 = round(mean(betaN2_commWide_N50[betaN2_commWide_N50$b==10, 9] < 0.05), 3)
power_CommWideN50_P16 = round(mean(betaN2_commWide_N50[betaN2_commWide_N50$b==16, 9] < 0.05), 3)

## 3.1f. communality=WIDE; N=100
power_CommWideN100_P6 = round(mean(betaN2_commWide_N100[betaN2_commWide_N100$b==6, 9] < 0.05), 3)
power_CommWideN100_P10 = round(mean(betaN2_commWide_N100[betaN2_commWide_N100$b==10, 9] < 0.05), 3)
power_CommWideN100_P16 = round(mean(betaN2_commWide_N100[betaN2_commWide_N100$b==16, 9] < 0.05), 3)

## 3.1f. communality=WIDE; N=300
power_CommWideN300_P6 = round(mean(betaN2_commWide_N300[betaN2_commWide_N300$b==6, 9] < 0.05), 3)
power_CommWideN300_P10 = round(mean(betaN2_commWide_N300[betaN2_commWide_N300$b==10, 9] < 0.05), 3)
power_CommWideN300_P16 = round(mean(betaN2_commWide_N300[betaN2_commWide_N300$b==16, 9] < 0.05), 3)

## 3.1g. communality=WIDE; N=3000
power_CommWideN3000_P6 = round(mean(betaN2_commWide_N3000[betaN2_commWide_N3000$b==6, 9] < 0.05), 3)
power_CommWideN3000_P10 = round(mean(betaN2_commWide_N3000[betaN2_commWide_N3000$b==10, 9] < 0.05), 3)
power_CommWideN3000_P16 = round(mean(betaN2_commWide_N3000[betaN2_commWide_N3000$b==16, 9] < 0.05), 3)

## 3.1h. communality=HIGH; N=50
power_CommHighN50_P6 = round(mean(betaN2_commHigh_N50[betaN2_commHigh_N50$b==6, 9] < 0.05), 3)
power_CommHighN50_P10 = round(mean(betaN2_commHigh_N50[betaN2_commHigh_N50$b==10, 9] < 0.05), 3)
power_CommHighN50_P16 = round(mean(betaN2_commHigh_N50[betaN2_commHigh_N50$b==16, 9] < 0.05), 3)

## 3.1i. communality=HIGH; N=100
power_CommHighN100_P6 = round(mean(betaN2_commHigh_N100[betaN2_commHigh_N100$b==6, 9] < 0.05), 3)
power_CommHighN100_P10 = round(mean(betaN2_commHigh_N100[betaN2_commHigh_N100$b==10, 9] < 0.05), 3)
power_CommHighN100_P16 = round(mean(betaN2_commHigh_N100[betaN2_commHigh_N100$b==16, 9] < 0.05), 3)

## 3.1j. communality=HIGH; N=300
power_CommHighN300_P6 = round(mean(betaN2_commHigh_N300[betaN2_commHigh_N300$b==6, 9] < 0.05), 3)
power_CommHighN300_P10 = round(mean(betaN2_commHigh_N300[betaN2_commHigh_N300$b==10, 9] < 0.05), 3)
power_CommHighN300_P16 = round(mean(betaN2_commHigh_N300[betaN2_commHigh_N300$b==16, 9] < 0.05), 3)

## 3.1k. communality=HIGH; N=3000
power_CommHighN3000_P6 = round(mean(betaN2_commHigh_N3000[betaN2_commHigh_N3000$b==6, 9] < 0.05), 3)
power_CommHighN3000_P10 = round(mean(betaN2_commHigh_N3000[betaN2_commHigh_N3000$b==10, 9] < 0.05), 3)
power_CommHighN3000_P16 = round(mean(betaN2_commHigh_N3000[betaN2_commHigh_N3000$b==16, 9] < 0.05), 3)


## --- Education ---
betaEDU1 = bias_out[(bias_out$paramHeader=="F1.ON")&(bias_out$param=="EDU"), ]
betaEDU2 = bias_out[(bias_out$paramHeader=="F2.ON")&(bias_out$param=="EDU"), ]

## divide the data by communality
## first factor
betaEDU1_commLow_N50 = betaEDU1[(betaEDU1$a==50)&(betaEDU1$d==1), ]
betaEDU1_commLow_N100 = betaEDU1[(betaEDU1$a==100)&(betaEDU1$d==1), ]
betaEDU1_commLow_N300 = betaEDU1[(betaEDU1$a==300)&(betaEDU1$d==1), ]
betaEDU1_commLow_N3000 = betaEDU1[(betaEDU1$a==3000)&(betaEDU1$d==1), ]

betaEDU1_commWide_N50 = betaEDU1[(betaEDU1$a==50)&(betaEDU1$d==2), ]
betaEDU1_commWide_N100 = betaEDU1[(betaEDU1$a==100)&(betaEDU1$d==2), ]
betaEDU1_commWide_N300 = betaEDU1[(betaEDU1$a==300)&(betaEDU1$d==2), ]
betaEDU1_commWide_N3000 = betaEDU1[(betaEDU1$a==3000)&(betaEDU1$d==2), ]

betaEDU1_commHigh_N50 = betaEDU1[(betaEDU1$a==50)&(betaEDU1$d==3), ]
betaEDU1_commHigh_N100 = betaEDU1[(betaEDU1$a==100)&(betaEDU1$d==3), ]
betaEDU1_commHigh_N300 = betaEDU1[(betaEDU1$a==300)&(betaEDU1$d==3), ]
betaEDU1_commHigh_N3000 = betaEDU1[(betaEDU1$a==3000)&(betaEDU1$d==3), ]


## second factor
betaEDU2_commLow_N50 = betaEDU2[(betaEDU2$a==50)&(betaEDU2$d==1), ]
betaEDU2_commLow_N100 = betaEDU2[(betaEDU2$a==100)&(betaEDU2$d==1), ]
betaEDU2_commLow_N300 = betaEDU2[(betaEDU2$a==300)&(betaEDU2$d==1), ]
betaEDU2_commLow_N3000 = betaEDU2[(betaEDU2$a==3000)&(betaEDU2$d==1), ]

betaEDU2_commWide_N50 = betaEDU2[(betaEDU2$a==50)&(betaEDU2$d==2), ]
betaEDU2_commWide_N100 = betaEDU2[(betaEDU2$a==100)&(betaEDU2$d==2), ]
betaEDU2_commWide_N300 = betaEDU2[(betaEDU2$a==300)&(betaEDU2$d==2), ]
betaEDU2_commWide_N3000 = betaEDU2[(betaEDU2$a==3000)&(betaEDU2$d==2), ]

betaEDU2_commHigh_N50 = betaEDU2[(betaEDU2$a==50)&(betaEDU2$d==3), ]
betaEDU2_commHigh_N100 = betaEDU2[(betaEDU2$a==100)&(betaEDU2$d==3), ]
betaEDU2_commHigh_N300 = betaEDU2[(betaEDU2$a==300)&(betaEDU2$d==3), ]
betaEDU2_commHigh_N3000 = betaEDU2[(betaEDU2$a==3000)&(betaEDU2$d==3), ]


## --- power analysis on EDU ---

## 4. Number of indicators

## 4a. communality=LOW; N=50
power_CommLowN50_P6 = round(mean(betaEDU1_commLow_N50[betaEDU1_commLow_N50$b==6, 9] < 0.05), 3)
power_CommLowN50_P10 = round(mean(betaEDU1_commLow_N50[betaEDU1_commLow_N50$b==10, 9] < 0.05), 3)
power_CommLowN50_P16 = round(mean(betaEDU1_commLow_N50[betaEDU1_commLow_N50$b==16, 9] < 0.05), 3)

## 4b. communality=LOW; N=100
power_CommLowN100_P6 = round(mean(betaEDU1_commLow_N100[betaEDU1_commLow_N100$b==6, 9] < 0.05), 3)
power_CommLowN100_P10 = round(mean(betaEDU1_commLow_N100[betaEDU1_commLow_N100$b==10, 9] < 0.05), 3)
power_CommLowN100_P16 = round(mean(betaEDU1_commLow_N100[betaEDU1_commLow_N100$b==16, 9] < 0.05), 3)

## 4c. communality=LOW; N=300
power_CommLowN300_P6 = round(mean(betaEDU1_commLow_N300[betaEDU1_commLow_N300$b==6, 9] < 0.05), 3)
power_CommLowN300_P10 = round(mean(betaEDU1_commLow_N300[betaEDU1_commLow_N300$b==10, 9] < 0.05), 3)
power_CommLowN300_P16 = round(mean(betaEDU1_commLow_N300[betaEDU1_commLow_N300$b==16, 9] < 0.05), 3)

## 4d. communality=LOW; N=3000
power_CommLowN3000_P6 = round(mean(betaEDU1_commLow_N3000[betaEDU1_commLow_N3000$b==6, 9] < 0.05), 3)
power_CommLowN3000_P10 = round(mean(betaEDU1_commLow_N3000[betaEDU1_commLow_N3000$b==10, 9] < 0.05), 3)
power_CommLowN3000_P16 = round(mean(betaEDU1_commLow_N3000[betaEDU1_commLow_N3000$b==16, 9] < 0.05), 3)

## 4e. communality=WIDE; N=50
power_CommWideN50_P6 = round(mean(betaEDU1_commWide_N50[betaEDU1_commWide_N50$b==6, 9] < 0.05), 3)
power_CommWideN50_P10 = round(mean(betaEDU1_commWide_N50[betaEDU1_commWide_N50$b==10, 9] < 0.05), 3)
power_CommWideN50_P16 = round(mean(betaEDU1_commWide_N50[betaEDU1_commWide_N50$b==16, 9] < 0.05), 3)

## 4f. communality=WIDE; N=100
power_CommWideN100_P6 = round(mean(betaEDU1_commWide_N100[betaEDU1_commWide_N100$b==6, 9] < 0.05), 3)
power_CommWideN100_P10 = round(mean(betaEDU1_commWide_N100[betaEDU1_commWide_N100$b==10, 9] < 0.05), 3)
power_CommWideN100_P16 = round(mean(betaEDU1_commWide_N100[betaEDU1_commWide_N100$b==16, 9] < 0.05), 3)

## 4g. communality=WIDE; N=300
power_CommWideN300_P6 = round(mean(betaEDU1_commWide_N300[betaEDU1_commWide_N300$b==6, 9] < 0.05), 3)
power_CommWideN300_P10 = round(mean(betaEDU1_commWide_N300[betaEDU1_commWide_N300$b==10, 9] < 0.05), 3)
power_CommWideN300_P16 = round(mean(betaEDU1_commWide_N300[betaEDU1_commWide_N300$b==16, 9] < 0.05), 3)

## 4h. communality=WIDE; N=3000
power_CommWideN3000_P6 = round(mean(betaEDU1_commWide_N3000[betaEDU1_commWide_N3000$b==6, 9] < 0.05), 3)
power_CommWideN3000_P10 = round(mean(betaEDU1_commWide_N3000[betaEDU1_commWide_N3000$b==10, 9] < 0.05), 3)
power_CommWideN3000_P16 = round(mean(betaEDU1_commWide_N3000[betaEDU1_commWide_N3000$b==16, 9] < 0.05), 3)

## 4i. communality=HIGH; N=50
power_CommHighN50_P6 = round(mean(betaEDU1_commHigh_N50[betaEDU1_commHigh_N50$b==6, 9] < 0.05), 3)
power_CommHighN50_P10 = round(mean(betaEDU1_commHigh_N50[betaEDU1_commHigh_N50$b==10, 9] < 0.05), 3)
power_CommHighN50_P16 = round(mean(betaEDU1_commHigh_N50[betaEDU1_commHigh_N50$b==16, 9] < 0.05), 3)

## 4j. communality=HIGH; N=100
power_CommHighN100_P6 = round(mean(betaEDU1_commHigh_N100[betaEDU1_commHigh_N100$b==6, 9] < 0.05), 3)
power_CommHighN100_P10 = round(mean(betaEDU1_commHigh_N100[betaEDU1_commHigh_N100$b==10, 9] < 0.05), 3)
power_CommHighN100_P16 = round(mean(betaEDU1_commHigh_N100[betaEDU1_commHigh_N100$b==16, 9] < 0.05), 3)

## 4k. communality=HIGH; N=300
power_CommHighN300_P6 = round(mean(betaEDU1_commHigh_N300[betaEDU1_commHigh_N300$b==6, 9] < 0.05), 3)
power_CommHighN300_P10 = round(mean(betaEDU1_commHigh_N300[betaEDU1_commHigh_N300$b==10, 9] < 0.05), 3)
power_CommHighN300_P16 = round(mean(betaEDU1_commHigh_N300[betaEDU1_commHigh_N300$b==16, 9] < 0.05), 3)

## 4l. communality=HIGH; N=3000
power_CommHighN3000_P6 = round(mean(betaEDU1_commHigh_N3000[betaEDU1_commHigh_N3000$b==6, 9] < 0.05), 3)
power_CommHighN3000_P10 = round(mean(betaEDU1_commHigh_N3000[betaEDU1_commHigh_N3000$b==10, 9] < 0.05), 3)
power_CommHighN3000_P16 = round(mean(betaEDU1_commHigh_N3000[betaEDU1_commHigh_N3000$b==16, 9] < 0.05), 3)

## --- factor 2 ---

## 4.1.a. communality=LOW; N=50
power_CommLowN50_P6 = round(mean(betaEDU2_commLow_N50[betaEDU2_commLow_N50$b==6, 9] < 0.05), 3)
power_CommLowN50_P10 = round(mean(betaEDU2_commLow_N50[betaEDU2_commLow_N50$b==10, 9] < 0.05), 3)
power_CommLowN50_P16 = round(mean(betaEDU2_commLow_N50[betaEDU2_commLow_N50$b==16, 9] < 0.05), 3)

## 4.1.b. communality=LOW; N=100
power_CommLowN100_P6 = round(mean(betaEDU2_commLow_N100[betaEDU2_commLow_N100$b==6, 9] < 0.05), 3)
power_CommLowN100_P10 = round(mean(betaEDU2_commLow_N100[betaEDU2_commLow_N100$b==10, 9] < 0.05), 3)
power_CommLowN100_P16 = round(mean(betaEDU2_commLow_N100[betaEDU2_commLow_N100$b==16, 9] < 0.05), 3)

## 4.1.c. communality=LOW; N=300
power_CommLowN300_P6 = round(mean(betaEDU2_commLow_N300[betaEDU2_commLow_N300$b==6, 9] < 0.05), 3)
power_CommLowN300_P10 = round(mean(betaEDU2_commLow_N300[betaEDU2_commLow_N300$b==10, 9] < 0.05), 3)
power_CommLowN300_P16 = round(mean(betaEDU2_commLow_N300[betaEDU2_commLow_N300$b==16, 9] < 0.05), 3)

## 4.1.d. communality=LOW; N=3000
power_CommLowN3000_P6 = round(mean(betaEDU2_commLow_N3000[betaEDU2_commLow_N3000$b==6, 9] < 0.05), 3)
power_CommLowN3000_P10 = round(mean(betaEDU2_commLow_N3000[betaEDU2_commLow_N3000$b==10, 9] < 0.05), 3)
power_CommLowN3000_P16 = round(mean(betaEDU2_commLow_N3000[betaEDU2_commLow_N3000$b==16, 9] < 0.05), 3)

## 4.1.e. communality=WIDE; N=50
power_CommWideN50_P6 = round(mean(betaEDU2_commWide_N50[betaEDU2_commWide_N50$b==6, 9] < 0.05), 3)
power_CommWideN50_P10 = round(mean(betaEDU2_commWide_N50[betaEDU2_commWide_N50$b==10, 9] < 0.05), 3)
power_CommWideN50_P16 = round(mean(betaEDU2_commWide_N50[betaEDU2_commWide_N50$b==16, 9] < 0.05), 3)

## 4.1.f. communality=WIDE; N=100
power_CommWideN100_P6 = round(mean(betaEDU2_commWide_N100[betaEDU2_commWide_N100$b==6, 9] < 0.05), 3)
power_CommWideN100_P10 = round(mean(betaEDU2_commWide_N100[betaEDU2_commWide_N100$b==10, 9] < 0.05), 3)
power_CommWideN100_P16 = round(mean(betaEDU2_commWide_N100[betaEDU2_commWide_N100$b==16, 9] < 0.05), 3)

## 4.1.g. communality=WIDE; N=300
power_CommWideN300_P6 = round(mean(betaEDU2_commWide_N300[betaEDU2_commWide_N300$b==6, 9] < 0.05), 3)
power_CommWideN300_P10 = round(mean(betaEDU2_commWide_N300[betaEDU2_commWide_N300$b==10, 9] < 0.05), 3)
power_CommWideN300_P16 = round(mean(betaEDU2_commWide_N300[betaEDU2_commWide_N300$b==16, 9] < 0.05), 3)

## 4.1.h. communality=WIDE; N=3000
power_CommWideN3000_P6 = round(mean(betaEDU2_commWide_N3000[betaEDU2_commWide_N3000$b==6, 9] < 0.05), 3)
power_CommWideN3000_P10 = round(mean(betaEDU2_commWide_N3000[betaEDU2_commWide_N3000$b==10, 9] < 0.05), 3)
power_CommWideN3000_P16 = round(mean(betaEDU2_commWide_N3000[betaEDU2_commWide_N3000$b==16, 9] < 0.05), 3)

## 4.1.i. communality=HIGH; N=50
power_CommHighN50_P6 = round(mean(betaEDU2_commHigh_N50[betaEDU2_commHigh_N50$b==6, 9] < 0.05), 3)
power_CommHighN50_P10 = round(mean(betaEDU2_commHigh_N50[betaEDU2_commHigh_N50$b==10, 9] < 0.05), 3)
power_CommHighN50_P16 = round(mean(betaEDU2_commHigh_N50[betaEDU2_commHigh_N50$b==16, 9] < 0.05), 3)

## 4.1.j. communality=HIGH; N=100
power_CommHighN100_P6 = round(mean(betaEDU2_commHigh_N100[betaEDU2_commHigh_N100$b==6, 9] < 0.05), 3)
power_CommHighN100_P10 = round(mean(betaEDU2_commHigh_N100[betaEDU2_commHigh_N100$b==10, 9] < 0.05), 3)
power_CommHighN100_P16 = round(mean(betaEDU2_commHigh_N100[betaEDU2_commHigh_N100$b==16, 9] < 0.05), 3)

## 4.1.k. communality=HIGH; N=300
power_CommHighN300_P6 = round(mean(betaEDU2_commHigh_N300[betaEDU2_commHigh_N300$b==6, 9] < 0.05), 3)
power_CommHighN300_P10 = round(mean(betaEDU2_commHigh_N300[betaEDU2_commHigh_N300$b==10, 9] < 0.05), 3)
power_CommHighN300_P16 = round(mean(betaEDU2_commHigh_N300[betaEDU2_commHigh_N300$b==16, 9] < 0.05), 3)

## 4.1.l. communality=HIGH; N=3000
power_CommHighN3000_P6 = round(mean(betaEDU2_commHigh_N3000[betaEDU2_commHigh_N3000$b==6, 9] < 0.05), 3)
power_CommHighN3000_P10 = round(mean(betaEDU2_commHigh_N3000[betaEDU2_commHigh_N3000$b==10, 9] < 0.05), 3)
power_CommHighN3000_P16 = round(mean(betaEDU2_commHigh_N3000[betaEDU2_commHigh_N3000$b==16, 9] < 0.05), 3)



