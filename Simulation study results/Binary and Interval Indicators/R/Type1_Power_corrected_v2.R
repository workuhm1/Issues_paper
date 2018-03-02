

## -------------------------------------------------------------  ##
## Mark suggested to modify the Type-I error C.I.s based on the   ##
##    5-way interaction data, i.e., 324 possible combinations     ##
##      [04 SEP 2015, Leiden]                                     ##
## -------------------------------------------------------------  ##


# import combined data
bias_out = read.table("bias_out_combined.txt", header=TRUE, sep=";")


## ----------------
## Type-I Analysis 
## ----------------

## extract the gamma's
gammaGEN1 = bias_out[(bias_out$paramHeader=="F1.ON")&(bias_out$param=="GEN"), ]
gammaGEN2 = bias_out[(bias_out$paramHeader=="F2.ON")&(bias_out$param=="GEN"), ]

## --- Age ---
## extract the gamma's
gammaAGE1 = bias_out[(bias_out$paramHeader=="F1.ON")&(bias_out$param=="AGE"), ]
gammaAGE2 = bias_out[(bias_out$paramHeader=="F2.ON")&(bias_out$param=="AGE"), ]

## power: gamma=0.95
powerGammaBig <- bias_out[bias_out$gammaTRUE == 0.95, ]

## power: gamma=-0.30
powerGammaModerate <- bias_out[bias_out$gammaTRUE == -0.3, ]

## power: gamma=0.10
powerGammaSmall <- bias_out[bias_out$gammaTRUE == 0.1, ]





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


## ------
## ----------------------------------------------------------------- ##
## Predictor: Age corresponding to the first factor, i.e., gammaAGE1 ##
## ----------------------------------------------------------------- ##
## ------

## --- BLR

## Table-1 [ROW-1]: ToIs=BLR; Commun=LOW; and, Corr=0.0
row1_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-2]: ToIs=BLR; Commun=LOW; and, Corr=0.4
row2_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-3]: ToIs=BLR; Commun=LOW; and, Corr=0.8
row3_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0.8), ]

## Table-1 [ROW-4]: ToIs=BLR; Commun=WIDE; and, Corr=0
row4_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-5]: ToIs=BLR; Commun=WIDE; and, Corr=0.4
row5_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-6]: ToIs=BLR; Commun=WIDE; and, Corr=0.8
row6_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0.8), ]

## Table-1 [ROW-7]: ToIs=BLR; Commun=HIGH; and, Corr=0.0
row7_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-8]: ToIs=BLR; Commun=HIGH; and, Corr=0.4
row8_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-9]: ToIs=BLR; Commun=HIGH; and, Corr=0.8
row9_AGE1 <- gammaAGE1[(gammaAGE1$e == 1)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0.8), ]


## --- BMR

## Table-1 [ROW-10]: ToIs=BMR; Commun=LOW; and, Corr=0.0
row10_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-11]: ToIs=BMR; Commun=LOW; and, Corr=0.4
row11_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-12]: ToIs=BMR; Commun=LOW; and, Corr=0.8
row12_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0.8), ]

## Table-1 [ROW-13]: ToIs=BMR; Commun=WIDE; and, Corr=0.0
row13_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-14]: ToIs=BMR; Commun=WIDE; and, Corr=0.4
row14_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-15]: ToIs=BMR; Commun=WIDE; and, Corr=0.8
row15_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0.8), ]

## Table-1 [ROW-16]: ToIs=BMR; Commun=HIGH; and, Corr=0.0
row16_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-17]: ToIs=BMR; Commun=HIGH; and, Corr=0.4
row17_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-18]: ToIs=BMR; Commun=HIGH; and, Corr=0.8
row18_AGE1 <- gammaAGE1[(gammaAGE1$e == 2)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0.8), ]


## --- Interval

## Table-1 [ROW-19]: ToIs=Interval; Commun=LOW; and, Corr=0.0
row19_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-20]: ToIs=Interval; Commun=LOW; and, Corr=0.4
row20_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-21]: ToIs=Interval; Commun=LOW; and, Corr=0.8
row21_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 1)&(gammaAGE1$c == 0.8), ]

## Table-1 [ROW-22]: ToIs=Interval; Commun=WIDE; and, Corr=0.0
row22_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-23]: ToIs=Interval; Commun=WIDE; and, Corr=0.4
row23_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-24]: ToIs=Interval; Commun=WIDE; and, Corr=0.8
row24_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 2)&(gammaAGE1$c == 0.8), ]

## Table-1 [ROW-25]: ToIs=Interval; Commun=HIGH; and, Corr=0
row25_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0), ]

## Table-1 [ROW-26]: ToIs=Interval; Commun=HIGH; and, Corr=0.4
row26_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0.4), ]

## Table-1 [ROW-27]: ToIs=Interval; Commun=HIGH; and, Corr=0.8
row27_AGE1 <- gammaAGE1[(gammaAGE1$e == 3)&(gammaAGE1$d == 3)&(gammaAGE1$c == 0.8), ]




## -----------------------------------  ##
## calculate the C.Is for Type-I error  ##
## -----------------------------------  ##


##--- ROW-1

## AGE1-ROW-2: N=50, no estimates because those models were NOT converged!!!
##    table(row1_AGE1$a): N=100(freq=88); 300(289); 3000 (300)

## AGE1-ROW-1: N=100
typeI_N100_P6 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 100)&(row1_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 100)&(row1_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 100)&(row1_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-1:N=300
typeI_N300_P6 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 300)&(row1_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 300)&(row1_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 300)&(row1_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-1:N=3000
typeI_N3000_P6 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 3000)&(row1_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 3000)&(row1_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row1_AGE1[((row1_AGE1$a == 3000)&(row1_AGE1$b == 16)), 10] < 0.05)


##--- ROW-2

## AGE1-ROW-2: N=50, no estimates because those models were NOT converged!!!
##    table(row1_AGE1$a): N=100(freq=87); 300(287); 3000(300) 

## AGE1-ROW-2:N=100
typeI_N100_P6 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 100)&(row2_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 100)&(row2_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 100)&(row2_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-2:N=300
typeI_N300_P6 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 300)&(row2_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 300)&(row2_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 300)&(row2_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-2:N=3000
typeI_N3000_P6 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 3000)&(row2_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 3000)&(row2_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row2_AGE1[((row2_AGE1$a == 3000)&(row2_AGE1$b == 16)), 10] < 0.05)


##--- ROW-3

## AGE1-ROW-3: N=50, no estimates because those models were NOT converged!!!
##    table(row1_AGE1$a): N=100(freq=94); 300(285); 3000(300) 

## AGE1-ROW-3:N=100
typeI_N100_P6 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 100)&(row3_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 100)&(row3_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 100)&(row3_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-3:N=300
typeI_N300_P6 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 300)&(row3_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 300)&(row3_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 300)&(row3_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-3:N=3000
typeI_N3000_P6 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 3000)&(row3_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 3000)&(row3_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row3_AGE1[((row3_AGE1$a == 3000)&(row3_AGE1$b == 16)), 10] < 0.05)


##--- ROW-4

## AGE1-ROW-4: N=50, no estimates because those models were NOT converged!!!
##    table(row4_AGE1$a): N=100(freq=32); 300(223); 3000 (299)

## AGE1-ROW-4: N=100
typeI_N100_P6 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 100)&(row4_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 100)&(row4_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 100)&(row4_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-4: N=300
typeI_N300_P6 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 300)&(row4_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 300)&(row4_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 300)&(row4_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-4: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 3000)&(row4_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 3000)&(row4_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row4_AGE1[((row4_AGE1$a == 3000)&(row4_AGE1$b == 16)), 10] < 0.05)


##--- ROW-5

## AGE1-ROW-5: N=50, no estimates because those models were NOT converged!!!
##    table(row5_AGE1$a): N=100(freq=45); 300(226); 3000 (299)

## AGE1-ROW-5: N=100
typeI_N100_P6 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 100)&(row5_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 100)&(row5_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 100)&(row5_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-5: N=300
typeI_N300_P6 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 300)&(row5_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 300)&(row5_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 300)&(row5_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-5: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 3000)&(row5_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 3000)&(row5_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row5_AGE1[((row5_AGE1$a == 3000)&(row5_AGE1$b == 16)), 10] < 0.05)


##--- ROW-6

## AGE1-ROW-6: N=50
typeI_N50_P6 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 50)&(row6_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 50)&(row6_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 50)&(row6_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-6: N=100
typeI_N100_P6 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 100)&(row6_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 100)&(row6_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 100)&(row6_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-6: N=300
typeI_N300_P6 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 300)&(row6_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 300)&(row6_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 300)&(row6_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-6: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 3000)&(row6_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 3000)&(row6_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row6_AGE1[((row6_AGE1$a == 3000)&(row6_AGE1$b == 16)), 10] < 0.05)


##--- ROW-7

## AGE1-ROW-7: N=50
typeI_N50_P6 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 50)&(row7_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 50)&(row7_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 50)&(row7_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-7: N=100
typeI_N100_P6 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 100)&(row7_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 100)&(row7_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 100)&(row7_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-7: N=300
typeI_N300_P6 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 300)&(row7_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 300)&(row7_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 300)&(row7_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-7: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 3000)&(row7_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 3000)&(row7_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row7_AGE1[((row7_AGE1$a == 3000)&(row7_AGE1$b == 16)), 10] < 0.05)


##--- ROW-8

## AGE1-ROW-8: N=50
typeI_N50_P6 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 50)&(row8_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 50)&(row8_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 50)&(row8_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-8: N=100
typeI_N100_P6 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 100)&(row8_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 100)&(row8_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 100)&(row8_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-8: N=300
typeI_N300_P6 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 300)&(row8_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 300)&(row8_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 300)&(row8_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-8: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 3000)&(row8_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 3000)&(row8_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row8_AGE1[((row8_AGE1$a == 3000)&(row8_AGE1$b == 16)), 10] < 0.05)


##--- ROW-9

## AGE1-ROW-9: N=50
typeI_N50_P6 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 50)&(row9_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 50)&(row9_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 50)&(row9_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-9: N=100
typeI_N100_P6 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 100)&(row9_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 100)&(row9_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 100)&(row9_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-9: N=300
typeI_N300_P6 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 300)&(row9_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 300)&(row9_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 300)&(row9_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-9: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 3000)&(row9_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 3000)&(row9_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row9_AGE1[((row9_AGE1$a == 3000)&(row9_AGE1$b == 16)), 10] < 0.05)


## -----------------  
## --- calculate BMR
## -----------------

##--- ROW-10

## AGE1-ROW-10: N=50
typeI_N50_P6 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 50)&(row10_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 50)&(row10_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 50)&(row10_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-10: N=100
typeI_N100_P6 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 100)&(row10_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 100)&(row10_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 100)&(row10_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-10: N=300
typeI_N300_P6 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 300)&(row10_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 300)&(row10_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 300)&(row10_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-10: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 3000)&(row10_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 3000)&(row10_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row10_AGE1[((row10_AGE1$a == 3000)&(row10_AGE1$b == 16)), 10] < 0.05)


##--- ROW-11

## AGE1-ROW-11: N=50
typeI_N50_P6 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 50)&(row11_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 50)&(row11_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 50)&(row11_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-11: N=100
typeI_N100_P6 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 100)&(row11_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 100)&(row11_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 100)&(row11_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-11: N=300
typeI_N300_P6 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 300)&(row11_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 300)&(row11_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 300)&(row11_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-11: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 3000)&(row11_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 3000)&(row11_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row11_AGE1[((row11_AGE1$a == 3000)&(row11_AGE1$b == 16)), 10] < 0.05)


##--- ROW-12

## AGE1-ROW-12: N=50
typeI_N50_P6 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 50)&(row12_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 50)&(row12_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 50)&(row12_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-12: N=100
typeI_N100_P6 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 100)&(row12_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 100)&(row12_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 100)&(row12_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-12: N=300
typeI_N300_P6 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 300)&(row12_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 300)&(row12_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 300)&(row12_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-12: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 3000)&(row12_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 3000)&(row12_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row12_AGE1[((row12_AGE1$a == 3000)&(row12_AGE1$b == 16)), 10] < 0.05)


##--- ROW-13

## AGE1-ROW-13: N=50
typeI_N50_P6 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 50)&(row13_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 50)&(row13_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 50)&(row13_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-13: N=100
typeI_N100_P6 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 100)&(row13_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 100)&(row13_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 100)&(row13_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-13: N=300
typeI_N300_P6 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 300)&(row13_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 300)&(row13_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 300)&(row13_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-13: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 3000)&(row13_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 3000)&(row13_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row13_AGE1[((row13_AGE1$a == 3000)&(row13_AGE1$b == 16)), 10] < 0.05)


##--- ROW-14

## AGE1-ROW-14: N=50
typeI_N50_P6 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 50)&(row14_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 50)&(row14_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 50)&(row14_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-14: N=100
typeI_N100_P6 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 100)&(row14_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 100)&(row14_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 100)&(row14_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-14: N=300
typeI_N300_P6 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 300)&(row14_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 300)&(row14_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 300)&(row14_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-14: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 3000)&(row14_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 3000)&(row14_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row14_AGE1[((row14_AGE1$a == 3000)&(row14_AGE1$b == 16)), 10] < 0.05)


##--- ROW-15

## AGE1-ROW-15: N=50
typeI_N50_P6 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 50)&(row15_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 50)&(row15_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 50)&(row15_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-15: N=100
typeI_N100_P6 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 100)&(row15_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 100)&(row15_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 100)&(row15_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-15: N=300
typeI_N300_P6 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 300)&(row15_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 300)&(row15_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 300)&(row15_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-15: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 3000)&(row15_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 3000)&(row15_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row15_AGE1[((row15_AGE1$a == 3000)&(row15_AGE1$b == 16)), 10] < 0.05)


##--- ROW-16

## AGE1-ROW-16: N=50
typeI_N50_P6 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 50)&(row16_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 50)&(row16_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 50)&(row16_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-16: N=100
typeI_N100_P6 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 100)&(row16_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 100)&(row16_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 100)&(row16_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-16: N=300
typeI_N300_P6 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 300)&(row16_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 300)&(row16_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 300)&(row16_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-16: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 3000)&(row16_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 3000)&(row16_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row16_AGE1[((row16_AGE1$a == 3000)&(row16_AGE1$b == 16)), 10] < 0.05)


##--- ROW-17

## AGE1-ROW-17: N=50
typeI_N50_P6 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 50)&(row17_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 50)&(row17_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 50)&(row17_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-17: N=100
typeI_N100_P6 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 100)&(row17_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 100)&(row17_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 100)&(row17_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-17: N=300
typeI_N300_P6 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 300)&(row17_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 300)&(row17_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 300)&(row17_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-17: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 3000)&(row17_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 3000)&(row17_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row17_AGE1[((row17_AGE1$a == 3000)&(row17_AGE1$b == 16)), 10] < 0.05)


##--- ROW-18

## AGE1-ROW-18: N=50
typeI_N50_P6 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 50)&(row18_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 50)&(row18_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 50)&(row18_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-18: N=100
typeI_N100_P6 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 100)&(row18_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 100)&(row18_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 100)&(row18_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-18: N=300
typeI_N300_P6 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 300)&(row18_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 300)&(row18_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 300)&(row18_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-18: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 3000)&(row18_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 3000)&(row18_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row18_AGE1[((row18_AGE1$a == 3000)&(row18_AGE1$b == 16)), 10] < 0.05)



## ----------------------
## --- calculate Interval
## ----------------------

##--- ROW-19

## AGE1-ROW-19: N=50
typeI_N50_P6 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 50)&(row19_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 50)&(row19_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 50)&(row19_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-19: N=100
typeI_N100_P6 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 100)&(row19_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 100)&(row19_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 100)&(row19_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-19: N=300
typeI_N300_P6 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 300)&(row19_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 300)&(row19_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 300)&(row19_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-19: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 3000)&(row19_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 3000)&(row19_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row19_AGE1[((row19_AGE1$a == 3000)&(row19_AGE1$b == 16)), 10] < 0.05)


##--- ROW-20

## AGE1-ROW-20: N=50
typeI_N50_P6 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 50)&(row20_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 50)&(row20_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 50)&(row20_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-20: N=100
typeI_N100_P6 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 100)&(row20_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 100)&(row20_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 100)&(row20_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-20: N=300
typeI_N300_P6 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 300)&(row20_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 300)&(row20_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 300)&(row20_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-20: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 3000)&(row20_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 3000)&(row20_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row20_AGE1[((row20_AGE1$a == 3000)&(row20_AGE1$b == 16)), 10] < 0.05)


##--- ROW-21

## AGE1-ROW-21: N=50
typeI_N50_P6 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 50)&(row21_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 50)&(row21_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 50)&(row21_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-21: N=100
typeI_N100_P6 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 100)&(row21_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 100)&(row21_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 100)&(row21_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-21: N=300
typeI_N300_P6 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 300)&(row21_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 300)&(row21_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 300)&(row21_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-21: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 3000)&(row21_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 3000)&(row21_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row21_AGE1[((row21_AGE1$a == 3000)&(row21_AGE1$b == 16)), 10] < 0.05)


##--- ROW-22

## AGE1-ROW-22: N=50
typeI_N50_P6 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 50)&(row22_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 50)&(row22_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 50)&(row22_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-22: N=100
typeI_N100_P6 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 100)&(row22_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 100)&(row22_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 100)&(row22_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-22: N=300
typeI_N300_P6 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 300)&(row22_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 300)&(row22_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 300)&(row22_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-22: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 3000)&(row22_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 3000)&(row22_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row22_AGE1[((row22_AGE1$a == 3000)&(row22_AGE1$b == 16)), 10] < 0.05)


##--- ROW-23

## AGE1-ROW-23: N=50
typeI_N50_P6 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 50)&(row23_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 50)&(row23_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 50)&(row23_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-23: N=100
typeI_N100_P6 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 100)&(row23_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 100)&(row23_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 100)&(row23_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-23: N=300
typeI_N300_P6 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 300)&(row23_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 300)&(row23_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 300)&(row23_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-23: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 3000)&(row23_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 3000)&(row23_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row23_AGE1[((row23_AGE1$a == 3000)&(row23_AGE1$b == 16)), 10] < 0.05)


##--- ROW-24

## AGE1-ROW-24: N=50
typeI_N50_P6 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 50)&(row24_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 50)&(row24_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 50)&(row24_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-24: N=100
typeI_N100_P6 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 100)&(row24_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 100)&(row24_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 100)&(row24_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-24: N=300
typeI_N300_P6 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 300)&(row24_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 300)&(row24_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 300)&(row24_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-24: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 3000)&(row24_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 3000)&(row24_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row24_AGE1[((row24_AGE1$a == 3000)&(row24_AGE1$b == 16)), 10] < 0.05)



##--- ROW-25

## AGE1-ROW-25: N=50
typeI_N50_P6 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 50)&(row25_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 50)&(row25_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 50)&(row25_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-25: N=100
typeI_N100_P6 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 100)&(row25_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 100)&(row25_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 100)&(row25_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-25: N=300
typeI_N300_P6 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 300)&(row25_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 300)&(row25_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 300)&(row25_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-25: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 3000)&(row25_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 3000)&(row25_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row25_AGE1[((row25_AGE1$a == 3000)&(row25_AGE1$b == 16)), 10] < 0.05)


##--- ROW-26

## AGE1-ROW-26: N=50
typeI_N50_P6 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 50)&(row26_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 50)&(row26_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 50)&(row26_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-26: N=100
typeI_N100_P6 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 100)&(row26_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 100)&(row26_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 100)&(row26_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-26: N=300
typeI_N300_P6 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 300)&(row26_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 300)&(row26_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 300)&(row26_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-26: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 3000)&(row26_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 3000)&(row26_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row26_AGE1[((row26_AGE1$a == 3000)&(row26_AGE1$b == 16)), 10] < 0.05)


##--- ROW-27

## AGE1-ROW-27: N=50
typeI_N50_P6 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 50)&(row27_AGE1$b == 6)), 10] < 0.05)
typeI_N50_P10 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 50)&(row27_AGE1$b == 10)), 10] < 0.05)
typeI_N50_P16 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 50)&(row27_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-27: N=100
typeI_N100_P6 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 100)&(row27_AGE1$b == 6)), 10] < 0.05)
typeI_N100_P10 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 100)&(row27_AGE1$b == 10)), 10] < 0.05)
typeI_N100_P16 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 100)&(row27_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-27: N=300
typeI_N300_P6 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 300)&(row27_AGE1$b == 6)), 10] < 0.05)
typeI_N300_P10 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 300)&(row27_AGE1$b == 10)), 10] < 0.05)
typeI_N300_P16 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 300)&(row27_AGE1$b == 16)), 10] < 0.05)

## AGE1-ROW-27: N=3000
typeI_N3000_P6 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 3000)&(row27_AGE1$b == 6)), 10] < 0.05)
typeI_N3000_P10 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 3000)&(row27_AGE1$b == 10)), 10] < 0.05)
typeI_N3000_P16 <- agresti_coull_CI(row27_AGE1[((row27_AGE1$a == 3000)&(row27_AGE1$b == 16)), 10] < 0.05)




## --------------------------
## Calculate Power: gammaSmall
## --------------------------

## --- BLR

## Table-1 [ROW-1]: ToIs=BLR; Commun=LOW; and, Corr=0.0
row1_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-2]: ToIs=BLR; Commun=LOW; and, Corr=0.4
row2_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-3]: ToIs=BLR; Commun=LOW; and, Corr=0.8
row3_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0.8), ]

## Table-1 [ROW-4]: ToIs=BLR; Commun=WIDE; and, Corr=0
row4_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-5]: ToIs=BLR; Commun=WIDE; and, Corr=0.4
row5_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-6]: ToIs=BLR; Commun=WIDE; and, Corr=0.8
row6_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0.8), ]

## Table-1 [ROW-7]: ToIs=BLR; Commun=HIGH; and, Corr=0.0
row7_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-8]: ToIs=BLR; Commun=HIGH; and, Corr=0.4
row8_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-9]: ToIs=BLR; Commun=HIGH; and, Corr=0.8
row9_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 1)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0.8), ]


## --- BMR

## Table-1 [ROW-10]: ToIs=BMR; Commun=LOW; and, Corr=0.0
row10_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-11]: ToIs=BMR; Commun=LOW; and, Corr=0.4
row11_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-12]: ToIs=BMR; Commun=LOW; and, Corr=0.8
row12_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0.8), ]

## Table-1 [ROW-13]: ToIs=BMR; Commun=WIDE; and, Corr=0.0
row13_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-14]: ToIs=BMR; Commun=WIDE; and, Corr=0.4
row14_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-15]: ToIs=BMR; Commun=WIDE; and, Corr=0.8
row15_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0.8), ]

## Table-1 [ROW-16]: ToIs=BMR; Commun=HIGH; and, Corr=0.0
row16_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-17]: ToIs=BMR; Commun=HIGH; and, Corr=0.4
row17_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-18]: ToIs=BMR; Commun=HIGH; and, Corr=0.8
row18_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 2)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0.8), ]


## --- Interval

## Table-1 [ROW-19]: ToIs=Interval; Commun=LOW; and, Corr=0.0
row19_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-20]: ToIs=Interval; Commun=LOW; and, Corr=0.4
row20_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-21]: ToIs=Interval; Commun=LOW; and, Corr=0.8
row21_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 1)&(powerGammaSmall$c == 0.8), ]

## Table-1 [ROW-22]: ToIs=Interval; Commun=WIDE; and, Corr=0.0
row22_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-23]: ToIs=Interval; Commun=WIDE; and, Corr=0.4
row23_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-24]: ToIs=Interval; Commun=WIDE; and, Corr=0.8
row24_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 2)&(powerGammaSmall$c == 0.8), ]

## Table-1 [ROW-25]: ToIs=Interval; Commun=HIGH; and, Corr=0
row25_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0), ]

## Table-1 [ROW-26]: ToIs=Interval; Commun=HIGH; and, Corr=0.4
row26_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0.4), ]

## Table-1 [ROW-27]: ToIs=Interval; Commun=HIGH; and, Corr=0.8
row27_gammaSmall <- powerGammaSmall[(powerGammaSmall$e == 3)&(powerGammaSmall$d == 3)&(powerGammaSmall$c == 0.8), ]





## ---------------------
## gammaModerate: Power
## ---------------------

## --- BLR

## Table-1 [ROW-1]: ToIs=BLR; Commun=LOW; and, Corr=0.0
row1_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-2]: ToIs=BLR; Commun=LOW; and, Corr=0.4
row2_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-3]: ToIs=BLR; Commun=LOW; and, Corr=0.8
row3_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0.8), ]

## Table-1 [ROW-4]: ToIs=BLR; Commun=WIDE; and, Corr=0
row4_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-5]: ToIs=BLR; Commun=WIDE; and, Corr=0.4
row5_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-6]: ToIs=BLR; Commun=WIDE; and, Corr=0.8
row6_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0.8), ]

## Table-1 [ROW-7]: ToIs=BLR; Commun=HIGH; and, Corr=0.0
row7_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-8]: ToIs=BLR; Commun=HIGH; and, Corr=0.4
row8_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-9]: ToIs=BLR; Commun=HIGH; and, Corr=0.8
row9_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 1)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0.8), ]


## --- BMR

## Table-1 [ROW-10]: ToIs=BMR; Commun=LOW; and, Corr=0.0
row10_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-11]: ToIs=BMR; Commun=LOW; and, Corr=0.4
row11_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-12]: ToIs=BMR; Commun=LOW; and, Corr=0.8
row12_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0.8), ]

## Table-1 [ROW-13]: ToIs=BMR; Commun=WIDE; and, Corr=0.0
row13_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-14]: ToIs=BMR; Commun=WIDE; and, Corr=0.4
row14_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-15]: ToIs=BMR; Commun=WIDE; and, Corr=0.8
row15_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0.8), ]

## Table-1 [ROW-16]: ToIs=BMR; Commun=HIGH; and, Corr=0.0
row16_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-17]: ToIs=BMR; Commun=HIGH; and, Corr=0.4
row17_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-18]: ToIs=BMR; Commun=HIGH; and, Corr=0.8
row18_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 2)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0.8), ]


## --- Interval

## Table-1 [ROW-19]: ToIs=Interval; Commun=LOW; and, Corr=0.0
row19_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-20]: ToIs=Interval; Commun=LOW; and, Corr=0.4
row20_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-21]: ToIs=Interval; Commun=LOW; and, Corr=0.8
row21_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 1)&(powerGammaModerate$c == 0.8), ]

## Table-1 [ROW-22]: ToIs=Interval; Commun=WIDE; and, Corr=0.0
row22_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-23]: ToIs=Interval; Commun=WIDE; and, Corr=0.4
row23_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-24]: ToIs=Interval; Commun=WIDE; and, Corr=0.8
row24_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 2)&(powerGammaModerate$c == 0.8), ]

## Table-1 [ROW-25]: ToIs=Interval; Commun=HIGH; and, Corr=0
row25_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0), ]

## Table-1 [ROW-26]: ToIs=Interval; Commun=HIGH; and, Corr=0.4
row26_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0.4), ]

## Table-1 [ROW-27]: ToIs=Interval; Commun=HIGH; and, Corr=0.8
row27_gammaModerate <- powerGammaModerate[(powerGammaModerate$e == 3)&(powerGammaModerate$d == 3)&(powerGammaModerate$c == 0.8), ]



## --------------------------
## Calculate Power: gammaBig
## --------------------------

## --- BLR

## Table-1 [ROW-1]: ToIs=BLR; Commun=LOW; and, Corr=0.0
row1_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-2]: ToIs=BLR; Commun=LOW; and, Corr=0.4
row2_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-3]: ToIs=BLR; Commun=LOW; and, Corr=0.8
row3_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0.8), ]

## Table-1 [ROW-4]: ToIs=BLR; Commun=WIDE; and, Corr=0
row4_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-5]: ToIs=BLR; Commun=WIDE; and, Corr=0.4
row5_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-6]: ToIs=BLR; Commun=WIDE; and, Corr=0.8
row6_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0.8), ]

## Table-1 [ROW-7]: ToIs=BLR; Commun=HIGH; and, Corr=0.0
row7_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-8]: ToIs=BLR; Commun=HIGH; and, Corr=0.4
row8_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-9]: ToIs=BLR; Commun=HIGH; and, Corr=0.8
row9_gammaBig <- powerGammaBig[(powerGammaBig$e == 1)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0.8), ]


## --- BMR

## Table-1 [ROW-10]: ToIs=BMR; Commun=LOW; and, Corr=0.0
row10_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-11]: ToIs=BMR; Commun=LOW; and, Corr=0.4
row11_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-12]: ToIs=BMR; Commun=LOW; and, Corr=0.8
row12_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0.8), ]

## Table-1 [ROW-13]: ToIs=BMR; Commun=WIDE; and, Corr=0.0
row13_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-14]: ToIs=BMR; Commun=WIDE; and, Corr=0.4
row14_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-15]: ToIs=BMR; Commun=WIDE; and, Corr=0.8
row15_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0.8), ]

## Table-1 [ROW-16]: ToIs=BMR; Commun=HIGH; and, Corr=0.0
row16_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-17]: ToIs=BMR; Commun=HIGH; and, Corr=0.4
row17_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-18]: ToIs=BMR; Commun=HIGH; and, Corr=0.8
row18_gammaBig <- powerGammaBig[(powerGammaBig$e == 2)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0.8), ]


## --- Interval

## Table-1 [ROW-19]: ToIs=Interval; Commun=LOW; and, Corr=0.0
row19_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-20]: ToIs=Interval; Commun=LOW; and, Corr=0.4
row20_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-21]: ToIs=Interval; Commun=LOW; and, Corr=0.8
row21_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 1)&(powerGammaBig$c == 0.8), ]

## Table-1 [ROW-22]: ToIs=Interval; Commun=WIDE; and, Corr=0.0
row22_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-23]: ToIs=Interval; Commun=WIDE; and, Corr=0.4
row23_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-24]: ToIs=Interval; Commun=WIDE; and, Corr=0.8
row24_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 2)&(powerGammaBig$c == 0.8), ]

## Table-1 [ROW-25]: ToIs=Interval; Commun=HIGH; and, Corr=0
row25_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0), ]

## Table-1 [ROW-26]: ToIs=Interval; Commun=HIGH; and, Corr=0.4
row26_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0.4), ]

## Table-1 [ROW-27]: ToIs=Interval; Commun=HIGH; and, Corr=0.8
row27_gammaBig <- powerGammaBig[(powerGammaBig$e == 3)&(powerGammaBig$d == 3)&(powerGammaBig$c == 0.8), ]



## -----------------------
## --- gammaSmall: BLR ---
## -----------------------

##--- ROW-1

## gammaSmall-ROW-1: N=50, no estimates because those models were NOT converged!!!
##    table(row1_gammaSmall$a): N=100(freq=88); 300(289); 3000 (300)

## powerSmall-ROW-1: N=100
power_N100_P6 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 100)&(row1_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 100)&(row1_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 100)&(row1_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-1: N=300
power_N300_P6 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 300)&(row1_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 300)&(row1_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 300)&(row1_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-1: N=3000
power_N3000_P6 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 3000)&(row1_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 3000)&(row1_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row1_gammaSmall[((row1_gammaSmall$a == 3000)&(row1_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-2

## gammaSmall-ROW-2: N=50, no estimates because those models were NOT converged!!!
##    table(row2_gammaSmall$a): N=100(freq=87); 300(287); 3000 (300)

## powerSmall-ROW-2: N=100
power_N100_P6 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 100)&(row2_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 100)&(row2_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 100)&(row2_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-2: N=300
power_N300_P6 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 300)&(row2_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 300)&(row2_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 300)&(row2_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-2: N=3000
power_N3000_P6 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 3000)&(row2_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 3000)&(row2_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row2_gammaSmall[((row2_gammaSmall$a == 3000)&(row2_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-3

## Neurot-ROW-3: N=50, no estimates because those models were NOT converged!!!
##    table(row3_gammaSmall$a): N=100(freq=94); 300(235); 3000 (236)

## powerSmall-ROW-3: N=100
power_N100_P6 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 100)&(row3_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 100)&(row3_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 100)&(row3_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-3: N=300
power_N300_P6 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 300)&(row3_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 300)&(row3_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 300)&(row3_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-3: N=3000
power_N3000_P6 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 3000)&(row3_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 3000)&(row3_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row3_gammaSmall[((row3_gammaSmall$a == 3000)&(row3_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-4

## gammaSmall-ROW-4: N=50, no estimates because those models were NOT converged!!!
##    table(row4_gammaSmall$a): N=100(freq=32); 300(223); 3000 (299)

## powerSmall-ROW-4: N=100
power_N100_P6 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 100)&(row4_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 100)&(row4_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 100)&(row4_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-4: N=300
power_N300_P6 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 300)&(row4_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 300)&(row4_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 300)&(row4_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-4: N=3000
power_N3000_P6 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 3000)&(row4_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 3000)&(row4_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row4_gammaSmall[((row4_gammaSmall$a == 3000)&(row4_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-5

## gammaSmall-ROW-5: N=50, no estimates because those models were NOT converged!!!
##    table(row5_gammaSmall$a): N=100(freq=45); 300(226); 3000 (299)

## powerSmall-ROW-5: N=100
power_N100_P6 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 100)&(row5_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 100)&(row5_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 100)&(row5_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-5: N=300
power_N300_P6 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 300)&(row5_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 300)&(row5_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 300)&(row5_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-5: N=3000
power_N3000_P6 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 3000)&(row5_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 3000)&(row5_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row5_gammaSmall[((row5_gammaSmall$a == 3000)&(row5_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-6

## powerSmall-ROW-6: N=50
power_N50_P6 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 50)&(row6_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 50)&(row6_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 50)&(row6_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-6: N=100
power_N100_P6 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 100)&(row6_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 100)&(row6_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 100)&(row6_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-6: N=300
power_N300_P6 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 300)&(row6_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 300)&(row6_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 300)&(row6_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-6: N=3000
power_N3000_P6 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 3000)&(row6_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 3000)&(row6_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row6_gammaSmall[((row6_gammaSmall$a == 3000)&(row6_gammaSmall$b == 16)), 10] < 0.05), 2)



##--- ROW-7

## powerSmall-ROW-7: N=50
power_N50_P6 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 50)&(row7_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 50)&(row7_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 50)&(row7_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-7: N=100
power_N100_P6 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 100)&(row7_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 100)&(row7_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 100)&(row7_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-7: N=300
power_N300_P6 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 300)&(row7_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 300)&(row7_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 300)&(row7_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-7: N=3000
power_N3000_P6 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 3000)&(row7_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 3000)&(row7_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row7_gammaSmall[((row7_gammaSmall$a == 3000)&(row7_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-8

## powerSmall-ROW-8: N=50
power_N50_P6 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 50)&(row8_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 50)&(row8_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 50)&(row8_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-8: N=100
power_N100_P6 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 100)&(row8_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 100)&(row8_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 100)&(row8_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-8: N=300
power_N300_P6 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 300)&(row8_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 300)&(row8_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 300)&(row8_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-8: N=3000
power_N3000_P6 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 3000)&(row8_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 3000)&(row8_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row8_gammaSmall[((row8_gammaSmall$a == 3000)&(row8_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-9

## powerSmall-ROW-9: N=50
power_N50_P6 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 50)&(row9_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 50)&(row9_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 50)&(row9_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-9: N=100
power_N100_P6 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 100)&(row9_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 100)&(row9_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 100)&(row9_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-9: N=300
power_N300_P6 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 300)&(row9_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 300)&(row9_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 300)&(row9_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-9: N=3000
power_N3000_P6 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 3000)&(row9_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 3000)&(row9_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row9_gammaSmall[((row9_gammaSmall$a == 3000)&(row9_gammaSmall$b == 16)), 10] < 0.05), 2)




## -----------------------
## --- gammaSmall: BMR ---
## -----------------------

##--- ROW-10

## powerSmall-ROW-10: N=50
power_N50_P6 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 50)&(row10_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 50)&(row10_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 50)&(row10_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-10: N=100
power_N100_P6 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 100)&(row10_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 100)&(row10_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 100)&(row10_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-10: N=300
power_N300_P6 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 300)&(row10_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 300)&(row10_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 300)&(row10_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-10: N=3000
power_N3000_P6 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 3000)&(row10_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 3000)&(row10_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row10_gammaSmall[((row10_gammaSmall$a == 3000)&(row10_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-11

## powerSmall-ROW-11: N=50
power_N50_P6 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 50)&(row11_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 50)&(row11_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 50)&(row11_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-11: N=100
power_N100_P6 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 100)&(row11_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 100)&(row11_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 100)&(row11_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-11: N=300
power_N300_P6 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 300)&(row11_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 300)&(row11_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 300)&(row11_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-11: N=3000
power_N3000_P6 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 3000)&(row11_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 3000)&(row11_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row11_gammaSmall[((row11_gammaSmall$a == 3000)&(row11_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-12

## powerSmall-ROW-12: N=50
power_N50_P6 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 50)&(row12_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 50)&(row12_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 50)&(row12_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-12: N=100
power_N100_P6 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 100)&(row12_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 100)&(row12_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 100)&(row12_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-12: N=300
power_N300_P6 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 300)&(row12_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 300)&(row12_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 300)&(row12_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-12: N=3000
power_N3000_P6 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 3000)&(row12_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 3000)&(row12_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row12_gammaSmall[((row12_gammaSmall$a == 3000)&(row12_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-13

## powerSmall-ROW-13: N=50
power_N50_P6 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 50)&(row13_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 50)&(row13_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 50)&(row13_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-13: N=100
power_N100_P6 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 100)&(row13_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 100)&(row13_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 100)&(row13_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-13: N=300
power_N300_P6 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 300)&(row13_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 300)&(row13_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 300)&(row13_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-13: N=3000
power_N3000_P6 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 3000)&(row13_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 3000)&(row13_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row13_gammaSmall[((row13_gammaSmall$a == 3000)&(row13_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-14

## powerSmall-ROW-14: N=50
power_N50_P6 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 50)&(row14_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 50)&(row14_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 50)&(row14_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-14: N=100
power_N100_P6 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 100)&(row14_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 100)&(row14_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 100)&(row14_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-14: N=300
power_N300_P6 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 300)&(row14_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 300)&(row14_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 300)&(row14_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-14: N=3000
power_N3000_P6 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 3000)&(row14_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 3000)&(row14_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row14_gammaSmall[((row14_gammaSmall$a == 3000)&(row14_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-15

## powerSmall-ROW-15: N=50
power_N50_P6 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 50)&(row15_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 50)&(row15_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 50)&(row15_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-15: N=100
power_N100_P6 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 100)&(row15_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 100)&(row15_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 100)&(row15_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-15: N=300
power_N300_P6 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 300)&(row15_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 300)&(row15_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 300)&(row15_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-15: N=3000
power_N3000_P6 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 3000)&(row15_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 3000)&(row15_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row15_gammaSmall[((row15_gammaSmall$a == 3000)&(row15_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-16

## powerSmall-ROW-16: N=50
power_N50_P6 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 50)&(row16_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 50)&(row16_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 50)&(row16_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-16: N=100
power_N100_P6 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 100)&(row16_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 100)&(row16_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 100)&(row16_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-16: N=300
power_N300_P6 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 300)&(row16_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 300)&(row16_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 300)&(row16_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-16: N=3000
power_N3000_P6 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 3000)&(row16_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 3000)&(row16_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row16_gammaSmall[((row16_gammaSmall$a == 3000)&(row16_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-17

## powerSmall-ROW-17: N=50
power_N50_P6 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 50)&(row17_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 50)&(row17_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 50)&(row17_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-17: N=100
power_N100_P6 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 100)&(row17_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 100)&(row17_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 100)&(row17_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-17: N=300
power_N300_P6 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 300)&(row17_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 300)&(row17_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 300)&(row17_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-17: N=3000
power_N3000_P6 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 3000)&(row17_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 3000)&(row17_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row17_gammaSmall[((row17_gammaSmall$a == 3000)&(row17_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-18

## powerSmall-ROW-18: N=50
power_N50_P6 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 50)&(row18_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 50)&(row18_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 50)&(row18_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-18: N=100
power_N100_P6 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 100)&(row18_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 100)&(row18_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 100)&(row18_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-18: N=300
power_N300_P6 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 300)&(row18_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 300)&(row18_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 300)&(row18_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-18: N=3000
power_N3000_P6 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 3000)&(row18_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 3000)&(row18_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row18_gammaSmall[((row18_gammaSmall$a == 3000)&(row18_gammaSmall$b == 16)), 10] < 0.05), 2)



## ----------------------------
## --- gammaSmall: Interval ---
## ----------------------------

##--- ROW-19

## powerSmall-ROW-19: N=50
power_N50_P6 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 50)&(row19_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 50)&(row19_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 50)&(row19_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-19: N=100
power_N100_P6 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 100)&(row19_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 100)&(row19_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 100)&(row19_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-19: N=300
power_N300_P6 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 300)&(row19_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 300)&(row19_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 300)&(row19_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-19: N=3000
power_N3000_P6 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 3000)&(row19_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 3000)&(row19_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row19_gammaSmall[((row19_gammaSmall$a == 3000)&(row19_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-20

## powerSmall-ROW-20: N=50
power_N50_P6 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 50)&(row20_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 50)&(row20_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 50)&(row20_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-20: N=100
power_N100_P6 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 100)&(row20_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 100)&(row20_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 100)&(row20_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-20: N=300
power_N300_P6 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 300)&(row20_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 300)&(row20_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 300)&(row20_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-20: N=3000
power_N3000_P6 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 3000)&(row20_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 3000)&(row20_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row20_gammaSmall[((row20_gammaSmall$a == 3000)&(row20_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-21

## powerSmall-ROW-21: N=50
power_N50_P6 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 50)&(row21_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 50)&(row21_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 50)&(row21_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-21: N=100
power_N100_P6 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 100)&(row21_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 100)&(row21_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 100)&(row21_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-21: N=300
power_N300_P6 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 300)&(row21_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 300)&(row21_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 300)&(row21_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-21: N=3000
power_N3000_P6 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 3000)&(row21_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 3000)&(row21_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row21_gammaSmall[((row21_gammaSmall$a == 3000)&(row21_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-22

## powerSmall-ROW-22: N=50
power_N50_P6 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 50)&(row22_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 50)&(row22_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 50)&(row22_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-22: N=100
power_N100_P6 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 100)&(row22_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 100)&(row22_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 100)&(row22_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-22: N=300
power_N300_P6 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 300)&(row22_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 300)&(row22_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 300)&(row22_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-22: N=3000
power_N3000_P6 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 3000)&(row22_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 3000)&(row22_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row22_gammaSmall[((row22_gammaSmall$a == 3000)&(row22_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-23

## powerSmall-ROW-23: N=50
power_N50_P6 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 50)&(row23_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 50)&(row23_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 50)&(row23_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-23: N=100
power_N100_P6 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 100)&(row23_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 100)&(row23_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 100)&(row23_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-23: N=300
power_N300_P6 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 300)&(row23_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 300)&(row23_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 300)&(row23_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-23: N=3000
power_N3000_P6 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 3000)&(row23_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 3000)&(row23_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row23_gammaSmall[((row23_gammaSmall$a == 3000)&(row23_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-24

## powerSmall-ROW-24: N=50
power_N50_P6 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 50)&(row24_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 50)&(row24_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 50)&(row24_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-24: N=100
power_N100_P6 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 100)&(row24_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 100)&(row24_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 100)&(row24_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-24: N=300
power_N300_P6 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 300)&(row24_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 300)&(row24_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 300)&(row24_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-24: N=3000
power_N3000_P6 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 3000)&(row24_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 3000)&(row24_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row24_gammaSmall[((row24_gammaSmall$a == 3000)&(row24_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-25

## powerSmall-ROW-25: N=50
power_N50_P6 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 50)&(row25_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 50)&(row25_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 50)&(row25_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-25: N=100
power_N100_P6 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 100)&(row25_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 100)&(row25_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 100)&(row25_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-25: N=300
power_N300_P6 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 300)&(row25_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 300)&(row25_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 300)&(row25_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-25: N=3000
power_N3000_P6 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 3000)&(row25_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 3000)&(row25_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row25_gammaSmall[((row25_gammaSmall$a == 3000)&(row25_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-26

## powerSmall-ROW-26: N=50
power_N50_P6 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 50)&(row26_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 50)&(row26_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 50)&(row26_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-26: N=100
power_N100_P6 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 100)&(row26_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 100)&(row26_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 100)&(row26_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-26: N=300
power_N300_P6 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 300)&(row26_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 300)&(row26_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 300)&(row26_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-26: N=3000
power_N3000_P6 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 3000)&(row26_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 3000)&(row26_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row26_gammaSmall[((row26_gammaSmall$a == 3000)&(row26_gammaSmall$b == 16)), 10] < 0.05), 2)


##--- ROW-27

## powerSmall-ROW-27: N=50
power_N50_P6 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 50)&(row27_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 50)&(row27_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 50)&(row27_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-27: N=100
power_N100_P6 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 100)&(row27_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 100)&(row27_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 100)&(row27_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-27: N=300
power_N300_P6 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 300)&(row27_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 300)&(row27_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 300)&(row27_gammaSmall$b == 16)), 10] < 0.05), 2)

## powerSmall-ROW-27: N=3000
power_N3000_P6 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 3000)&(row27_gammaSmall$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 3000)&(row27_gammaSmall$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row27_gammaSmall[((row27_gammaSmall$a == 3000)&(row27_gammaSmall$b == 16)), 10] < 0.05), 2)




## --------------------------
## --- gammaModerate: BLR ---
## --------------------------

##--- ROW-1

## gammaModerate-ROW-1: N=50, no estimates because those models were NOT converged!!!
##    table(row1_gammaModerate$a): N=100(freq=88); 300(289); 3000 (300)

## powerModerate-ROW-1: N=100
power_N100_P6 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 100)&(row1_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 100)&(row1_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 100)&(row1_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-1: N=300
power_N300_P6 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 300)&(row1_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 300)&(row1_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 300)&(row1_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-1: N=3000
power_N3000_P6 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 3000)&(row1_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 3000)&(row1_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row1_gammaModerate[((row1_gammaModerate$a == 3000)&(row1_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-2

## gammaModerate-ROW-2: N=50, no estimates because those models were NOT converged!!!
##    table(row2_gammaModerate$a): N=100(freq=87); 300(287); 3000 (300)

## powerModerate-ROW-2: N=100
power_N100_P6 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 100)&(row2_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 100)&(row2_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 100)&(row2_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-2: N=300
power_N300_P6 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 300)&(row2_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 300)&(row2_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 300)&(row2_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-2: N=3000
power_N3000_P6 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 3000)&(row2_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 3000)&(row2_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row2_gammaModerate[((row2_gammaModerate$a == 3000)&(row2_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-3

## gammaModerate-ROW-3: N=50, no estimates because those models were NOT converged!!!
##    table(row3_gammaModerate$a): N=100(freq=94); 300(285); 3000 (300)

## powerModerate-ROW-3: N=100
power_N100_P6 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 100)&(row3_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 100)&(row3_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 100)&(row3_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-3: N=300
power_N300_P6 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 300)&(row3_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 300)&(row3_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 300)&(row3_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-3: N=3000
power_N3000_P6 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 3000)&(row3_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 3000)&(row3_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row3_gammaModerate[((row3_gammaModerate$a == 3000)&(row3_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-4

## gammaModerate-ROW-4: N=50, no estimates because those models were NOT converged!!!
##    table(row4_gammaModerate$a): N=100(freq=32); 300(223); 3000 (299)

## powerModerate-ROW-4: N=100
power_N100_P6 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 100)&(row4_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 100)&(row4_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 100)&(row4_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-4: N=300
power_N300_P6 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 300)&(row4_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 300)&(row4_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 300)&(row4_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-4: N=3000
power_N3000_P6 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 3000)&(row4_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 3000)&(row4_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row4_gammaModerate[((row4_gammaModerate$a == 3000)&(row4_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-5

## gammaModerate-ROW-5: N=50, no estimates because those models were NOT converged!!!
##    table(row4_gammaModerate$a): N=100(freq=45); 300(226); 3000 (299)

## powerModerate-ROW-5: N=100
power_N100_P6 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 100)&(row5_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 100)&(row5_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 100)&(row5_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-5: N=300
power_N300_P6 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 300)&(row5_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 300)&(row5_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 300)&(row5_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-5: N=3000
power_N3000_P6 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 3000)&(row5_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 3000)&(row5_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row5_gammaModerate[((row5_gammaModerate$a == 3000)&(row5_gammaModerate$b == 16)), 10] < 0.05), 2)



##--- ROW-6

## powerModerate-ROW-6: N=50
power_N50_P6 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 50)&(row6_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 50)&(row6_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 50)&(row6_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-6: N=100
power_N100_P6 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 100)&(row6_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 100)&(row6_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 100)&(row6_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-6: N=300
power_N300_P6 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 300)&(row6_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 300)&(row6_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 300)&(row6_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-6: N=3000
power_N3000_P6 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 3000)&(row6_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 3000)&(row6_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row6_gammaModerate[((row6_gammaModerate$a == 3000)&(row6_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-7

## powerModerate-ROW-7: N=50
power_N50_P6 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 50)&(row7_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 50)&(row7_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 50)&(row7_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-7: N=100
power_N100_P6 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 100)&(row7_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 100)&(row7_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 100)&(row7_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-7: N=300
power_N300_P6 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 300)&(row7_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 300)&(row7_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 300)&(row7_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-7: N=3000
power_N3000_P6 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 3000)&(row7_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 3000)&(row7_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row7_gammaModerate[((row7_gammaModerate$a == 3000)&(row7_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-8

## powerModerate-ROW-8: N=50
power_N50_P6 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 50)&(row8_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 50)&(row8_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 50)&(row8_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-8: N=100
power_N100_P6 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 100)&(row8_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 100)&(row8_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 100)&(row8_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-8: N=300
power_N300_P6 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 300)&(row8_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 300)&(row8_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 300)&(row8_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-8: N=3000
power_N3000_P6 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 3000)&(row8_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 3000)&(row8_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row8_gammaModerate[((row8_gammaModerate$a == 3000)&(row8_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-9

## powerModerate-ROW-9: N=50
power_N50_P6 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 50)&(row9_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 50)&(row9_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 50)&(row9_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-9: N=100
power_N100_P6 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 100)&(row9_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 100)&(row9_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 100)&(row9_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-9: N=300
power_N300_P6 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 300)&(row9_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 300)&(row9_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 300)&(row9_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-9: N=3000
power_N3000_P6 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 3000)&(row9_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 3000)&(row9_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row9_gammaModerate[((row9_gammaModerate$a == 3000)&(row9_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-10

## powerModerate-ROW-10: N=50
power_N50_P6 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 50)&(row10_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 50)&(row10_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 50)&(row10_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-10: N=100
power_N100_P6 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 100)&(row10_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 100)&(row10_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 100)&(row10_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-10: N=300
power_N300_P6 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 300)&(row10_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 300)&(row10_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 300)&(row10_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-10: N=3000
power_N3000_P6 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 3000)&(row10_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 3000)&(row10_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row10_gammaModerate[((row10_gammaModerate$a == 3000)&(row10_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-11

## powerModerate-ROW-11: N=50
power_N50_P6 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 50)&(row11_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 50)&(row11_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 50)&(row11_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-11: N=100
power_N100_P6 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 100)&(row11_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 100)&(row11_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 100)&(row11_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-11: N=300
power_N300_P6 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 300)&(row11_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 300)&(row11_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 300)&(row11_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-11: N=3000
power_N3000_P6 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 3000)&(row11_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 3000)&(row11_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row11_gammaModerate[((row11_gammaModerate$a == 3000)&(row11_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-12

## powerModerate-ROW-12: N=50
power_N50_P6 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 50)&(row12_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 50)&(row12_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 50)&(row12_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-12: N=100
power_N100_P6 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 100)&(row12_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 100)&(row12_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 100)&(row12_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-12: N=300
power_N300_P6 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 300)&(row12_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 300)&(row12_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 300)&(row12_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-12: N=3000
power_N3000_P6 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 3000)&(row12_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 3000)&(row12_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row12_gammaModerate[((row12_gammaModerate$a == 3000)&(row12_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-13

## powerModerate-ROW-13: N=50
power_N50_P6 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 50)&(row13_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 50)&(row13_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 50)&(row13_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-13: N=100
power_N100_P6 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 100)&(row13_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 100)&(row13_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 100)&(row13_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-13: N=300
power_N300_P6 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 300)&(row13_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 300)&(row13_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 300)&(row13_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-13: N=3000
power_N3000_P6 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 3000)&(row13_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 3000)&(row13_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row13_gammaModerate[((row13_gammaModerate$a == 3000)&(row13_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-14

## powerModerate-ROW-14: N=50
power_N50_P6 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 50)&(row14_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 50)&(row14_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 50)&(row14_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-14: N=100
power_N100_P6 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 100)&(row14_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 100)&(row14_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 100)&(row14_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-14: N=300
power_N300_P6 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 300)&(row14_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 300)&(row14_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 300)&(row14_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-14: N=3000
power_N3000_P6 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 3000)&(row14_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 3000)&(row14_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row14_gammaModerate[((row14_gammaModerate$a == 3000)&(row14_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-15

## powerModerate-ROW-15: N=50
power_N50_P6 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 50)&(row15_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 50)&(row15_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 50)&(row15_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-15: N=100
power_N100_P6 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 100)&(row15_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 100)&(row15_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 100)&(row15_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-15: N=300
power_N300_P6 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 300)&(row15_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 300)&(row15_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 300)&(row15_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-15: N=3000
power_N3000_P6 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 3000)&(row15_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 3000)&(row15_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row15_gammaModerate[((row15_gammaModerate$a == 3000)&(row15_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-16

## powerModerate-ROW-16: N=50
power_N50_P6 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 50)&(row16_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 50)&(row16_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 50)&(row16_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-16: N=100
power_N100_P6 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 100)&(row16_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 100)&(row16_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 100)&(row16_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-16: N=300
power_N300_P6 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 300)&(row16_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 300)&(row16_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 300)&(row16_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-16: N=3000
power_N3000_P6 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 3000)&(row16_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 3000)&(row16_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row16_gammaModerate[((row16_gammaModerate$a == 3000)&(row16_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-17

## powerModerate-ROW-17: N=50
power_N50_P6 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 50)&(row17_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 50)&(row17_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 50)&(row17_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-17: N=100
power_N100_P6 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 100)&(row17_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 100)&(row17_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 100)&(row17_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-17: N=300
power_N300_P6 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 300)&(row17_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 300)&(row17_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 300)&(row17_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-17: N=3000
power_N3000_P6 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 3000)&(row17_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 3000)&(row17_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row17_gammaModerate[((row17_gammaModerate$a == 3000)&(row17_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-18

## powerModerate-ROW-18: N=50
power_N50_P6 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 50)&(row18_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 50)&(row18_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 50)&(row18_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-18: N=100
power_N100_P6 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 100)&(row18_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 100)&(row18_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 100)&(row18_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-18: N=300
power_N300_P6 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 300)&(row18_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 300)&(row18_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 300)&(row18_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-18: N=3000
power_N3000_P6 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 3000)&(row18_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 3000)&(row18_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row18_gammaModerate[((row18_gammaModerate$a == 3000)&(row18_gammaModerate$b == 16)), 10] < 0.05), 2)



## -------------------------------
## --- gammaModerate: Interval ---
## -------------------------------

##--- ROW-19

## powerModerate-ROW-19: N=50
power_N50_P6 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 50)&(row19_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 50)&(row19_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 50)&(row19_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-19: N=100
power_N100_P6 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 100)&(row19_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 100)&(row19_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 100)&(row19_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-19: N=300
power_N300_P6 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 300)&(row19_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 300)&(row19_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 300)&(row19_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-19: N=3000
power_N3000_P6 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 3000)&(row19_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 3000)&(row19_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row19_gammaModerate[((row19_gammaModerate$a == 3000)&(row19_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-20

## powerModerate-ROW-20: N=50
power_N50_P6 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 50)&(row20_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 50)&(row20_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 50)&(row20_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-20: N=100
power_N100_P6 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 100)&(row20_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 100)&(row20_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 100)&(row20_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-20: N=300
power_N300_P6 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 300)&(row20_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 300)&(row20_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 300)&(row20_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-20: N=3000
power_N3000_P6 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 3000)&(row20_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 3000)&(row20_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row20_gammaModerate[((row20_gammaModerate$a == 3000)&(row20_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-21

## powerModerate-ROW-21: N=50
power_N50_P6 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 50)&(row21_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 50)&(row21_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 50)&(row21_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-21: N=100
power_N100_P6 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 100)&(row21_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 100)&(row21_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 100)&(row21_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-21: N=300
power_N300_P6 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 300)&(row21_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 300)&(row21_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 300)&(row21_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-21: N=3000
power_N3000_P6 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 3000)&(row21_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 3000)&(row21_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row21_gammaModerate[((row21_gammaModerate$a == 3000)&(row21_gammaModerate$b == 16)), 10] < 0.05), 2)



##--- ROW-22

## powerModerate-ROW-22: N=50
power_N50_P6 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 50)&(row22_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 50)&(row22_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 50)&(row22_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-22: N=100
power_N100_P6 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 100)&(row22_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 100)&(row22_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 100)&(row22_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-22: N=300
power_N300_P6 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 300)&(row22_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 300)&(row22_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 300)&(row22_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-22: N=3000
power_N3000_P6 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 3000)&(row22_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 3000)&(row22_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row22_gammaModerate[((row22_gammaModerate$a == 3000)&(row22_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-23

## powerModerate-ROW-23: N=50
power_N50_P6 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 50)&(row23_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 50)&(row23_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 50)&(row23_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-23: N=100
power_N100_P6 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 100)&(row23_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 100)&(row23_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 100)&(row23_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-23: N=300
power_N300_P6 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 300)&(row23_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 300)&(row23_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 300)&(row23_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-23: N=3000
power_N3000_P6 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 3000)&(row23_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 3000)&(row23_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row23_gammaModerate[((row23_gammaModerate$a == 3000)&(row23_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-24

## powerModerate-ROW-24: N=50
power_N50_P6 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 50)&(row24_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 50)&(row24_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 50)&(row24_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-24: N=100
power_N100_P6 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 100)&(row24_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 100)&(row24_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 100)&(row24_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-24: N=300
power_N300_P6 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 300)&(row24_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 300)&(row24_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 300)&(row24_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-24: N=3000
power_N3000_P6 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 3000)&(row24_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 3000)&(row24_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row24_gammaModerate[((row24_gammaModerate$a == 3000)&(row24_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-25

## powerModerate-ROW-25: N=50
power_N50_P6 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 50)&(row25_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 50)&(row25_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 50)&(row25_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-25: N=100
power_N100_P6 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 100)&(row25_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 100)&(row25_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 100)&(row25_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-25: N=300
power_N300_P6 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 300)&(row25_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 300)&(row25_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 300)&(row25_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-25: N=3000
power_N3000_P6 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 3000)&(row25_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 3000)&(row25_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row25_gammaModerate[((row25_gammaModerate$a == 3000)&(row25_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-26

## powerModerate-ROW-26: N=50
power_N50_P6 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 50)&(row26_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 50)&(row26_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 50)&(row26_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-26: N=100
power_N100_P6 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 100)&(row26_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 100)&(row26_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 100)&(row26_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-26: N=300
power_N300_P6 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 300)&(row26_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 300)&(row26_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 300)&(row26_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-26: N=3000
power_N3000_P6 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 3000)&(row26_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 3000)&(row26_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row26_gammaModerate[((row26_gammaModerate$a == 3000)&(row26_gammaModerate$b == 16)), 10] < 0.05), 2)


##--- ROW-27

## powerModerate-ROW-27: N=50
power_N50_P6 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 50)&(row27_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 50)&(row27_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 50)&(row27_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-27: N=100
power_N100_P6 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 100)&(row27_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 100)&(row27_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 100)&(row27_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-27: N=300
power_N300_P6 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 300)&(row27_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 300)&(row27_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 300)&(row27_gammaModerate$b == 16)), 10] < 0.05), 2)

## powerModerate-ROW-27: N=3000
power_N3000_P6 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 3000)&(row27_gammaModerate$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 3000)&(row27_gammaModerate$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row27_gammaModerate[((row27_gammaModerate$a == 3000)&(row27_gammaModerate$b == 16)), 10] < 0.05), 2)



## ---------------------
## --- gammaBig: BLR ---
## ---------------------

##--- ROW-1

## Neurot-ROW-1: N=50, no estimates because those models were NOT converged!!!
##    table(row1_gammaBig$a): N=100(freq=88); 300(289); 3000 (300)

## powerBig-ROW-1: N=100
power_N100_P6 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 100)&(row1_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 100)&(row1_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 100)&(row1_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-1: N=300
power_N300_P6 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 300)&(row1_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 300)&(row1_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 300)&(row1_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-1: N=3000
power_N3000_P6 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 3000)&(row1_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 3000)&(row1_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row1_gammaBig[((row1_gammaBig$a == 3000)&(row1_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-2

## gammaBig-ROW-2: N=50, no estimates because those models were NOT converged!!!
##    table(row2_gammaBig$a): N=100(freq=87); 300(287); 3000 (300)

## powerBig-ROW-2: N=100
power_N100_P6 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 100)&(row2_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 100)&(row2_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 100)&(row2_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-2: N=300
power_N300_P6 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 300)&(row2_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 300)&(row2_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 300)&(row2_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-2: N=3000
power_N3000_P6 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 3000)&(row2_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 3000)&(row2_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row2_gammaBig[((row2_gammaBig$a == 3000)&(row2_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-3

## gammaBig-ROW-3: N=50, no estimates because those models were NOT converged!!!
##    table(row3_gammaBig$a): N=100(freq=94); 300(285); 3000 (300)

## powerBig-ROW-3: N=100
power_N100_P6 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 100)&(row3_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 100)&(row3_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 100)&(row3_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-3: N=300
power_N300_P6 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 300)&(row3_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 300)&(row3_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 300)&(row3_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-3: N=3000
power_N3000_P6 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 3000)&(row3_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 3000)&(row3_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row3_gammaBig[((row3_gammaBig$a == 3000)&(row3_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-4

## gammaBig-ROW-4: N=50, no estimates because those models were NOT converged!!!
##    table(row4_gammaBig$a): N=100(freq=32); 300(223); 3000 (299)

## powerBig-ROW-4: N=100
power_N100_P6 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 100)&(row4_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 100)&(row4_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 100)&(row4_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-4: N=300
power_N300_P6 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 300)&(row4_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 300)&(row4_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 300)&(row4_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-4: N=3000
power_N3000_P6 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 3000)&(row4_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 3000)&(row4_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row4_gammaBig[((row4_gammaBig$a == 3000)&(row4_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-5

## gammaBig-ROW-5: N=50, no estimates because those models were NOT converged!!!
##    table(row5_gammaBig$a): N=100(freq=45); 300(226); 3000 (299)

## powerBig-ROW-5: N=100
power_N100_P6 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 100)&(row5_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 100)&(row5_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 100)&(row5_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-5: N=300
power_N300_P6 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 300)&(row5_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 300)&(row5_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 300)&(row5_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-5: N=3000
power_N3000_P6 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 3000)&(row5_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 3000)&(row5_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row5_gammaBig[((row5_gammaBig$a == 3000)&(row5_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-6

## powerBig-ROW-6: N=50
power_N50_P6 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 50)&(row6_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 50)&(row6_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 50)&(row6_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-6: N=100
power_N100_P6 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 100)&(row6_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 100)&(row6_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 100)&(row6_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-6: N=300
power_N300_P6 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 300)&(row6_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 300)&(row6_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 300)&(row6_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-6: N=3000
power_N3000_P6 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 3000)&(row6_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 3000)&(row6_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row6_gammaBig[((row6_gammaBig$a == 3000)&(row6_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-7

## powerBig-ROW-7: N=50
power_N50_P6 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 50)&(row7_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 50)&(row7_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 50)&(row7_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-7: N=100
power_N100_P6 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 100)&(row7_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 100)&(row7_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 100)&(row7_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-7: N=300
power_N300_P6 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 300)&(row7_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 300)&(row7_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 300)&(row7_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-7: N=3000
power_N3000_P6 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 3000)&(row7_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 3000)&(row7_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row7_gammaBig[((row7_gammaBig$a == 3000)&(row7_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-8

## powerBig-ROW-8: N=50
power_N50_P6 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 50)&(row8_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 50)&(row8_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 50)&(row8_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-8: N=100
power_N100_P6 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 100)&(row8_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 100)&(row8_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 100)&(row8_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-8: N=300
power_N300_P6 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 300)&(row8_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 300)&(row8_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 300)&(row8_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-8: N=3000
power_N3000_P6 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 3000)&(row8_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 3000)&(row8_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row8_gammaBig[((row8_gammaBig$a == 3000)&(row8_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-9

## powerBig-ROW-9: N=50
power_N50_P6 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 50)&(row9_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 50)&(row9_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 50)&(row9_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-9: N=100
power_N100_P6 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 100)&(row9_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 100)&(row9_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 100)&(row9_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-9: N=300
power_N300_P6 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 300)&(row9_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 300)&(row9_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 300)&(row9_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-9: N=3000
power_N3000_P6 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 3000)&(row9_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 3000)&(row9_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row9_gammaBig[((row9_gammaBig$a == 3000)&(row9_gammaBig$b == 16)), 10] < 0.05), 2)


## ---------------------
## --- gammaBig: BMR ---
## ---------------------

##--- ROW-10

## powerBig-ROW-10: N=50
power_N50_P6 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 50)&(row10_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 50)&(row10_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 50)&(row10_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-10: N=100
power_N100_P6 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 100)&(row10_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 100)&(row10_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 100)&(row10_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-10: N=300
power_N300_P6 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 300)&(row10_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 300)&(row10_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 300)&(row10_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-10: N=3000
power_N3000_P6 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 3000)&(row10_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 3000)&(row10_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row10_gammaBig[((row10_gammaBig$a == 3000)&(row10_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-11

## powerBig-ROW-11: N=50
power_N50_P6 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 50)&(row11_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 50)&(row11_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 50)&(row11_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-11: N=100
power_N100_P6 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 100)&(row11_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 100)&(row11_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 100)&(row11_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-11: N=300
power_N300_P6 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 300)&(row11_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 300)&(row11_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 300)&(row11_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-11: N=3000
power_N3000_P6 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 3000)&(row11_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 3000)&(row11_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row11_gammaBig[((row11_gammaBig$a == 3000)&(row11_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-12

## powerBig-ROW-12: N=50
power_N50_P6 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 50)&(row12_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 50)&(row12_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 50)&(row12_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-12: N=100
power_N100_P6 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 100)&(row12_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 100)&(row12_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 100)&(row12_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-12: N=300
power_N300_P6 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 300)&(row12_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 300)&(row12_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 300)&(row12_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-12: N=3000
power_N3000_P6 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 3000)&(row12_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 3000)&(row12_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row12_gammaBig[((row12_gammaBig$a == 3000)&(row12_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-13

## powerBig-ROW-13: N=50
power_N50_P6 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 50)&(row13_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 50)&(row13_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 50)&(row13_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-13: N=100
power_N100_P6 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 100)&(row13_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 100)&(row13_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 100)&(row13_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-13: N=300
power_N300_P6 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 300)&(row13_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 300)&(row13_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 300)&(row13_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-13: N=3000
power_N3000_P6 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 3000)&(row13_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 3000)&(row13_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row13_gammaBig[((row13_gammaBig$a == 3000)&(row13_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-14

## powerBig-ROW-14: N=50
power_N50_P6 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 50)&(row14_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 50)&(row14_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 50)&(row14_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-14: N=100
power_N100_P6 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 100)&(row14_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 100)&(row14_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 100)&(row14_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-14: N=300
power_N300_P6 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 300)&(row14_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 300)&(row14_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 300)&(row14_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-14: N=3000
power_N3000_P6 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 3000)&(row14_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 3000)&(row14_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row14_gammaBig[((row14_gammaBig$a == 3000)&(row14_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-15

## powerBig-ROW-15: N=50
power_N50_P6 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 50)&(row15_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 50)&(row15_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 50)&(row15_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-15: N=100
power_N100_P6 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 100)&(row15_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 100)&(row15_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 100)&(row15_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-15: N=300
power_N300_P6 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 300)&(row15_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 300)&(row15_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 300)&(row15_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-15: N=3000
power_N3000_P6 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 3000)&(row15_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 3000)&(row15_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row15_gammaBig[((row15_gammaBig$a == 3000)&(row15_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-16

## powerBig-ROW-16: N=50
power_N50_P6 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 50)&(row16_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 50)&(row16_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 50)&(row16_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-16: N=100
power_N100_P6 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 100)&(row16_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 100)&(row16_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 100)&(row16_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-16: N=300
power_N300_P6 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 300)&(row16_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 300)&(row16_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 300)&(row16_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-16: N=3000
power_N3000_P6 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 3000)&(row16_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 3000)&(row16_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row16_gammaBig[((row16_gammaBig$a == 3000)&(row16_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-17

## powerBig-ROW-17: N=50
power_N50_P6 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 50)&(row17_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 50)&(row17_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 50)&(row17_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-17: N=100
power_N100_P6 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 100)&(row17_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 100)&(row17_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 100)&(row17_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-17: N=300
power_N300_P6 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 300)&(row17_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 300)&(row17_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 300)&(row17_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-17: N=3000
power_N3000_P6 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 3000)&(row17_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 3000)&(row17_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row17_gammaBig[((row17_gammaBig$a == 3000)&(row17_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-18

## powerBig-ROW-18: N=50
power_N50_P6 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 50)&(row18_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 50)&(row18_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 50)&(row18_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-18: N=100
power_N100_P6 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 100)&(row18_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 100)&(row18_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 100)&(row18_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-18: N=300
power_N300_P6 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 300)&(row18_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 300)&(row18_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 300)&(row18_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-18: N=3000
power_N3000_P6 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 3000)&(row18_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 3000)&(row18_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row18_gammaBig[((row18_gammaBig$a == 3000)&(row18_gammaBig$b == 16)), 10] < 0.05), 2)


## --------------------------
## --- gammaBig: Interval ---
## --------------------------

##--- ROW-19

## powerBig-ROW-19: N=50
power_N50_P6 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 50)&(row19_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 50)&(row19_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 50)&(row19_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-19: N=100
power_N100_P6 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 100)&(row19_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 100)&(row19_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 100)&(row19_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-19: N=300
power_N300_P6 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 300)&(row19_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 300)&(row19_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 300)&(row19_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-19: N=3000
power_N3000_P6 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 3000)&(row19_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 3000)&(row19_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row19_gammaBig[((row19_gammaBig$a == 3000)&(row19_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-20

## powerBig-ROW-20: N=50
power_N50_P6 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 50)&(row20_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 50)&(row20_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 50)&(row20_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-20: N=100
power_N100_P6 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 100)&(row20_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 100)&(row20_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 100)&(row20_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-20: N=300
power_N300_P6 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 300)&(row20_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 300)&(row20_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 300)&(row20_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-20: N=3000
power_N3000_P6 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 3000)&(row20_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 3000)&(row20_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row20_gammaBig[((row20_gammaBig$a == 3000)&(row20_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-21

## powerBig-ROW-21: N=50
power_N50_P6 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 50)&(row21_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 50)&(row21_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 50)&(row21_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-21: N=100
power_N100_P6 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 100)&(row21_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 100)&(row21_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 100)&(row21_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-21: N=300
power_N300_P6 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 300)&(row21_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 300)&(row21_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 300)&(row21_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-21: N=3000
power_N3000_P6 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 3000)&(row21_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 3000)&(row21_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row21_gammaBig[((row21_gammaBig$a == 3000)&(row21_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-22

## powerBig-ROW-22: N=50
power_N50_P6 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 50)&(row22_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 50)&(row22_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 50)&(row22_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-22: N=100
power_N100_P6 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 100)&(row22_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 100)&(row22_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 100)&(row22_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-22: N=300
power_N300_P6 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 300)&(row22_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 300)&(row22_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 300)&(row22_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-22: N=3000
power_N3000_P6 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 3000)&(row22_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 3000)&(row22_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row22_gammaBig[((row22_gammaBig$a == 3000)&(row22_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-23

## powerBig-ROW-23: N=50
power_N50_P6 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 50)&(row23_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 50)&(row23_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 50)&(row23_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-23: N=100
power_N100_P6 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 100)&(row23_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 100)&(row23_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 100)&(row23_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-23: N=300
power_N300_P6 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 300)&(row23_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 300)&(row23_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 300)&(row23_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-23: N=3000
power_N3000_P6 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 3000)&(row23_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 3000)&(row23_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row23_gammaBig[((row23_gammaBig$a == 3000)&(row23_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-24

## powerBig-ROW-24: N=50
power_N50_P6 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 50)&(row24_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 50)&(row24_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 50)&(row24_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-24: N=100
power_N100_P6 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 100)&(row24_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 100)&(row24_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 100)&(row24_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-24: N=300
power_N300_P6 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 300)&(row24_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 300)&(row24_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 300)&(row24_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-24: N=3000
power_N3000_P6 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 3000)&(row24_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 3000)&(row24_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row24_gammaBig[((row24_gammaBig$a == 3000)&(row24_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-25

## powerBig-ROW-25: N=50
power_N50_P6 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 50)&(row25_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 50)&(row25_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 50)&(row25_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-25: N=100
power_N100_P6 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 100)&(row25_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 100)&(row25_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 100)&(row25_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-25: N=300
power_N300_P6 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 300)&(row25_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 300)&(row25_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 300)&(row25_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-25: N=3000
power_N3000_P6 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 3000)&(row25_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 3000)&(row25_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row25_gammaBig[((row25_gammaBig$a == 3000)&(row25_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-26

## powerBig-ROW-26: N=50
power_N50_P6 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 50)&(row26_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 50)&(row26_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 50)&(row26_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-26: N=100
power_N100_P6 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 100)&(row26_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 100)&(row26_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 100)&(row26_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-26: N=300
power_N300_P6 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 300)&(row26_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 300)&(row26_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 300)&(row26_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-26: N=3000
power_N3000_P6 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 3000)&(row26_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 3000)&(row26_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row26_gammaBig[((row26_gammaBig$a == 3000)&(row26_gammaBig$b == 16)), 10] < 0.05), 2)


##--- ROW-27

## powerBig-ROW-27: N=50
power_N50_P6 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 50)&(row27_gammaBig$b == 6)), 10] < 0.05), 2)
power_N50_P10 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 50)&(row27_gammaBig$b == 10)), 10] < 0.05), 2)
power_N50_P16 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 50)&(row27_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-27: N=100
power_N100_P6 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 100)&(row27_gammaBig$b == 6)), 10] < 0.05), 2)
power_N100_P10 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 100)&(row27_gammaBig$b == 10)), 10] < 0.05), 2)
power_N100_P16 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 100)&(row27_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-27: N=300
power_N300_P6 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 300)&(row27_gammaBig$b == 6)), 10] < 0.05), 2)
power_N300_P10 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 300)&(row27_gammaBig$b == 10)), 10] < 0.05), 2)
power_N300_P16 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 300)&(row27_gammaBig$b == 16)), 10] < 0.05), 2)

## powerBig-ROW-27: N=3000
power_N3000_P6 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 3000)&(row27_gammaBig$b == 6)), 10] < 0.05), 2)
power_N3000_P10 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 3000)&(row27_gammaBig$b == 10)), 10] < 0.05), 2)
power_N3000_P16 <- round(mean(row27_gammaBig[((row27_gammaBig$a == 3000)&(row27_gammaBig$b == 16)), 10] < 0.05), 2)




