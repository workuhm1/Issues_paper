

##----------------------------------------------------------------------##
## Summarize simulation results focusing on measures of accuracy,       ##
## i.e., rate of Type-I error and Power of the Test.                    ##
##----------------------------------------------------------------------##
## Developer: Hailemichael M. Worku (aka, Haile)                        ##
## Date:      01-April-2017                                             ##
##----------------------------------------------------------------------##

## Environmental settings
setwd(dir = "/Volumes/HM_WORKU/LU/Issues paper/Simulation study_20170319/Codes_20170319/R_both (Binary and Interval Indicators)")

## Load library, needed for creating LaTex tables
# library(xtable)
library(tables)
library(Hmisc)
library(knitr)
library(tab)

## ------------------------------------------------- ##
## <!-- 0. Import simulation results (combined).  -- ##
## ------------------------------------------------- ##
sim_results <- read.table("bias_out_combined.txt", header=TRUE, sep=";")


## ------------------------------------------------- ##
## <!-- 1. Summarize Type-1 error rate results.   -- ##
## ------------------------------------------------- ##

agresti_coull_CI = function(dat){
## function for calculating Type-1 error proportion using Agresti and Coull method
  
  require(binom)
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

  
calculate_type1Rate <- function(target_stat, factor_nbr, which_pred) {
## function for calculating Type-1 error rate
  
  ## Get factor number
  factor_nbr_in <- paste("F", factor_nbr, ".", "ON", sep = '')
  
  ## Calculate Type-I error rate
  out_type1 <- data.frame(target            = character(), 
                          factor            = character(),
                          predictor         = character(),
                          row_nbr           = numeric(),
                          type_of_indicator = numeric(),
                          communality       = numeric(),
                          corr_f1f2         = numeric(),
                          sample_size       = numeric(),
                          Num_Ind_per_fact  = numeric(),
                          alpha_bar         = integer(),
                          alpha_lower       = integer(),
                          alpha_upper       = integer(), stringsAsFactors = FALSE)
  
  ## Get analysis dataset
  dset <- sim_results[(sim_results$paramHeader==factor_nbr_in) & (sim_results$param == which_pred), ]
  
  ## Initialize 
  row_nbr <- 0
  
  for (e_r in c(1, 2, 3)) {
    for (d_r in c(1, 2, 3)) {
      for (c_r in c(0.0, 0.4, 0.8)) {
        
        ## Get row-specific analysis dataset 
        dset_r <- dset[(dset$e == e_r)&(dset$d == d_r)&(dset$c == c_r), ]
        
        for (a_r in c(50, 100, 300, 3000)) {
          for (b_r in c(6, 10, 16)) {
            
            ## calculate cell specific Type-I error rate
            row_nbr <- row_nbr + 1
            
            type1_cell <- try(agresti_coull_CI(dset_r[((dset_r$a == a_r) & (dset_r$b == b_r)), 10] < 0.05),
                              silent = TRUE)
            
            if (class(type1_cell) == "try-error") {
              
              out_type1 <- rbind(out_type1, data.frame(target            = target_stat, 
                                                       factor            = factor_nbr_in,
                                                       predictor         = which_pred,
                                                       row_nbr           = row_nbr,
                                                       type_of_indicator = e_r,
                                                       communality       = d_r,
                                                       corr_f1f2         = c_r,
                                                       sample_size       = a_r,
                                                       Num_Ind_per_fact  = b_r,
                                                       alpha_bar         = 99,
                                                       alpha_lower       = -999,
                                                       alpha_upper       = 999) )
            } 
            
            if (class(type1_cell) != "try-error") { 
              
              out_type1 <- rbind(out_type1, data.frame(target            = target_stat, 
                                                       factor            = factor_nbr_in,
                                                       predictor         = which_pred,
                                                       row_nbr           = row_nbr,
                                                       type_of_indicator = e_r,
                                                       communality       = d_r,
                                                       corr_f1f2         = c_r,
                                                       sample_size       = a_r,
                                                       Num_Ind_per_fact  = b_r,
                                                       alpha_bar         = format(round(type1_cell$mean, 2), nsmall = 2),
                                                       alpha_lower       = round(type1_cell$lower, 3),
                                                       alpha_upper       = round(type1_cell$upper, 3)) )
            }
          }
        }
        
      }
    }
  }
  
  return(out_type1)
  
}

## ---------------------------------------------------------- ##
## AGE (on first factor): Calculate Type-1 error rate.        ##
## ---------------------------------------------------------- ##
out_type1_AGE <- calculate_type1Rate(target_stat = "Type-1", factor_nbr = 1, which_pred = "AGE")

## Format variable values
out_type1_AGE$type_of_indicator <- factor(out_type1_AGE$type_of_indicator, levels = c('1', '2', '3'),
                                          labels = c('BLR', 'BMR', 'Interval'))
out_type1_AGE$communality <- factor(out_type1_AGE$communality, levels = c('1', '2', '3'),
                                    labels = c('Low', 'Wide', 'High'))
out_type1_AGE$corr_f1f2 <- factor(out_type1_AGE$corr_f1f2, levels = c(0.0, 0.4, 0.8),
                                          labels = c('Independence', 'Moderate', 'Strong'))
out_type1_AGE$sample_size <- factor(out_type1_AGE$sample_size, levels = c(50, 100, 300, 3000),
                                  labels = c('50', '100', '300', '3000'))
out_type1_AGE$Num_Ind_per_fact <- factor(out_type1_AGE$Num_Ind_per_fact, levels = c(6, 10, 16),
                                  labels = c('6', '10', '16'))

## Export result
write.csv(out_type1_AGE, file = "results_type1_AGE.csv", row.names = FALSE)

## ----------------------- ##
## Create LaTex table      ##
## ------------------------##
# out_latex <- out_type1_AGE[ ,-c(1,2,3,4)]
# write.csv(out_type1_AGE, file = "results_type1_AGE_filtered.csv", row.names = FALSE)

# tbl <- ftable(out_latex$type_of_indicator, out_latex$communality, out_latex$corr_f1f2, out_latex$sample_size, out_latex$Num_Ind_per_fact, 
#               row.vars = c(1,2,3), dnn=c("Type of Indicators", "Communality", "Correlation between factors", "Sample Size", "Number of Indicators"))

# xtable(out_latex)

# latex(object = out_latex, 
#       cgroup = c("Sample Size", "Number of Indicators"), n.cgroup = c(4,5),
#       rgroup = c("Type of Indicators", "Communality", "Correlation between factors"), n.rgroup = c(1,2,3))
# 
# tabular(table = type_of_indicator + communality ~ alpha_bar,
#         data  = out_latex,
#         latex = TRUE)

# kable(head(iris), format = "latex")

# kable(x = out_latex, format = "latex")

## <!-- Not able to create LaTex table ---> ##
out_latex[(out_latex$type_of_indicator=="Interval") & (out_latex$communality=="High") & (out_latex$corr_f1f2=="Strong"), ]


## ------------------------------------------------- ##
## <!-- 2. Summarize Power of Test results.       -- ##
## ------------------------------------------------- ##

calculate_powerTest <- function(target_stat, effect_size) {
  
  ## Calculate Power of Test
  out_powerTest <- data.frame(target            = character(), 
                              factor            = character(),
                              predictor         = character(),
                              gammaTrue         = integer(),
                              row_nbr           = numeric(),
                              type_of_indicator = numeric(),
                              communality       = numeric(),
                              corr_f1f2         = numeric(),
                              sample_size       = numeric(),
                              Num_Ind_per_fact  = numeric(),
                              power_bar         = integer(), stringsAsFactors = FALSE)
  
  ## Get analysis dataset
  dset <- sim_results[sim_results$gammaTRUE == effect_size, ]
  
  ## Initialize 
  row_nbr <- 0
  
  for (e_r in c(1, 2, 3)) {
    for (d_r in c(1, 2, 3)) {
      for (c_r in c(0.0, 0.4, 0.8)) {
        
        ## Get row-specific analysis dataset 
        dset_r <- dset[(dset$e == e_r)&(dset$d == d_r)&(dset$c == c_r), ]
        
        for (a_r in c(50, 100, 300, 3000)) {
          for (b_r in c(6, 10, 16)) {
            
            ## calculate cell specific Power of Test
            row_nbr <- row_nbr + 1
            
            power_cell <- mean(dset_r[((dset_r$a == a_r)&(dset_r$b == b_r)), 10] < 0.05)
            
            if (is.na(power_cell)) {
              
              out_powerTest <- rbind(out_powerTest, data.frame(target            = target_stat, 
                                                               factor            = unique(dset_r$paramHeader),
                                                               predictor         = unique(dset_r$param),
                                                               gammaTrue         = unique(dset_r$gammaTRUE),
                                                               row_nbr           = row_nbr,
                                                               type_of_indicator = e_r,
                                                               communality       = d_r,
                                                               corr_f1f2         = c_r,
                                                               sample_size       = a_r,
                                                               Num_Ind_per_fact  = b_r,
                                                               power_bar         = 99 ))
            } 
            
            if (!is.na(power_cell)) { 
              
              out_powerTest <- rbind(out_powerTest, data.frame(target            = target_stat, 
                                                               factor            = unique(dset_r$paramHeader),
                                                               predictor         = unique(dset_r$param),
                                                               gammaTrue         = unique(dset_r$gammaTRUE),
                                                               row_nbr           = row_nbr,
                                                               type_of_indicator = e_r,
                                                               communality       = d_r,
                                                               corr_f1f2         = c_r,
                                                               sample_size       = a_r,
                                                               Num_Ind_per_fact  = b_r,
                                                               power_bar         = format(round(power_cell, 2), nsmall = 2) ))
            }
          }
        }
        
      }
    }
  }
  
  return(out_powerTest)
  
}


## ----------------------------------------------------------------- ##
## conscientious (on second factor): Calculate Power of Test.        ##
## ----------------------------------------------------------------- ##
out_power_CON <- calculate_powerTest(target_stat = "Power", effect_size = 0.10)

## ------------------------------ ##
## Export summarized results.     ##
## -------------------------------##
export_result <- function(out_object, out_name) {
  
  out <- out_object
  
  ## Format variable values
  out$type_of_indicator <- factor(out$type_of_indicator, levels = c('1', '2', '3'),
                                  labels = c('BLR', 'BMR', 'Interval'))
  
  out$communality <- factor(out$communality, levels = c('1', '2', '3'),
                            labels = c('Low', 'Wide', 'High'))
  
  out$corr_f1f2 <- factor(out$corr_f1f2, levels = c(0.0, 0.4, 0.8),
                          labels = c('Independence', 'Moderate', 'Strong'))
  
  ## Export result
  write.csv(out, file = paste(out_name, '.csv', sep = ''),
            row.names = FALSE)
  
}

export_result(out_object = out_power_CON, out_name = "results_power_CON")


## ----------------------------------------------------------------- ##
## Extraversion (on first factor): Calculate Power of Test.          ##
## ----------------------------------------------------------------- ##
out_power_EXT <- calculate_powerTest(target_stat = "Power", effect_size = -0.30)

export_result(out_object = out_power_EXT, out_name = "results_power_EXT")


## ----------------------------------------------------------------- ##
## Neuroticism (on second factor): Calculate Power of Test.          ##
## ----------------------------------------------------------------- ##
out_power_NEU <- calculate_powerTest(target_stat = "Power", effect_size = 0.95)

export_result(out_object = out_power_NEU, out_name = "results_power_NEU")




