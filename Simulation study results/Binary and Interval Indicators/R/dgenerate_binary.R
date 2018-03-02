dgenerateBinary = function(N, P, rho_lat, lambda, prob, seed){
  # function for data generation
  # input
  #   N: sample size
  #   P: number of variables per factor
  #   rho.lat: correlation between the latent variables
  #   lambda: the range of factor loadings
  #   \copyright Worku, H. M. et al. October 2014.
  #--------------------------------------------------------------------------
  #--------------------------------------------------------------------------

  ## required packages
  require(mvtnorm)

  ld = length(lambda)
#   ld = 1      ## since comm=HIGH is skipped during the simulation
  le = length(prob)
#   le = 1

  for(a in N){
    for(b in P){
      for(c in rho_lat){
        for(d in 1:ld){
          for(e in 1:le){


            ## ----------------------------------------------------------
            ## STEP-1: generate the exogenous variables, X ~ N(0, sigma)
            ## ----------------------------------------------------------

            ## regression coefficients for exogenous variables, gamma (Mplus notation)
            gammaEDU1=-0.10; gammaEDU2=-0.20; gammaGEN1=gammaGEN2=0; gammaAGE1=gammaAGE2=0
            gammaN1=1.00; gammaN2=0.95; gammaE1=-0.30; gammaE2=-0.25;
            gammaA1=gammaA2=0; gammaC1=0; gammaC2=0.10; gammaO1=gammaO2=0

            ## predictors
            GEN = rbinom(a, 1, 0.67)
            X = rmvnorm(a, mean=matrix(0, 7, 1), sigma=diag(7), method="chol")
            colnames(X) = c("EDU", "AGE", "N", "E", "A", "C", "O")
            data0 = data.frame(GEN, X)


            ## --------------------------------
            ## STEP-2: simulate factor scores
            ## --------------------------------

            ## eta
            eta1 = gammaGEN1*data0$GEN + gammaEDU1*data0$EDU + gammaAGE1*data0$AGE + gammaN1*data0$N + gammaE1*data0$E + gammaA1*data0$A + gammaC1*data0$C + gammaO1*data0$O 
            eta2 = gammaGEN2*data0$GEN + gammaEDU2*data0$EDU + gammaAGE2*data0$AGE + gammaN2*data0$N + gammaE2*data0$E + gammaA2*data0$A + gammaC2*data0$C + gammaO2*data0$O

            ## specify the mean structure of the factor scores, i.e., predictors + residuals
            mu1 = eta1 + rnorm(a)
            mu2 = eta2 + rnorm(a)

            ## specify the variance-covariance structure of the factor scores   
            # variance
            psi11 = var(mu1)
            psi22 = var(mu2)
            psi_off  = c * sqrt(psi11*psi22)            

            psi_mat = matrix(c(psi11, psi_off, psi_off, psi22),2,2,byrow=T)            

#             theta = rmvnorm(a, mean=matrix(c(mean(XB1), mean(XB2)),2,1), psi_mat, method="chol")
            theta = matrix(0, a, 2)
            for(i in 1:a){
              theta[i, ] = rmvnorm(1, mean=matrix(c(mu1[i], mu2[i]),2,1), sigma=psi_mat, method="chol")
            }


            ## calculate the residuals for the indicators, i.e., unique factors
            # factor loadings
            lam_d = lambda[[d]]
#             lam_d = lambda      ## since comm=HIGH is skipped during the simulation  

            lam1 = round(runif(b/2, min=lam_d[1], max=lam_d[2]), 3)
            lam2 = round(runif(b/2, min=lam_d[1], max=lam_d[2]), 3)

            lambda_mat =  rbind(cbind(matrix(lam1,b/2,1), matrix(0,b/2,1)), 
                                cbind(matrix(0,b/2,1), matrix(lam2,b/2,1)))

#             phi_mat2 = diag(b) - diag(c(lam1^2, lam2^2))
            Phi_tmp = diag(lambda_mat%*%psi_mat%*%t(lambda_mat))
            phi_mat = diag(b) - diag(Phi_tmp)

            # unique factors
            epsil = rmvnorm(a, mean=matrix(0, b, 1), sigma=phi_mat, method="chol")            


            ## ----------------------------------------------------------
            ## STEP-3: Finally, obtain the 'continues' manifest variables
            ## ----------------------------------------------------------

            Ylat = theta%*%t(lambda_mat) + epsil

#             Ylat = rmvnorm(a, mean=matrix(c(0),b,1), Cov_Y, method="chol")


            ## ----------------------------------------------------------
            ## STEP-4: Finally, obtain the 'binary' manifest variables
            ## ----------------------------------------------------------

            # obtain the treshold, tau
            pi = prob[[e]]
#             pi = prob         ## since dichot=rare is skipped
            tau = numeric(b)
            YBin = matrix(c(0), a, b, byrow=FALSE)

            for(k in 1:b){
              tau[k] = qnorm(p = runif(1,min=pi[1],max=pi[2]), 
                          mean = mean(Ylat[ ,k]), sd = sd(Ylat[ ,k]), 
                          lower.tail = FALSE, log.p = FALSE)

              # Dichotomize the underlying variate, Ylat            
              YBin[ ,k] = ifelse(Ylat[ ,k] >= tau[k], 1, 0)
            }

            colnames(YBin) = paste("Y", 1:b, sep="")
            YBin = as.data.frame(YBin)

            XX = data0

            ## export
            fn <- paste(a,"-",b,"-",c,"-",d,"-",e,"-",seed,".RData",sep="")
            fp <- file.path(paste("./dataBinary/", fn))
            
            save(a, b, c, d, e, YBin, XX, theta, seed,
                 lam1 = lam1, lam2 = lam2, psi_off = psi_off, psi11, psi22,
                 file=fp)
          }
        }
      }
    }
  }
}