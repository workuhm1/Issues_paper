dualMIMIC <- function(summ, file, res,
                      seed, a, b, c, d, 
                      theta, lam1, lam2, psi_off){
  
  
  ## save recovery
  if((res$err1==TRUE) | (res$err2==TRUE)){    
    #i = i+1; cat("i=",i,"\n")
    
    res$rho1 = NA
    res$rho2 = NA   
    
    res$flipped = "NOT-IMP"
    
  }else{
    #j = j+1; cat("j=",j,"\n")
    
    ## ------------------------------------------------ 
    ## bias-MSE-coverage for # lambda & psi12_off=corr    
    ## ------------------------------------------------
    
    parms_Est = summ$parameters$ci.unstandardized[c(1:b,(b+17)),
                                                  c("param" ,"est","low2.5","up2.5")] 
    
    parms_Est2 = summ$parameters$unstandardized[c((b+1):(b+16)),
                                                c("paramHeader", "param", "est", "pval")]
    
    ## ---------
    ## Recovery
    ## ---------
    
    fEst1 = summ$savedata[ ,(b+9)]
    fEst2 = summ$savedata[ ,(b+11)]
    
    # Merge 
    F1 = cbind(theta[,1], fEst1)
    F2 = cbind(theta[,2], fEst2)
    
    # correlation
    res$rho1 = round(cor(F1)[1,2], 3)
    res$rho2 = round(cor(F2)[1,2], 3)
    
    res$flipped = "NOT-flipped"
    
    ## --------------------------------------------------
    ## Checks if the factor loadings are flipped: START
    ## --------------------------------------------------
    
    latParm1 = parms_Est[1:(b/2), ]
    latParm2 = parms_Est[((b/2)+1):b, ]
    
    ## number of negatives (nn)
    nn1 = sum(latParm1$est < 0)
    nn2 = sum(latParm2$est < 0)
    
    bb = b/2
    
    if(bb == 3){
      if(nn1 >= (bb-1)){
        
        # flip the lambda estimate
        tmp_est = parms_Est[1:(b/2), "est"]
        parms_Est[1:(b/2), "est"] = -tmp_est
        
        # flip CIs
        tmp_low2.5 = parms_Est[1:(b/2), "low2.5"]
        tmp_up2.5  = parms_Est[1:(b/2), "up2.5"]
        parms_Est[1:(b/2), "low2.5"] = -tmp_up2.5
        parms_Est[1:(b/2), "up2.5"]  = -tmp_low2.5
        
        # correlation
        tmp_rho1 = res$rho1
        res$rho1 = -tmp_rho1
        #         tmp_rho2 = res$rho2
        #         res$rho2 = -tmp_rho2
        
        res$flipped = "flipped"
      }
      
      if (nn2 >= (bb-1)){
        
        # flip the lambda estimate
        tmp_est = parms_Est[((b/2)+1):b, "est"]
        parms_Est[((b/2)+1):b, "est"] = -tmp_est
        
        # flip CIs
        tmp_low2.5 = parms_Est[((b/2)+1):b, "low2.5"]
        tmp_up2.5  = parms_Est[((b/2)+1):b, "up2.5"]
        parms_Est[((b/2)+1):b, "low2.5"] = -tmp_up2.5
        parms_Est[((b/2)+1):b, "up2.5"]  = -tmp_low2.5
        
        # correlation
        #         tmp_rho1 = res$rho1
        #         res$rho1 = -tmp_rho1
        tmp_rho2 = res$rho2
        res$rho2 = -tmp_rho2
        
        res$flipped = "flipped"
        
      }
    }
    else if(bb == 5){
      if(nn1 >= (bb-2)){
        
        # flip the lambda estimate
        tmp_est = parms_Est[1:(b/2), "est"]
        parms_Est[1:(b/2), "est"] = -tmp_est
        
        # flip CIs
        tmp_low2.5 = parms_Est[1:(b/2), "low2.5"]
        tmp_up2.5  = parms_Est[1:(b/2), "up2.5"]
        parms_Est[1:(b/2), "low2.5"] = -tmp_up2.5
        parms_Est[1:(b/2), "up2.5"]  = -tmp_low2.5
        
        # correlation
        tmp_rho1 = res$rho1
        res$rho1 = -tmp_rho1
        #       tmp_rho2 = res$rho2
        #       res$rho2 = -tmp_rho2
        
        res$flipped = "flipped"
      }
      
      if (nn2 >= (bb-2)){
        
        # flip the lambda estimate
        tmp_est = parms_Est[((b/2)+1):b, "est"]
        parms_Est[((b/2)+1):b, "est"] = -tmp_est
        
        # flip CIs
        tmp_low2.5 = parms_Est[((b/2)+1):b, "low2.5"]
        tmp_up2.5  = parms_Est[((b/2)+1):b, "up2.5"]
        parms_Est[((b/2)+1):b, "low2.5"] = -tmp_up2.5
        parms_Est[((b/2)+1):b, "up2.5"]  = -tmp_low2.5
        
        # correlation
        #       tmp_rho1 = res$rho1
        #       res$rho1 = -tmp_rho1
        tmp_rho2 = res$rho2
        res$rho2 = -tmp_rho2
        
        res$flipped = "flipped"
        
      }
    }
    else{
      if(nn1 >= (bb-3)){
        
        # flip the lambda estimate
        tmp_est = parms_Est[1:(b/2), "est"]
        parms_Est[1:(b/2), "est"] = -tmp_est
        
        # flip CIs
        tmp_low2.5 = parms_Est[1:(b/2), "low2.5"]
        tmp_up2.5  = parms_Est[1:(b/2), "up2.5"]
        parms_Est[1:(b/2), "low2.5"] = -tmp_up2.5
        parms_Est[1:(b/2), "up2.5"]  = -tmp_low2.5
        
        # correlation
        tmp_rho1 = res$rho1
        res$rho1 = -tmp_rho1
        #       tmp_rho2 = res$rho2
        #       res$rho2 = -tmp_rho2
        
        res$flipped = "flipped"
      }
      
      if(nn2 >= (bb-3)){
        
        # flip the lambda estimate
        tmp_est = parms_Est[((b/2)+1):b, "est"]
        parms_Est[((b/2)+1):b, "est"] = -tmp_est
        
        # flip CIs
        tmp_low2.5 = parms_Est[((b/2)+1):b, "low2.5"]
        tmp_up2.5  = parms_Est[((b/2)+1):b, "up2.5"]
        parms_Est[((b/2)+1):b, "low2.5"] = -tmp_up2.5
        parms_Est[((b/2)+1):b, "up2.5"]  = -tmp_low2.5
        
        # correlation
        #       tmp_rho1 = res$rho1
        #       res$rho1 = -tmp_rho1
        tmp_rho2 = res$rho2
        res$rho2 = -tmp_rho2
        
        res$flipped = "flipped"
        
      }
    }
    
    # TRUE values
    parms_True = rbind(as.matrix(lam1),as.matrix(lam2),as.matrix(psi_off))
    colnames(parms_True) = "True"
    
    ## merge
    parms = data.frame(parms_Est, parms_True)    
    
    res2 = data.frame(seed, a, b, c, d, parms)  
    
    ## export
    # write.table(x=res2, file=paste("SEM", "_", "lambda", "_", file), row.names=FALSE)
    fn <- paste("SEM", "_", "lambda", "_", file)
    fp <- file.path(paste("./MIMIC/lambda/", fn))
    
    write.table(x=res2, file=fp, row.names=FALSE)
    
    ## save the gamma's, the regression weights of exogeneous variables
    gammaTRUE = matrix(c(0.00,-0.10,0.00,1.00,-0.30,0.00,0.00,0.00,
                         0.00,-0.20,0.00,0.95,-0.25,0.00,0.10,0.00))
    gammaTRUE = data.frame(gammaTRUE=gammaTRUE)
    
    res3 = data.frame(seed, a, b, c, d, parms_Est2, gammaTRUE)  
    
    ## export
    # write.table(x=res3, file=paste("SEM", "_", "beta", "_", file), row.names=FALSE)
    fn2 <- paste("SEM", "_", "gamma", "_", file)
    fp2 <- file.path(paste("./MIMIC/gamma/", fn2))
    
    write.table(x=res3, file=fp2, row.names=FALSE)
    
  }
  
  #   unlink("data.dat")
  #   unlink("data.inp")
  #   unlink("data.out")
  
  ## export
  # write.table(x=res, file=paste("SEM", "_", "FI", "_", file), row.names=FALSE)
  fn3 <- paste("SEM", "_", "FI", "_", file)
  fp3 <- file.path(paste("./MIMIC/FI/", fn3))
  
  write.table(x=res, file=fp3, row.names=FALSE)
  
}