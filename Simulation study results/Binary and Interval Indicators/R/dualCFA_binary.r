dualCFA <- function(summ, file, res,
                    seed, a, b, c, d, e,
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
    
    parms_Est = summ$parameters$ci.unstandardized[1:(b+1),
                                                  c("param" ,"est","low2.5","up2.5")] 
    
    ## ---------
    ## Recovery
    ## ---------
    
    fEst1 = summ$savedata[ ,(b+1)]
    fEst2 = summ$savedata[ ,(b+2)]
    
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
        #         fEst1 <- -fEst1
        #         fEst2 <- -fEst2
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
        #         fEst1 <- -fEst1
        #         fEst2 <- -fEst2
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
        #       fEst1 <- -fEst1
        #       fEst2 <- -fEst2
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
        #       fEst1 <- -fEst1
        #       fEst2 <- -fEst2
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
        #       fEst1 <- -fEst1
        #       fEst2 <- -fEst2
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
        #       fEst1 <- -fEst1
        #       fEst2 <- -fEst2
      }
    }
    
    # TRUE values
    parms_True = rbind(as.matrix(lam1),as.matrix(lam2),as.matrix(psi_off))
    colnames(parms_True) = "True"
    
    ## merge
    parms = data.frame(parms_Est, parms_True)  
    res2 = data.frame(seed, a, b, c, d, e, parms)
    
    fn <- paste("CFA", "_", "lambda", "_", file)
    fp <- file.path(paste("./CFA/lambda/", fn))
    
    write.table(x=res2, file=fp, row.names=FALSE)
    
  }
  
  fn2 <- paste("CFA", "_", "FI", "_", file)
  fp2 <- file.path(paste("./CFA/FI/", fn2))
  
  write.table(x=res, file=fp2, row.names=FALSE)
}


