dualModel <- function(dataY, dataXX, file, seed){
  ## fits both CFA and MIMIC model on the simulated data
  
  require(MplusAutomation)
  
  Y_names <- colnames(dataY)  
  X_names <- colnames(dataXX)  
  XY_names <- c(Y_names, X_names)
  dataXY <- data.frame(dataY, dataXX)
  
  
  ## The CFA (Measurement) Model
  model_spec1 = paste("\n f1 BY", paste(Y_names[1],"*",sep=""),
                      paste(Y_names[2:(b/2)],collapse=" "),
                      ";\n", 
                      "f2 BY", paste(Y_names[(b/2)+1],"*",sep=""),
                      paste(Y_names[((b/2)+2):b],collapse=" "),
                      ";\n", 
                      "f1@1",
                      ";\n", 
                      "f2@1",
                      ";", sep=" ")
  
  ## The MIMIC model, with exogenous variables
  model_spec2 = paste("\n f1 BY", paste(Y_names[1],"*",sep=""),
                      paste(Y_names[2:(b/2)],collapse=" "),
                      ";\n", 
                      "f2 BY", paste(Y_names[(b/2)+1],"*",sep=""),
                      paste(Y_names[((b/2)+2):b],collapse=" "),
                      ";\n", 
                      "f1@1",
                      ";\n", 
                      "f2@1",
                      ";", 
                      "\n f1 ON", paste(X_names[1:length(X_names)],collapse=" "),
                      ";\n", 
                      "f2 ON", paste(X_names[1:length(X_names)],collapse=" "),
                      ";", 
                      sep=" ")  
  
  ## Specify the CFA model
  mod1 <- mplusObject(
    TITLE = "fit the CFA Model;",
    ANALYSIS = "
      ESTIMATOR = MLR;
      ITERATIONS = 10000;
      STARTS = 2000;",
    MODEL = 
      model_spec1,
    OUTPUT = "
      TECH4 STANDARDIZED MODINDICES(4) RES CINTERVAL;",
    SAVEDATA = "
      FILE is scores.txt;
      SAVE = fscores;",
    usevariables = Y_names,
    rdata = dataY)
  
  ## Specify the MIMIC model in Mplus
  mod2 <- mplusObject(
    TITLE = "fit the MIMIC Model;",
    ANALYSIS = "
      ESTIMATOR = MLR;
      ITERATIONS = 10000;
      STARTS = 2000;",
    MODEL = 
      model_spec2,
    OUTPUT = "
      TECH4 STANDARDIZED MODINDICES(4) RES CINTERVAL;",
    SAVEDATA = "
      FILE is scores.txt;
      SAVE = fscores;",
    usevariables = XY_names,
    rdata = dataXY)
  
  ## fit CFA model
  fit1 <- mplusModeler(object=mod1, 
                       dataout="dataCFA.dat",  
                       modelout = "CFA.inp", 
                       run = 1L)  
  
  ## extract the results of CFA
  summCFA = readModels()  
  unlink("CFA.out")
  unlink("CFA.inp")
  unlink("dataCFA.dat")
  
  ## fit the MIMIC model
  fit2 <- mplusModeler(object=mod2, 
                       dataout="dataMIMIC.dat",  
                       modelout = "MIMIC.inp", 
                       run = 1L)  
  
  ## extract the results from the SEM
  summMIMIC = readModels()  
  unlink("MIMIC.out")
  unlink("MIMIC.inp")
  unlink("dataMIMIC.dat")
  
  ## save outputs, CFA
  res_CFA <- data.frame(seed=seed, a=a, b=b, c=c, d=d, 
                        rho1=numeric(1), rho2=numeric(1), 
                        err1=logical(1), err2=logical(1),
                        flipped=character(1))
  
  ## save outputs, MIMIC
  res_MIMIC <- data.frame(seed=seed, a=a, b=b, c=c, d=d, 
                          rho1=numeric(1), rho2=numeric(1), 
                          err1=logical(1), err2=logical(1),
                          flipped=character(1))  
  
  ## heywood cases: if no standard errors (not positive definite covariance matrix)
  prms_unstd_CFA = summCFA$parameters$unstandardized
  res_CFA$err1 = (!any(colnames(prms_unstd_CFA) == "se"))
  
  prms_unstd_MIMIC = summMIMIC$parameters$unstandardized
  res_MIMIC$err1 = (!any(colnames(prms_unstd_MIMIC) == "se"))
  
  
  ## nonconvergence
  res_CFA$err2 = any(is.na(summCFA$savedata))
  res_MIMIC$err2 = any(is.na(summMIMIC$savedata))
  
  ## calculate the correlation btn theta_CFA and theta_MIMIC
  if(((res_CFA$err1 == FALSE) & (res_CFA$err2 == FALSE))&
     ((res_MIMIC$err1 == FALSE) & (res_MIMIC$err2 == FALSE))){
    
    res_corrCFAMIMIC <- data.frame(seed=seed, a=a, b=b, c=c, d=d, 
                                   rho1_CFAMIMIC=numeric(1), rho2_CFAMIMIC=numeric(1))  
    
    ## CFA
    f1_CFA = summCFA$savedata[ ,(b+1)]
    f2_CFA = summCFA$savedata[ ,(b+3)]
    
    ## MIMIC
    f1_MIMIC = summMIMIC$savedata[ ,(b+9)]
    f2_MIMIC = summMIMIC$savedata[ ,(b+11)]
    
    ## correlations
    res_corrCFAMIMIC$rho1_CFAMIMIC <- round(cor(f1_CFA, f1_MIMIC), 3)
    res_corrCFAMIMIC$rho2_CFAMIMIC <- round(cor(f2_CFA, f2_MIMIC), 3)
    
    ## export the results
    fn <- paste("corrCFAMIMIC", "_", file)
    fp <- file.path(paste("./corr_ThetaCFA_ThetaMIMIC/", fn))
    
    write.table(x=res_corrCFAMIMIC, file=fp, row.names=FALSE)
  }
  
  ## call CFA
  dualCFA(summCFA, file, res_CFA,
          seed, a, b, c, d,
          theta, lam1, lam2, psi_off)
  
  ## call MIMIC
  dualMIMIC(summMIMIC, file, res_MIMIC,
            seed, a, b, c, d,
            theta, lam1, lam2, psi_off)
  
  
  
}