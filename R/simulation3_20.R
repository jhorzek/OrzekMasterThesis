#' Simulation study 3
#' 20 percent of cross-lagged parameters non-zero
#'
#'
#' @param sampleSize sampleSize
#' @param seed seed
#' @param wd working directory. The results will be saved here
#' @param total_repetitions total number of repetitions
#' @param crossEffect size of cross-lagged effects
#' @author Jannik Orzek
#' @import OpenMx caret laremm regsem tcltk
#'
#' @examples
#'
#' @export
#'
#'
simulation3_20 <- function(sampleSize, seed, wd, total_repetitions, crossEffect, autoEffect){
  setwd(wd)
  ##### Settings #####

  numLatent = 5
  Timepoints = 5
  burning = 5
  set.seed(seed)

  ##### Parameters #####

  Avalues <- diag(.5,nrow=5,ncol = 5)
  Avalues[1,2] <-crossEffect
  Avalues[2,3] <-crossEffect
  Avalues[3,4] <-crossEffect
  Avalues[4,5] <-crossEffect
  Avalues[5,1] <-crossEffect

  Afree <- matrix(TRUE, 5,5)
  Alabel <- matrix(c("a11", "a12", "a13", "a14", "a15", "a21", "a22", "a23", "a24", "a25", "a31", "a32", "a33", "a34", "a35", "a41", "a42", "a43", "a44", "a45", "a51", "a52", "a53", "a54", "a55"), nrow = 5, byrow = T)

  Svalues <- diag(1-.5^2,nrow=5,ncol = 5)
  Svalues[1,1] <- 1-.5^2 -crossEffect^2
  Svalues[2,2] <- 1-.5^2 -crossEffect^2
  Svalues[3,3] <- 1-.5^2 -crossEffect^2
  Svalues[4,4] <- 1-.5^2 -crossEffect^2
  Svalues[5,5] <- 1-.5^2 -crossEffect^2

  Sfree <- diag(TRUE, 5,5)
  Slabel <- matrix(c("s11", NA, NA, NA, NA,
                     NA, "s22", NA, NA, NA,
                     NA, NA, "s33", NA, NA,
                     NA, NA, NA, "s44", NA,
                     NA, NA, NA, NA, "s55"), nrow = 5, byrow = T)

  ##### select cross-lagged values for regularization #####

  penAelement <- matrix(1, 5,5)
  diag(penAelement) <- 0
  OBS_manifestVarNames <- c()
  OBS_latentVarNames <- c()
  for(i in (burning+1):(Timepoints+burning)){
    OBS_manifestVarNames <- c(OBS_manifestVarNames, paste("x", 1:numLatent,"_t",i, sep=""))
    OBS_latentVarNames <- c(OBS_latentVarNames, paste("eta", 1:numLatent,"_t",i, sep=""))
  }
  penAFull <- matrix(0,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  rownames(penAFull) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(penAFull) <- c(OBS_manifestVarNames, OBS_latentVarNames)

  for (row in seq(which(rownames(penAFull) == paste("eta1_t",burning+2,sep = "")), which(rownames(penAFull) == paste("eta1_t", (Timepoints+burning), sep = "")), by = numLatent)){

    penAFull[row:(row+numLatent-1),(row-numLatent):(row-1)] <- penAelement

  }


  ##### prepare result data sets #####
  total_repetitions = total_repetitions
  iteration = 1
  improper_solutions = 0
  error = 0

  overall_evaluation <- data.frame("AIC_false_nonzero" = rep(NA,total_repetitions),
                                   "AIC_false_zero"=rep(NA,total_repetitions) ,
                                   "FIML_AIC_false_nonzero" = rep(NA,total_repetitions),
                                   "FIML_AIC_false_zero"=rep(NA,total_repetitions) ,

                                   "BIC_false_nonzero"=rep(NA,total_repetitions),
                                   "BIC_false_zero"=rep(NA,total_repetitions),
                                   "FIML_BIC_false_nonzero"=rep(NA,total_repetitions),
                                   "FIML_BIC_false_zero"=rep(NA,total_repetitions),

                                   "CV_m2LL_false_nonzero"=rep(NA,total_repetitions),
                                   "CV_m2LL_false_zero"=rep(NA,total_repetitions),
                                   "FIML_CV_m2LL_false_nonzero"=rep(NA,total_repetitions),
                                   "FIML_CV_m2LL_false_zero"=rep(NA,total_repetitions)
  )

  parameter_Table <- matrix(NA, nrow = length(unique(Alabel[Afree])), ncol = total_repetitions)

  rownames(parameter_Table) <- c(unique(Alabel[Afree]))

  parameterValues <- list("Base_Model_cov" = parameter_Table,
                          "Base_Model_FIML" = parameter_Table,
                          "AIC"=parameter_Table, "BIC"=parameter_Table,
                          "FIML_BIC"=parameter_Table, "FIML_AIC"=parameter_Table,
                          "CV_m2LL" = parameter_Table,
                          "FIML_CV_m2LL" =parameter_Table)


  RMSE_table <- matrix(NA, nrow = 3, ncol = total_repetitions)
  rownames(RMSE_table) <- c("RMSE_Amat", "RMSE_Smat", "RMSE_comb")
  RMSE <- list("Base_Model_cov" = RMSE_table,
               "Base_Model_FIML" = RMSE_table,
               "AIC"=RMSE_table, "BIC"=RMSE_table,
               "FIML_BIC"=RMSE_table, "FIML_AIC"=RMSE_table,
               "CV_m2LL" = RMSE_table,
               "FIML_CV_m2LL" =RMSE_table)
  ##### start simulation

  # Progress bar
  global_pb <- tkProgressBar(title = "progress bar", min = 0,
                             max = total_repetitions, width = 300) # progress - bar

  while(iteration <= total_repetitions){

    # simulate ARCL-SEM:
    simObj <- simARCL_scaled_20(seed = seed, sampleSize = sampleSize, autoEffect = autoEffect, crossEffect = crossEffect)

    # extract and scale raw data
    full_raw_data <- simObj$mxARCL_FIML$data$observed

    # split dataset in train and test dataset
    Folds <- createFolds(c(1:sampleSize),2)

    train_raw_data <- full_raw_data[Folds$Fold1,]
    test_raw_data <- full_raw_data[Folds$Fold2,]

    train_raw_data <- scale(train_raw_data)
    test_raw_data <- scale(test_raw_data)

    full_raw_data <- scale(full_raw_data)

    simObj$mxARCL_FIML$data$observed <- full_raw_data
    simObj$mxARCL_cov$data$observed <- cov(full_raw_data)
    ##### build models with openMx:#####

    # covariance based:

    full_fitMxARCL_cov <- tryCatch(mxTryHard(simObj$mxARCL_cov),
                                   warning = function(w){
                                     print("warning: did not find proper solution")
                                     return(NA)
                                   },
                                   error = function(e){
                                     print("warning: did not find proper solution")
                                     return(NA)
                                   }
    )
    if(is.logical(full_fitMxARCL_cov)){
      error = error+ 1
      next
    }

    # FIML based:

    full_fitMxARCL_FIML <- tryCatch(mxTryHard(simObj$mxARCL_FIML),
                                    warning = function(w){
                                      print("warning: did not find proper solution")
                                      return(NA)
                                    },
                                    error = function(e){
                                      print("warning: did not find proper solution")
                                      return(NA)
                                    }
    )
    if(is.logical(full_fitMxARCL_FIML)){
      error = error+ 1
      next
    }

    ## with CV

    cov_train_data = tryCatch(mxData(observed = cov(train_raw_data), numObs = length(Folds$Fold1), type = "cov"),
                              warning = function(w){
                                print("warning: did not find proper solution")
                                return(NA)
                              },
                              error = function(e){
                                print("warning: did not find proper solution")
                                return(NA)
                              }
    )
    cov_test_data = tryCatch(mxData(cov(test_raw_data), numObs = length(Folds$Fold2), type = "cov"),
                             warning = function(w){
                               print("warning: did not find proper solution")
                               return(NA)
                             },
                             error = function(e){
                               print("warning: did not find proper solution")
                               return(NA)
                             }
    )

    # run Models:
    if(!is.logical(cov_train_data)&&!is.logical(cov_test_data)){
      CV_covMxModel <- simObj$mxARCL_cov
      CV_covMxModel$data <- cov_train_data
      CV_covMxModel <- tryCatch(mxTryHard(CV_covMxModel, silent = T),
                                warning = function(w){
                                  print("warning: did not find proper solution")
                                  return(NA)
                                },
                                error = function(e){
                                  print("warning: did not find proper solution")
                                  return(NA)
                                }
      )
    }else{
      error = error+ 1
      next
    }
    if(is.logical(CV_covMxModel)){
      error = error+ 1
      next
    }

    # FIML based:
    # create mxData objects:
    FIML_train_data = tryCatch(mxData(train_raw_data,type = "raw"),
                               warning = function(w){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               },
                               error = function(e){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               }
    )

    FIML_test_data = tryCatch(mxData(test_raw_data,  type = "raw"),
                              warning = function(w){
                                print("warning: did not find proper solution")
                                return(NA)
                              },
                              error = function(e){
                                print("warning: did not find proper solution")
                                return(NA)
                              }
    )

    # run Models:
    if(!is.logical(FIML_train_data)&& !is.logical(FIML_test_data)){
      FIML_CV_rawMxModel <- simObj$mxARCL_FIML
      FIML_CV_rawMxModel$data <- FIML_train_data
      FIML_CV_rawMxModel <- tryCatch(mxTryHard(FIML_CV_rawMxModel, silent = T),
                                     warning = function(w){
                                       print("warning: did not find proper solution")
                                       return(NA)
                                     },
                                     error = function(e){
                                       print("warning: did not find proper solution")
                                       return(NA)
                                     }
      )

    }else{
      error = error+ 1
      next
    }

    ## check the model fitting for errors #####
    if(is.logical(FIML_CV_rawMxModel)){
      error = error+ 1
      next
    }


    # save base parameters
    if(any(!names(parameterValues$Base_Model_cov[,iteration]) == colnames(getUniqueA(full_fitMxARCL_cov$A)))){
      print("Warning: dimension names do not match! Could not save parameters")
    }else{
      parameterValues$Base_Model_cov[,iteration] <- getUniqueA(full_fitMxARCL_cov$A)
      parameterValues$Base_Model_FIML[,iteration] <- getUniqueA(full_fitMxARCL_FIML$A)

      AmatEst = matrix(getUniqueA(full_fitMxARCL_cov$A), 5, byrow = F)
      StempMat <- getUniqueA(full_fitMxARCL_cov$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$Base_Model_cov[,iteration] <- computeRMSE(A_est = AmatEst,
                                                     S_est = SmatEst,
                                                     Apop = Avalues,
                                                     Spop = Svalues,
                                                     Afree = Afree,
                                                     Sfree = Sfree)

      AmatEst = matrix(getUniqueA(full_fitMxARCL_FIML$A), 5, byrow = F)
      StempMat <- getUniqueA(full_fitMxARCL_FIML$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$Base_Model_FIML[,iteration] <- computeRMSE(A_est = AmatEst,
                                                      S_est = SmatEst,
                                                      Apop = Avalues,
                                                      Spop = Svalues,
                                                      Afree = Afree,
                                                      Sfree = Sfree)

    }
    ############################################BIS HIER GEKOMMEN ###########################
    ###### computation without CV ######
    # OpenMx
    ##### Covariance Based Models ####
    ## start fitting with different penalty values:

    full_cov_reg_Model <- tryCatch(fitRegModels(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",
                                                fitfun = "FML",pen_on = "A",selectedA = penAFull,pen_start = 0,pen_end = .4,pen_stepsize = .01,fit_index = "BIC"),

                                   error = function(e){
                                     print("warning: did not find proper solution")
                                     return(NA)
                                   }
    )
    if(!is.logical(full_cov_reg_Model)){


      ## search best penalty values:
      valid_fitmeasures <- (full_cov_reg_Model$`fit measures`[,"convergence"]==0) == (full_cov_reg_Model$`fit measures`[,"negative variances"]==0)
      valid_fitmeasures <- full_cov_reg_Model$`fit measures`[valid_fitmeasures,]

      cov_minimum_AIC <- valid_fitmeasures[which(valid_fitmeasures[,"AIC"]==min(valid_fitmeasures[,"AIC"])), "penalty"] # get best penalty value

      cov_minimum_BIC <- valid_fitmeasures[which(valid_fitmeasures[,"BIC"]==min(valid_fitmeasures[,"BIC"])), "penalty"]

      ## get parameter values:

      full_cov_reg_Model_AIC <- createRegModel(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = cov_minimum_AIC,pen_on = "A",selectedA = penAFull)
      full_cov_reg_Model_AIC <- mxRun(full_cov_reg_Model_AIC, silent = T)


      full_cov_reg_Model_BIC <- createRegModel(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = cov_minimum_BIC,pen_on = "A",selectedA = penAFull)
      full_cov_reg_Model_BIC <- mxRun(full_cov_reg_Model_BIC, silent = T)

      #### saving parameters
      parameterValues$AIC[,iteration] <- getUniqueA(full_cov_reg_Model_AIC$BaseModel$A)
      parameterValues$BIC[,iteration] <- getUniqueA(full_cov_reg_Model_BIC$BaseModel$A)

      AmatEst = matrix(getUniqueA(full_cov_reg_Model_AIC$BaseModel$A), 5, byrow = F)
      StempMat <- getUniqueA(full_cov_reg_Model_AIC$BaseModel$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$AIC[,iteration] <- computeRMSE(A_est = AmatEst,
                                          S_est = SmatEst,
                                          Apop = Avalues,
                                          Spop = Svalues,
                                          Afree = Afree,
                                          Sfree = Sfree)

      AmatEst = matrix(getUniqueA(full_cov_reg_Model_BIC$BaseModel$A), 5, byrow = F)
      StempMat <- getUniqueA(full_cov_reg_Model_BIC$BaseModel$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$BIC[,iteration] <- computeRMSE(A_est = AmatEst,
                                          S_est = SmatEst,
                                          Apop = Avalues,
                                          Spop = Svalues,
                                          Afree = Afree,
                                          Sfree = Sfree)


      print("Covariance without CV successful")
    }else{
      error = error+ 1
      next
    }

    # with FIML

    ##### FIML #####

    ## start fitting with different penalty values:

    full_FIML_reg_Model <- tryCatch(fitRegModels(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel",
                                                 fitfun = "FIML",pen_on = "A",selectedA = penAFull,pen_start = 0,pen_end = .4,pen_stepsize = .01,fit_index = "BIC"),
                                    error = function(e){
                                      print("warning: did not find proper solution")
                                      return(NA)
                                    }
    )
    if(!is.logical(full_FIML_reg_Model)){



      ## search best penalty values:
      valid_fitmeasures <- (full_FIML_reg_Model$`fit measures`[,"convergence"]==0) == (full_FIML_reg_Model$`fit measures`[,"negative variances"]==0)
      valid_fitmeasures <- full_FIML_reg_Model$`fit measures`[valid_fitmeasures,]

      FIML_minimum_AIC <- valid_fitmeasures[which(valid_fitmeasures[,"AIC"]==min(valid_fitmeasures[,"AIC"])), "penalty"] # get best penalty value

      FIML_minimum_BIC <- valid_fitmeasures[which(valid_fitmeasures[,"BIC"]==min(valid_fitmeasures[,"BIC"])), "penalty"]

      # get parameters:
      full_FIML_reg_Model_AIC <- createRegModel(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel", fitfun = "FIML",pen_value = FIML_minimum_AIC,pen_on = "A",selectedA = penAFull)
      full_FIML_reg_Model_AIC <- mxRun(full_FIML_reg_Model_AIC, silent = T)

      full_FIML_reg_Model_BIC <- createRegModel(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel",fitfun = "FIML",pen_value = FIML_minimum_BIC,pen_on = "A",selectedA = penAFull)
      full_FIML_reg_Model_BIC <- mxRun(full_FIML_reg_Model_BIC, silent = T)


      # save parameters

      parameterValues$FIML_AIC[,iteration] <- getUniqueA(full_FIML_reg_Model_AIC$BaseModel$A)
      parameterValues$FIML_BIC[,iteration] <- getUniqueA(full_FIML_reg_Model_BIC$BaseModel$A)

      AmatEst = matrix(getUniqueA(full_FIML_reg_Model_AIC$BaseModel$A), 5, byrow = F)
      StempMat <- getUniqueA(full_FIML_reg_Model_AIC$BaseModel$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$FIML_AIC[,iteration] <- computeRMSE(A_est = AmatEst,
                                               S_est = SmatEst,
                                               Apop = Avalues,
                                               Spop = Svalues,
                                               Afree = Afree,
                                               Sfree = Sfree)

      AmatEst = matrix(getUniqueA(full_FIML_reg_Model_BIC$BaseModel$A), 5, byrow = F)
      StempMat <- getUniqueA(full_FIML_reg_Model_BIC$BaseModel$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$FIML_BIC[,iteration] <- computeRMSE(A_est = AmatEst,
                                               S_est = SmatEst,
                                               Apop = Avalues,
                                               Spop = Svalues,
                                               Afree = Afree,
                                               Sfree = Sfree)


      print("FIML without CV successful")
    }else{
      error = error+ 1
      next
    }



    ##### with CV ####

    # OpenMx
    ## Note: Cross validation is performed differently from Jacobucci (2016): it is assumed that the researcher only
    # has one sample of site sample_size and this sample is split into a training and a testing sample


    ### covariance based: ###

    CV_cov_reg_Model <- tryCatch(fitRegModels(model = CV_covMxModel, data_type = "cov",model_type = "mxModel",
                                              fitfun = "FML",pen_on = "A",selectedA = penAFull,pen_start = 0,pen_end = .4,pen_stepsize = .01,
                                              fit_index = "BIC", CV = T, Test_Sample = cov_test_data) ,

                                 error = function(e){
                                   print("warning: did not find proper solution")
                                   return(NA)
                                 }
    )
    if(!is.logical(CV_cov_reg_Model)){

      # get best penalty values:
      valid_fitmeasures <- (CV_cov_reg_Model$`fit measures`[,"convergence"]==0) == (CV_cov_reg_Model$`fit measures`[,"negative variances"]==0)
      valid_fitmeasures <- CV_cov_reg_Model$`fit measures`[valid_fitmeasures,]


      minimum_CV_cov_m2LL <- valid_fitmeasures[which(valid_fitmeasures[,"CV m2LL"]==min(valid_fitmeasures[,"CV m2LL"])),"penalty"]

      # getting parameters:
      cov_reg_Model_CV_m2LL <- createRegModel(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = minimum_CV_cov_m2LL,pen_on = "A",selectedA = penAFull)
      cov_reg_Model_CV_m2LL <- mxRun(cov_reg_Model_CV_m2LL, silent = T)


      # save parameters

      parameterValues$CV_m2LL[,iteration] <- getUniqueA(cov_reg_Model_CV_m2LL$BaseModel$A)

      AmatEst = matrix(getUniqueA(cov_reg_Model_CV_m2LL$BaseModel$A), 5, byrow = F)
      StempMat <- getUniqueA(cov_reg_Model_CV_m2LL$BaseModel$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$CV_m2LL[,iteration] <- computeRMSE(A_est = AmatEst,
                                              S_est = SmatEst,
                                              Apop = Avalues,
                                              Spop = Svalues,
                                              Afree = Afree,
                                              Sfree = Sfree)

      print("Covariance with CV successful")
    }else{
      error = error+ 1
      next
    }

    ## FIML based models ###

    FIML_CV_reg_Model <- tryCatch(fitRegModels(model = FIML_CV_rawMxModel, data_type = "raw",model_type = "mxModel",
                                               fitfun = "FIML",pen_on = "A",selectedA = penAFull,pen_start = 0,pen_end = .4,pen_stepsize = .01,
                                               fit_index = "BIC", CV = T, Test_Sample = FIML_test_data),


                                  error = function(e){
                                    print("warning: did not find proper solution")
                                    return(NA)
                                  }
    )

    if(!is.logical(FIML_CV_reg_Model)){

      FIML_minimum_CV_m2LL <- FIML_CV_reg_Model$`fit measures`[which(FIML_CV_reg_Model$`fit measures`[,"CV m2LL"]==min(FIML_CV_reg_Model$`fit measures`[,"CV m2LL"])),"penalty"]

      FIML_reg_Model_CV_m2LL <- createRegModel(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel",fitfun = "FIML",pen_value = FIML_minimum_CV_m2LL,pen_on = "A",selectedA = penAFull)
      FIML_reg_Model_CV_m2LL <- mxRun(FIML_reg_Model_CV_m2LL, silent = T)

      # save parameters

      parameterValues$FIML_CV_m2LL[,iteration] <- getUniqueA(FIML_reg_Model_CV_m2LL$BaseModel$A)

      AmatEst = matrix(getUniqueA(FIML_reg_Model_CV_m2LL$BaseModel$A), 5, byrow = F)
      StempMat <- getUniqueA(FIML_reg_Model_CV_m2LL$BaseModel$S)
      SmatEst = diag(StempMat[16:20])
      RMSE$FIML_CV_m2LL[,iteration] <- computeRMSE(A_est = AmatEst,
                                                   S_est = SmatEst,
                                                   Apop = Avalues,
                                                   Spop = Svalues,
                                                   Afree = Afree,
                                                   Sfree = Sfree)

      print("FIML with CV successful")
    }else{
      error = error+ 1
      next
    }
    # set progress bar:
    setTkProgressBar(global_pb, iteration, label=paste( round(iteration/total_repetitions*100, 0),
                                                        "% done"))

    iteration = iteration+1


  }

  # overall evaluataion:

  for(i in 1:total_repetitions){
    # a12 is non-zero in true model:
    overall_evaluation$AIC_false_zero[i] <- sum(abs(parameterValues$AIC[Alabel[Avalues > .001],i]) < .001)
    overall_evaluation$FIML_AIC_false_zero[i] <- sum(abs(parameterValues$FIML_AIC[Alabel[Avalues > .001],i]) < .001)
    overall_evaluation$BIC_false_zero[i] <- sum(abs(parameterValues$BIC[Alabel[Avalues > .001],i]) < .001)
    overall_evaluation$FIML_BIC_false_zero[i] <- sum(abs(parameterValues$FIML_BIC[Alabel[Avalues > .001],i]) < .001)
    overall_evaluation$CV_m2LL_false_zero[i] <- sum(abs(parameterValues$CV_m2LL[Alabel[Avalues > .001],i]) < .001)
    overall_evaluation$FIML_CV_m2LL_false_zero[i] <- sum(abs(parameterValues$FIML_CV_m2LL[Alabel[Avalues > .001],i]) < .001)

    # a21 is zero in true model

    overall_evaluation$AIC_false_nonzero[i] <- sum(abs(parameterValues$AIC[Alabel[Avalues < .001],i]) > .001)
    overall_evaluation$FIML_AIC_false_nonzero[i] <- sum(abs(parameterValues$FIML_AIC[Alabel[Avalues < .001],i]) > .001)
    overall_evaluation$BIC_false_nonzero[i] <- sum(abs(parameterValues$BIC[Alabel[Avalues < .001],i]) > .001)
    overall_evaluation$FIML_BIC_false_nonzero[i] <- sum(abs(parameterValues$FIML_BIC[Alabel[Avalues < .001],i]) > .001)
    overall_evaluation$CV_m2LL_false_nonzero[i] <- sum(abs(parameterValues$CV_m2LL[Alabel[Avalues < .001],i]) > .001)
    overall_evaluation$FIML_CV_m2LL_false_nonzero[i] <- sum(abs(parameterValues$FIML_CV_m2LL[Alabel[Avalues < .001],i]) > .001)

  }


  sumsOverallEval <- apply(overall_evaluation,2, sum, na.rm = TRUE)
  names(sumsOverallEval) <- colnames(overall_evaluation)
  perRunOverAll <- sumsOverallEval/total_repetitions


  # note: regularization makes fit worse. The reason is that both parameters, the non-zero and the true-zero one are regularized. The non-zero one
  # gets pulled away from its true value and this impacts the over-all RMSEA more than setting a parameter that is close to zero to zero
  save(overall_evaluation, parameterValues, RMSE, total_repetitions, Avalues, penAelement, Svalues, error, seed, crossEffect, autoEffect, file = paste("Simulation3_N", sampleSize, "_1_1000.RData", sep = ""))

}
