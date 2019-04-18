#' Simulation study 1
#'
#'
#'
#' @param sampleSize sampleSize
#' @param seed seed
#' @param wd working directory. The results will be saved here
#'
#' @author Jannik Orzek
#' @import OpenMx caret laremm regsem tcltk
#'
#' @examples
#'
#' @export
#'
#'

simulation1 <- function(sampleSize, seed, wd){
  setwd(wd)
  set.seed(seed)
  # in this model, each factor has one item that it doesn't share with the other items; this is done to
  # keep the items linked to one factor alone and to not let them switch between them
  # Model specification:

  # variables
  LV <- c("f1","f2", "f3")
  manif <- c("x1","x2","x3","x4","x5","x6", "x7","x8","x9","x10","x11", "x12")

  # Amat
  A_values <- matrix(0, nrow = length(c(manif,LV)), ncol = length(c(manif,LV)))
  rownames(A_values) <- c(manif,LV)
  colnames(A_values) <- c(manif,LV)

  A_values["x1", "f1"] <- .9
  A_values["x2", "f1"] <- .9
  A_values["x3", "f1"] <- .8
  A_values["x4", "f1"] <- .2

  A_values["x5", "f2"] <- .9
  A_values["x6", "f2"] <- .9
  A_values["x7", "f2"] <- .8
  A_values["x8", "f2"] <- .2

  A_values["x9", "f3"]  <- .9
  A_values["x10", "f3"] <- .9
  A_values["x11", "f3"] <- .8
  A_values["x12", "f3"] <- .2

  A_free <- matrix(FALSE, nrow = length(c(manif,LV)), ncol = length(c(manif,LV)))
  rownames(A_free) <- c(manif,LV)
  colnames(A_free) <- c(manif,LV)

  A_free[A_values > 0] <- TRUE

  # S matrix
  S_values <- matrix(0, nrow = length(c(manif,LV)), ncol = length(c(manif,LV)))
  rownames(S_values) <- c(manif,LV)
  colnames(S_values) <- c(manif,LV)

  S_values["f1", "f1"] <- 1
  S_values["x1", "x1"] <- .19
  S_values["x2", "x2"] <- .19
  S_values["x3", "x3"] <- .36
  S_values["x4", "x4"] <- .96

  S_values["f2", "f2"] <- 1
  S_values["x5", "x5"] <- .19
  S_values["x6", "x6"] <- .19
  S_values["x7", "x7"] <- .36
  S_values["x8", "x8"] <- .96

  S_values["f3", "f3"] <- 1
  S_values["x9", "x9"] <- .19
  S_values["x10", "x10"] <- .19
  S_values["x11", "x11"] <- .36
  S_values["x12", "x12"] <- .96

  S_free <- matrix(FALSE, nrow = length(c(manif,LV)), ncol = length(c(manif,LV)))
  rownames(S_free) <- c(manif,LV)
  colnames(S_free) <- c(manif,LV)
  S_free[S_values >0] <- TRUE

  # F matrix
  F_values <- matrix(0, nrow = length(c(manif)), ncol = length(c(manif,LV)))
  diag(F_values[1:length(c(manif)), 1:length(c(manif))]) <- 1

  F_free <- matrix(FALSE, nrow = length(c(manif)), ncol = length(c(manif,LV)))

  # build simulation model
  SimModel <- mxModel(model = "SimModel", type="RAM",
                      manifestVars = manif, latentVars = LV,
                      mxMatrix(type = "Full", nrow = length(c(manif,LV)), ncol = length(c(manif,LV)), free = A_free, values = A_values, name = "A"),
                      mxMatrix(type = "Full", nrow = length(c(manif,LV)), ncol = length(c(manif,LV)), free = S_free, values = S_values, name = "S"),
                      mxMatrix(type = "Full", nrow = length(c(manif)), ncol = length(c(manif,LV)), free = F_free, values = F_values, name = "F"),
                      mxMatrix(type = "Full", nrow = 1, ncol = length(c(manif,LV)), free= F, values = 0, name = "M")
  )



  sampleSize = sampleSize


  #### Model for Analysis ####

  ## in lavaan
  analys_Model <- "
f1 =~ NA*x1 + l12*x2 + l13*x3 + l14*x4          + l16*x6 + l17*x7 + l18*x8
f2 =~ NA*x5 + l26*x6 + l27*x7 + l28*x8          + l210*x10 + l211*x11 + l212*x12
f3 =~ NA*x9 + l310*x10 + l311*x11 + l312*x12    + l32*x2 + l33*x3 + l34*x4
f1 ~~1*f1
f2 ~~1*f2
f3 ~~1*f3
f1 ~~ 0*f3+0*f2
f2 ~~ 0*f3
"
  ## in OpenMx

  # loadings
  f1loadings <- mxPath(from = "f1", to = c("x1","x2","x3","x4","x6","x7","x8"), free = T)
  f2loadings <- mxPath(from = "f2", to = c("x5","x6","x7","x8","x10","x11","x12"), free = T)
  f3loadings <- mxPath(from = "f3", to = c("x9","x10","x11", "x12","x2","x3","x4"), free = T)
  ######PROBLEM: Loadings wecheseln ab und zu: manchmal werden Items 1-3 Faktor 3, Items 4-6 Faktor 1 und Items 7-9 Faktor 2
  # covariances
  cov_manif <- mxPath(from = manif, arrows = 2,values = 1, free = T)
  cov_LV <- mxPath(from = LV, arrows = 2, free =c(F), values = c(1),connect = "single")


  ##### start outer loop #####
  # settings
  total_repetitions = 1000
  iteration = 1
  improper_solutions = 0

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
                                   "FIML_CV_m2LL_false_zero"=rep(NA,total_repetitions),


                                   "regsem_BIC_false_nonzero"=rep(NA,total_repetitions),
                                   "regsem_BIC_false_zero"=rep(NA,total_repetitions),
                                   "regsem_AIC_false_nonzero"=rep(NA,total_repetitions),
                                   "regsem_AIC_false_zero"=rep(NA,total_repetitions),
                                   "regsem_CV_chisq_false_nonzero" = rep(NA,total_repetitions),
                                   "regsem_CV_chisq_false_zero" = rep(NA,total_repetitions)
  )


  # regularization:
  selectedA = "all"

  # matrix mit logicals for all free A parameters. Is used for selecting the estimated parameters
  A_estimated <- A_free
  A_estimated["x6", "f1"] <- TRUE
  A_estimated["x7", "f1"] <- TRUE
  A_estimated["x8", "f1"] <- TRUE
  A_estimated["x10", "f2"] <- TRUE
  A_estimated["x11", "f2"] <- TRUE
  A_estimated["x12", "f2"] <- TRUE
  A_estimated["x2", "f3"] <- TRUE
  A_estimated["x3", "f3"] <- TRUE
  A_estimated["x4", "f3"] <- TRUE

  S_estimated <- S_free
  S_estimated["f1", "f1"] <- FALSE
  S_estimated["f2", "f2"] <- FALSE
  S_estimated["f3", "f3"] <- FALSE


  parameter_Table <- matrix(NA, nrow = 21, ncol = total_repetitions)
  rownames(parameter_Table) <- c("f1->x1","f1->x2","f1->x3","f1->x4","f1->x6","f1->x7","f1->x8",
                                 "f2->x5","f2->x6","f2->x7","f2->x8","f2->x10","f2->x11","f2->x12",
                                 "f3->x2","f3->x3","f3->x4","f3->x9","f3->x10","f3->x11","f3->x12")

  parameterValues <- list("Base_Model_cov" = parameter_Table,
                          "Base_Model_FIML" = parameter_Table,
                          "AIC"=parameter_Table, "BIC"=parameter_Table,
                          "FIML_BIC"=parameter_Table, "FIML_AIC"=parameter_Table,
                          "CV_m2LL" = parameter_Table,
                          "FIML_CV_m2LL" =parameter_Table,
                          "regsem_BIC" =parameter_Table, "regsem_AIC" =parameter_Table,
                          "CV_regsem_chisqu" =parameter_Table)

  RMSE_table <- matrix(NA, nrow = 1, ncol = total_repetitions)
  rownames(RMSE_table) <- c("RMSE")
  RMSE <- list("Base_Model_cov" = RMSE_table,
               "Base_Model_FIML" = RMSE_table,
               "AIC"=RMSE_table, "BIC"=RMSE_table,
               "FIML_AIC"=RMSE_table, "FIML_BIC"=RMSE_table,
               "CV_m2LL" = RMSE_table,
               "FIML_CV_m2LL" =RMSE_table,
               "regsem_AIC" =RMSE_table, "regsem_BIC" =RMSE_table,
               "CV_regsem_chisqu" =RMSE_table)
  # start
  global_pb <- tkProgressBar(title = "progress bar", min = 0,
                             max = total_repetitions, width = 300) # progress - bar

  error = 0 # counts number of errors in simulation

  while(iteration <= total_repetitions){

    ### simulate dataset ####

    raw_data <- mxGenerateData(SimModel, sampleSize)
    colnames(raw_data) <- c("x1","x2","x3","x4","x5","x6", "x7","x8","x9","x10","x11", "x12")

    raw_data_base <- scale(raw_data)

    cov_data_base <- cov(raw_data_base)

    # split dataset in train and test dataset
    Folds <- createFolds(c(1:sampleSize),2)

    train_raw_data <- raw_data[Folds$Fold1,]
    test_raw_data <- raw_data[Folds$Fold2,]

    train_raw_data <- scale(train_raw_data)
    test_raw_data <- scale(test_raw_data)

    cov_train_data <- cov(train_raw_data)
    cov_test_data <- cov(test_raw_data)

    ## build models ##
    # lavaan:
    # full model
    lavFullModel <- tryCatch(cfa(analys_Model, sample.cov = cov_data_base, sample.nobs = sampleSize),
                             warning = function(w){
                               print("warning: did not find proper solution")
                               return(NA)
                             },
                             error = function(e){
                               print("warning: did not find proper solution")
                               return(NA)
                             }
    )
    if(is.logical(lavFullModel)){
      error = error+ 1
      next
    }
    #train model
    lavTrainModel <- tryCatch(cfa(analys_Model, sample.cov = cov_train_data, sample.nobs = length(Folds$Fold1)),
                              warning = function(w){
                                print("warning: did not find proper solution")
                                return(NA)
                              },
                              error = function(e){
                                print("warning: did not find proper solution")
                                return(NA)
                              }
    )

    if(is.logical(lavTrainModel)){
      error = error+ 1
      next
    }
    # OpenMx:
    # covariance based:
    mxCov_Data <- tryCatch(mxData(observed = cov_data_base, type = "cov",numObs = sampleSize),
                           warning = function(w){
                             print("warning: did not find proper solution")
                             return(NA)
                           },
                           error = function(e){
                             print("warning: did not find proper solution")
                             return(NA)
                           }
    )

    if(!is.logical(mxCov_Data)){
      covFullMxModel <- mxModel(name = "mxModel", type = "RAM", manifestVars = manif, latentVars = LV, data =
                                  mxCov_Data, f1loadings,f2loadings,f3loadings, cov_manif, cov_LV)
      covFullMxModel <- tryCatch(mxRun(covFullMxModel, silent = T),
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
    if(is.logical(covFullMxModel)){
      error = error+ 1
      next
    }
    #FIML based:

    raw_Data <- tryCatch(mxData(observed = raw_data_base, type = "raw"),
                         warning = function(w){
                           print("warning: did not find proper solution")
                           return(NA)
                         },
                         error = function(e){
                           print("warning: did not find proper solution")
                           return(NA)
                         }
    )
    if(!is.logical(raw_Data)){
      rawMxModel <- mxModel(name = "mxModel", type = "RAM", manifestVars = manif, latentVars = LV, data =
                              raw_Data, f1loadings,f2loadings,f3loadings, cov_manif, cov_LV,
                            mxPath(from = "one", to = manif, values = 0, free = F))
      rawMxModel <- tryCatch(mxRun(rawMxModel, silent = T),
                             warning = function(w){
                               print("warning: did not find proper solution")
                               return(NA)
                             },
                             error = function(e){
                               print("warning: did not find proper solution")
                               return(NA)
                             }
      )

      #rawMxModel_sat <- mxRefModels(rawMxModel, run = TRUE)$Saturated
    }else{
      error = error+ 1
      next
    }
    if(is.logical(rawMxModel)){
      error = error+ 1
      next
    }

    # CV based:
    # cov:

    mx_cov_train_data = tryCatch(mxData(observed = cov_train_data, numObs = length(Folds$Fold1), type = "cov"),
                                 warning = function(w){
                                   print("warning: did not find proper solution")
                                   return(NA)
                                 },
                                 error = function(e){
                                   print("warning: did not find proper solution")
                                   return(NA)
                                 }
    )
    mx_cov_test_data = tryCatch(mxData(cov_test_data, numObs = length(Folds$Fold2), type = "cov"),
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
    if(!is.logical(mx_cov_train_data)&&!is.logical(mx_cov_test_data)){
      CV_covMxModel <- mxModel(name = "mxModel", type = "RAM", manifestVars = manif, latentVars = LV, data =
                                 mx_cov_train_data, f1loadings,f2loadings,f3loadings, cov_manif, cov_LV)
      CV_covMxModel <- tryCatch(mxRun(CV_covMxModel, silent = T),
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
      FIML_CV_rawMxModel <- mxModel(name = "mxModel", type = "RAM", manifestVars = manif, latentVars = LV, data =
                                      FIML_train_data, f1loadings,f2loadings,f3loadings, cov_manif, cov_LV,
                                    mxPath(from = 'one', to = manif, values = 0,free = FALSE))
      FIML_CV_rawMxModel <- tryCatch(mxRun(FIML_CV_rawMxModel, silent = T),
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
    parameterValues$Base_Model_cov[,iteration] <- cvectorize(covFullMxModel$A$values[A_estimated])
    RMSE$Base_Model_cov[,iteration] <- computeRMSE(A_est = covFullMxModel$A$values, S_est = covFullMxModel$S$values,
                                                   Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]

    parameterValues$Base_Model_FIML[,iteration] <- cvectorize(rawMxModel$A$values[A_estimated])
    RMSE$Base_Model_FIML[,iteration] <- computeRMSE(A_est = rawMxModel$A$values, S_est = rawMxModel$S$values,
                                                    Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]
    ###### computation without CV ######
    # OpenMx
    ##### Covariance Based Models ####
    ## start fitting with different penalty values:

    cov_reg_Model <- tryCatch(fitRegModels(covFullMxModel, data_type = "cov",model_type = "mxModel",
                                           fitfun = "FML",pen_on = "A",selectedA = selectedA,pen_start = 0,pen_end = .4,pen_stepsize = .01,fit_index = "BIC"),

                              error = function(e){
                                print("warning: did not find proper solution")
                                return(NA)
                              }
    )
    if(!is.logical(cov_reg_Model)){


      ## search best penalty values:
      valid_fitmeasures <- (cov_reg_Model$`fit measures`[,"convergence"]==0) == (cov_reg_Model$`fit measures`[,"negative variances"]==0)
      valid_fitmeasures <- cov_reg_Model$`fit measures`[valid_fitmeasures,]

      cov_minimum_AIC <- valid_fitmeasures[which(valid_fitmeasures[,"AIC"]==min(valid_fitmeasures[,"AIC"])), "penalty"] # get best penalty value

      cov_minimum_BIC <- valid_fitmeasures[which(valid_fitmeasures[,"BIC"]==min(valid_fitmeasures[,"BIC"])), "penalty"]

      ## get parameter values:

      cov_reg_Model_AIC <- createRegModel(covFullMxModel, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = cov_minimum_AIC,pen_on = "A",selectedA = selectedA)
      cov_reg_Model_AIC <- mxRun(cov_reg_Model_AIC, silent = T)


      cov_reg_Model_BIC <- createRegModel(covFullMxModel, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = cov_minimum_BIC,pen_on = "A",selectedA = selectedA)
      cov_reg_Model_BIC <- mxRun(cov_reg_Model_BIC, silent = T)

      #### saving parameters
      parameterValues$AIC[,iteration] <- cvectorize(cov_reg_Model_AIC$BaseModel$A$values[A_estimated])
      parameterValues$BIC[,iteration] <- cvectorize(cov_reg_Model_BIC$BaseModel$A$values[A_estimated])

      # compute RMSE
      RMSE$AIC[,iteration] <- computeRMSE(A_est = cov_reg_Model_AIC$BaseModel$A$values, S_est = cov_reg_Model_AIC$BaseModel$S$values,
                                          Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]

      RMSE$BIC[,iteration] <- computeRMSE(A_est = cov_reg_Model_BIC$BaseModel$A$values, S_est = cov_reg_Model_BIC$BaseModel$S$values,
                                          Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]

      print("Covariance without CV successful")
    }else{
      error = error+ 1
      next
    }

    # with FIML

    ##### FIML #####

    ## start fitting with different penalty values:

    FIML_reg_Model <- tryCatch(fitRegModels(model = rawMxModel, data_type = "raw",model_type = "mxModel",
                                            fitfun = "FIML",pen_on = "A",selectedA = selectedA,pen_start = 0,pen_end = .4,pen_stepsize = .01,fit_index = "BIC"),
                               error = function(e){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               }
    )
    if(!is.logical(FIML_reg_Model)){



      ## search best penalty values:
      valid_fitmeasures <- (FIML_reg_Model$`fit measures`[,"convergence"]==0) == (FIML_reg_Model$`fit measures`[,"negative variances"]==0)
      valid_fitmeasures <- FIML_reg_Model$`fit measures`[valid_fitmeasures,]

      FIML_minimum_AIC <- valid_fitmeasures[which(valid_fitmeasures[,"AIC"]==min(valid_fitmeasures[,"AIC"])), "penalty"] # get best penalty value

      FIML_minimum_BIC <- valid_fitmeasures[which(valid_fitmeasures[,"BIC"]==min(valid_fitmeasures[,"BIC"])), "penalty"]

      # get parameters:
      FIML_reg_Model_AIC <- createRegModel(model = rawMxModel, data_type = "raw",model_type = "mxModel", fitfun = "FIML",pen_value = FIML_minimum_AIC,pen_on = "A",selectedA = selectedA)
      FIML_reg_Model_AIC <- mxRun(FIML_reg_Model_AIC, silent = T)

      FIML_reg_Model_BIC <- createRegModel(model = rawMxModel, data_type = "raw",model_type = "mxModel",fitfun = "FIML",pen_value = FIML_minimum_BIC,pen_on = "A",selectedA = selectedA)
      FIML_reg_Model_BIC <- mxRun(FIML_reg_Model_BIC, silent = T)


      # save parameters

      parameterValues$FIML_AIC[,iteration] <- cvectorize(FIML_reg_Model_AIC$BaseModel$A$values[A_estimated])
      parameterValues$FIML_BIC[,iteration] <- cvectorize(FIML_reg_Model_BIC$BaseModel$A$values[A_estimated])

      # compute RMSE
      RMSE$FIML_AIC[,iteration] <- computeRMSE(A_est = FIML_reg_Model_AIC$BaseModel$A$values, S_est = FIML_reg_Model_AIC$BaseModel$S$values,
                                               Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]

      RMSE$FIML_BIC[,iteration] <- computeRMSE(A_est = FIML_reg_Model_BIC$BaseModel$A$values, S_est = FIML_reg_Model_BIC$BaseModel$S$values,
                                               Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]


      print("FIML without CV successful")
    }else{
      error = error+ 1
      next
    }


    # regsem:

    regsem_out <- tryCatch(cv_regsem(lavFullModel, n.lambda = 40, pars_pen = "loadings", fit.ret=c("BIC", "AIC"), metric = "BIC"),
                           warning = function(w){
                             print("warning: did not find proper solution")
                             return(NA)
                           },
                           error = function(e){
                             print("warning: did not find proper solution")
                             return(NA)
                           }
    )
    if(!is.logical(regsem_out)){

      print("regsem without CV successful")

      best_lambda_aic <- regsem_out$fits[which(regsem_out$fits[,"AIC"]== min(regsem_out$fits[,"AIC"], na.rm = T)),"lambda"]
      best_lambda_bic <- regsem_out$fits[which(regsem_out$fits[,"BIC"]== min(regsem_out$fits[,"BIC"], na.rm = T)),"lambda"]

      # save parameters:
      # BIC
      bic_regsem <- regsem(lavFullModel, lambda = best_lambda_bic, pars_pen = "loadings")
      bic_regsem_param <- bic_regsem$out$pars
      parameterValues$regsem_BIC[,iteration] <- as.numeric(bic_regsem_param[c("f1 -> x1","f1 -> x2","f1 -> x3","f1 -> x4","f1 -> x6","f1 -> x7","f1 -> x8","f2 -> x5","f2 -> x6","f2 -> x7","f2 -> x8","f2 -> x10", "f2 -> x11", "f2 -> x12", "f3 -> x2","f3 -> x3","f3 -> x4","f3 -> x9","f3 -> x10","f3 -> x11","f3 -> x12")])

      # compute RMSE
      RMSE$regsem_BIC[,iteration] <- sqrt(
        (sum(
          (as.numeric(bic_regsem_param[c("f1 -> x1","f1 -> x2","f1 -> x3","f1 -> x4","f1 -> x6","f1 -> x7","f1 -> x8","f2 -> x5","f2 -> x6","f2 -> x7","f2 -> x8","f2 -> x10", "f2 -> x11", "f2 -> x12", "f3 -> x2","f3 -> x3","f3 -> x4","f3 -> x9","f3 -> x10","f3 -> x11","f3 -> x12")])
           - c(A_values[c(1:4, 6,7,8),13], A_values[c(5:8, 10,11,12),14], A_values[c(2,3,4, 9:12),15]))^2
        ) +
          sum(
            (as.numeric(bic_regsem_param[c("x1 ~~ x1","x2 ~~ x2","x3 ~~ x3","x4 ~~ x4","x5 ~~ x5","x6 ~~ x6","x7 ~~ x7","x8 ~~ x8","x9 ~~ x9","x10 ~~ x10","x11 ~~ x11","x12 ~~ x12")])
             - diag(S_values[1:length(manif), 1:length(manif)]))^2
          )
        )/length(bic_regsem_param))


      # AIC
      aic_regsem <- regsem(lavFullModel, lambda = best_lambda_aic, pars_pen = "loadings")

      aic_regsem_param <- aic_regsem$out$pars

      parameterValues$regsem_AIC[,iteration] <- as.numeric(aic_regsem_param[c("f1 -> x1","f1 -> x2","f1 -> x3","f1 -> x4","f1 -> x6","f1 -> x7","f1 -> x8","f2 -> x5","f2 -> x6","f2 -> x7","f2 -> x8","f2 -> x10", "f2 -> x11", "f2 -> x12", "f3 -> x2","f3 -> x3","f3 -> x4","f3 -> x9","f3 -> x10","f3 -> x11","f3 -> x12")])

      # compute RMSE
      RMSE$regsem_AIC[,iteration] <- sqrt(
        (sum(
          (as.numeric(aic_regsem_param[c("f1 -> x1","f1 -> x2","f1 -> x3","f1 -> x4","f1 -> x6","f1 -> x7","f1 -> x8","f2 -> x5","f2 -> x6","f2 -> x7","f2 -> x8","f2 -> x10", "f2 -> x11", "f2 -> x12", "f3 -> x2","f3 -> x3","f3 -> x4","f3 -> x9","f3 -> x10","f3 -> x11","f3 -> x12")])
           - c(A_values[c(1:4, 6,7,8),13], A_values[c(5:8, 10,11,12),14], A_values[c(2,3,4, 9:12),15]))^2
        ) +
          sum(
            (as.numeric(aic_regsem_param[c("x1 ~~ x1","x2 ~~ x2","x3 ~~ x3","x4 ~~ x4","x5 ~~ x5","x6 ~~ x6","x7 ~~ x7","x8 ~~ x8","x9 ~~ x9","x10 ~~ x10","x11 ~~ x11","x12 ~~ x12")])
             - diag(S_values[1:length(manif), 1:length(manif)]))^2
          )
        )/length(aic_regsem_param))

    }else{
      error = error+ 1
      next
    }



    ##### with CV ####

    # OpenMx
    ## Note: Cross validation is performed differently from Jacobucci (2016): it is assumed that the researcher only
    # has one sample of site sampleSize and this sample is split into a training and a testing sample


    ### covariance based: ###

    CV_cov_reg_Model <- tryCatch(fitRegModels(model = CV_covMxModel, data_type = "cov",model_type = "mxModel",
                                              fitfun = "FML",pen_on = "A",selectedA = selectedA,pen_start = 0,pen_end = .4,pen_stepsize = .01,
                                              fit_index = "CV_m2LL", CV = T, Test_Sample = mx_cov_test_data) ,

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
      cov_reg_Model_CV_m2LL <- createRegModel(covFullMxModel, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = minimum_CV_cov_m2LL,pen_on = "A",selectedA = selectedA)
      cov_reg_Model_CV_m2LL <- mxRun(cov_reg_Model_CV_m2LL, silent = T)

      # save parameters

      parameterValues$CV_m2LL[,iteration] <- cvectorize(cov_reg_Model_CV_m2LL$BaseModel$A$values[A_estimated])

      # compute RMSE
      RMSE$CV_m2LL[,iteration] <- computeRMSE(A_est = cov_reg_Model_CV_m2LL$BaseModel$A$values, S_est = cov_reg_Model_CV_m2LL$BaseModel$S$values,
                                              Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]



      print("Covariance with CV successful")
    }else{
      error = error+ 1
      next
    }

    ## FIML based models ###

    FIML_CV_reg_Model <- tryCatch(fitRegModels(model = FIML_CV_rawMxModel, data_type = "raw",model_type = "mxModel",
                                               fitfun = "FIML",pen_on = "A",selectedA = selectedA,pen_start = 0,pen_end = .4,pen_stepsize = .01,
                                               fit_index = "CV_m2LL", CV = T, Test_Sample = FIML_test_data),


                                  error = function(e){
                                    print("warning: did not find proper solution")
                                    return(NA)
                                  }
    )

    if(!is.logical(FIML_CV_reg_Model)){

      FIML_minimum_CV_m2LL <- FIML_CV_reg_Model$`fit measures`[which(FIML_CV_reg_Model$`fit measures`[,"CV m2LL"]==min(FIML_CV_reg_Model$`fit measures`[,"CV m2LL"])),"penalty"]

      FIML_reg_Model_CV_m2LL <- createRegModel(model = rawMxModel, data_type = "raw",model_type = "mxModel",fitfun = "FIML",pen_value = FIML_minimum_CV_m2LL,pen_on = "A",selectedA = selectedA)
      FIML_reg_Model_CV_m2LL <- mxRun(FIML_reg_Model_CV_m2LL, silent = T)

      # save parameters

      parameterValues$FIML_CV_m2LL[,iteration] <- cvectorize(FIML_reg_Model_CV_m2LL$BaseModel$A$values[A_estimated])

      RMSE$FIML_CV_m2LL[,iteration] <- computeRMSE(A_est = FIML_reg_Model_CV_m2LL$BaseModel$A$values, S_est = FIML_reg_Model_CV_m2LL$BaseModel$S$values,
                                                   Apop = A_values, Spop = S_values, Afree = A_estimated, Sfree = S_estimated)[3,1]


      print("FIML with CV successful")
    }else{
      error = error+ 1
      next
    }

    # regsem
    # note: a difficulty with regsem is that the variable order in the model implied covariance matrix is often different from the one in the data set provided
    # step 1: adjust variable order in cv-sample:
    reorderedCVSample <- test_raw_data[,rownames(fitted(lavTrainModel)$cov)]
    reorderedCVSampleCov <- cov(reorderedCVSample)

    CV_regsem_out <- tryCatch(cv_regsem(lavTrainModel, n.lambda = 40, pars_pen = "loadings", metric = "chisq", fit.ret2 = "test", test.cov = reorderedCVSampleCov, test.n.obs = nrow(reorderedCVSample)),
                              warning = function(w){
                                print("warning: did not find proper solution")
                                return(NA)
                              },
                              error = function(e){
                                print("warning: did not find proper solution")
                                return(NA)
                              }
    )
    if(!is.logical(CV_regsem_out)){

      # fit best model on full data:
      best_lambda_CV_chisq <- CV_regsem_out$fits[which(CV_regsem_out$fits[,"chisq"]== min(CV_regsem_out$fits[,"chisq"], na.rm = T)),"lambda"]

      # save parameters:
      # BIC
      CV_chisq_regsem <- regsem(lavFullModel, lambda = best_lambda_CV_chisq, pars_pen = "loadings")
      CV_chisq_regsem_param <- CV_chisq_regsem$out$pars
      parameterValues$CV_regsem_chisqu[,iteration] <- as.numeric(CV_chisq_regsem_param[c("f1 -> x1","f1 -> x2","f1 -> x3","f1 -> x4","f1 -> x6","f1 -> x7","f1 -> x8","f2 -> x5","f2 -> x6","f2 -> x7","f2 -> x8","f2 -> x10", "f2 -> x11", "f2 -> x12", "f3 -> x2","f3 -> x3","f3 -> x4","f3 -> x9","f3 -> x10","f3 -> x11","f3 -> x12")])

      # compute RMSE
      RMSE$CV_regsem_chisqu[,iteration] <- sqrt(
        (sum(
          (as.numeric(CV_chisq_regsem_param[c("f1 -> x1","f1 -> x2","f1 -> x3","f1 -> x4","f1 -> x6","f1 -> x7","f1 -> x8","f2 -> x5","f2 -> x6","f2 -> x7","f2 -> x8","f2 -> x10", "f2 -> x11", "f2 -> x12", "f3 -> x2","f3 -> x3","f3 -> x4","f3 -> x9","f3 -> x10","f3 -> x11","f3 -> x12")])
           - c(A_values[c(1:4, 6,7,8),13], A_values[c(5:8, 10,11,12),14], A_values[c(2,3,4, 9:12),15]))^2
        ) +
          sum(
            (as.numeric(CV_chisq_regsem_param[c("x1 ~~ x1","x2 ~~ x2","x3 ~~ x3","x4 ~~ x4","x5 ~~ x5","x6 ~~ x6","x7 ~~ x7","x8 ~~ x8","x9 ~~ x9","x10 ~~ x10","x11 ~~ x11","x12 ~~ x12")])
             - diag(S_values[1:length(manif), 1:length(manif)]))^2
          )
        )/length(CV_chisq_regsem_param))


      print("regsem with CV successful")

    }else{
      error = error+ 1
      next
    }








    # set progress bar:
    setTkProgressBar(global_pb, iteration, label=paste( round(iteration/total_repetitions*100, 0),
                                                        "% done"))

    iteration = iteration+1

    if(iteration%%5 == 0){
      print(iteration)
    }







  }

  setwd(wd)
  save(overall_evaluation, parameterValues, RMSE, sampleSize, error, seed, total_repetitions, file = paste("Simulation1_N", sampleSize, "_1_1000.RData", sep = ""))

}




