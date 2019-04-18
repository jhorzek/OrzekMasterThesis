#' Simulation study 2
#'
#'
#'
#' @param sampleSize sampleSize
#' @param gamma gamma
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
simulation2 <- function(sampleSize, seed, wd, gamma){
  setwd(wd)
  set.seed(seed)

  ##### Settings #####

  gamma = gamma
  sampleSize = sampleSize

  ##### prepare result data sets #####
  total_repetitions = 1000
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
                                   "FIML_CV_m2LL_false_zero"=rep(NA,total_repetitions),

                                   "regsem_AIC_false_nonzero"=rep(NA,total_repetitions),
                                   "regsem_AIC_false_zero"=rep(NA,total_repetitions),
                                   "regsem_BIC_false_nonzero"=rep(NA,total_repetitions),
                                   "regsem_BIC_false_zero"=rep(NA,total_repetitions),
                                   "regsem_CV_chisq_false_nonzero"=rep(NA,total_repetitions),
                                   "regsem_CV_chisq_false_zero"=rep(NA,total_repetitions)
  )

  parameter_Table <- matrix(NA, nrow = 1, ncol = total_repetitions)

  rownames(parameter_Table) <- "gamma"

  parameterValues <- list("Base_Model_cov" = parameter_Table,
                          "Base_Model_FIML" = parameter_Table,
                          "Base_Model_lavaan" = parameter_Table,
                          "AIC"= parameter_Table,
                          "BIC"= parameter_Table,
                          "FIML_BIC"= parameter_Table,
                          "FIML_AIC"= parameter_Table,
                          "CV_m2LL" = parameter_Table,
                          "FIML_CV_m2LL" = parameter_Table,
                          "regsem_AIC" = parameter_Table,
                          "regsem_BIC" = parameter_Table,
                          "regsem_CV_chisq" = parameter_Table)

  RMSE_table <- matrix(NA, nrow = 3, ncol = total_repetitions)
  rownames(RMSE_table) <- c("RMSE_Amat", "RMSE_Smat", "RMSE_comb")
  RMSE <- list("Base_Model_cov" = RMSE_table,
               "Base_Model_FIML" = RMSE_table,
               "Base_Model_lavaan" = RMSE_table,
               "AIC"= RMSE_table,
               "BIC"= RMSE_table,
               "FIML_BIC"= RMSE_table,
               "FIML_AIC"= RMSE_table,
               "CV_m2LL" = RMSE_table,
               "FIML_CV_m2LL" = RMSE_table,
               "regsem_AIC" = matrix(NA, nrow = 1, ncol = total_repetitions),
               "regsem_BIC" = matrix(NA, nrow = 1, ncol = total_repetitions),
               "regsem_CV_chisq" = matrix(NA, nrow = 1, ncol = total_repetitions))
  ##### start simulation

  # Progress bar
  global_pb <- tkProgressBar(title = "progress bar", min = 0,
                             max = total_repetitions, width = 300) # progress - bar

  while(iteration <= total_repetitions){

    # simulate data:
    SimModel <-paste('# latent variables
f1 =~ .9*x1 + .8*x2 + .7*x3 + .6*x4
f2 =~ .9*y1 + .8*y2 + .7*y3 + .6*y4

# regressions
f2 ~ ', gamma, '*f1

# covariances
f1 ~~ 1*f1
f2 ~~ (1- ', gamma, '^2)*f2

x1 ~~ (1-.9^2)*x1
x2 ~~ (1-.8^2)*x2
x3 ~~ (1-.7^2)*x3
x4 ~~ (1-.6^2)*x4
y1 ~~ (1-.9^2)*y1
y2 ~~ (1-.8^2)*y2
y3 ~~ (1-.7^2)*y3
y4 ~~ (1-.6^2)*y4
',sep = "")


    full_raw_data <- simulateData(SimModel, meanstructure = F, sample.nobs = sampleSize)
    full_raw_data <- scale(full_raw_data, center = T, scale = F)

    # split dataset in train and test dataset
    Folds <- createFolds(c(1:sampleSize),2)

    train_raw_data <- full_raw_data[Folds$Fold1,]
    test_raw_data <- full_raw_data[Folds$Fold2,]

    train_raw_data <- scale(train_raw_data)
    test_raw_data <- scale(test_raw_data)


    ##### define lavaan model #####
    FullLavaanModel <- paste('
                           # covariances
                           f1 ~~ 1*f1
                           f2 ~~ ',1- gamma^2,'*f2

                           # latent variables
                           f1 =~ NA*x1 + x2 + x3 + x4
                           f2 =~ NA*y1 + y2 + y3 + y4

                           # regressions
                           f2 ~ gamma*f1
                           ', sep = "")

    ##### define OpenMx model #####
    latent = c('f1', 'f2')
    manifest = c('x1', 'x2', 'x3', 'x4', 'y1', 'y2', 'y3', 'y4')


    f1path <- mxPath(from = "f1", to = c('x1', 'x2', 'x3', 'x4'), values = 1, free = T)
    f2path <- mxPath(from = "f2", to = c('y1', 'y2', 'y3', 'y4'), values = 1, free = T)
    lpath <- mxPath(from = "f1", to = "f2", values = 1, free = T, labels = "gamma")

    lcov <- mxPath(from = latent, values = c(1, 1-gamma^2), free = F, arrows = 2)
    mcov <- mxPath(from = manifest, values = 1, free = T, arrows = 2)

    #population parameters
    popAvalues <- matrix(c(rep(0,8), .9,0,
                           rep(0,8), .8,0,
                           rep(0,8), .7,0,
                           rep(0,8), .6,0,
                           rep(0,8), 0,.9,
                           rep(0,8), 0,.8,
                           rep(0,8), 0,.7,
                           rep(0,8), 0,.6,
                           rep(0,10),
                           rep(0,8), gamma,0), ncol = 10, byrow = T)
    rownames(popAvalues) <- c(manifest,latent)
    colnames(popAvalues) <- c(manifest,latent)


    popSvalues <- matrix(0, 10,10)
    diag(popSvalues) <- c(.19,.36,.51,.64,.19,.36,.51,.64,1,1-gamma^2)
    rownames(popSvalues) <- c(manifest,latent)
    colnames(popSvalues) <- c(manifest,latent)

    Afree <- !popAvalues == 0
    rownames(Afree) <- c(manifest,latent)
    colnames(Afree) <- c(manifest,latent)
    Afree["f2", "f1"] <- TRUE

    Sfree <- !popSvalues == 0
    rownames(Sfree) <- c(manifest,latent)
    colnames(Sfree) <- c(manifest,latent)
    Sfree["f1", "f1"] <- FALSE
    Sfree["f2", "f2"] <- FALSE

    regsem_pop_Param <- c(popAvalues["x1", "f1"], popAvalues["x2", "f1"], popAvalues["x3", "f1"], popAvalues["x4", "f1"],
                          popAvalues["y1", "f2"], popAvalues["y2", "f2"], popAvalues["y3", "f2"], popAvalues["y4", "f2"],
                          popAvalues["f2", "f1"],
                          popSvalues[Sfree])

    ##### build models with lavaan #####
    full_fitLavaan <- tryCatch(cfa(FullLavaanModel, sample.cov = cov(full_raw_data), sample.nobs = sampleSize),
                               warning = function(w){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               },
                               error = function(e){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               }
    )
    if(is.logical(full_fitLavaan)){
      error = error+ 1
      next
    }


    ##### build models with openMx:#####





    # covariance based:

    full_fitMxARCL_cov <- tryCatch(mxTryHard(mxModel(model = "Full_mxModel", type = "RAM", manifestVars = manifest, latentVars = latent,
                                                     f1path, f2path, lpath, lcov, mcov,
                                                     mxData(observed = cov(full_raw_data), type = "cov", numObs = sampleSize))),
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

    full_fitMxARCL_FIML <- tryCatch(mxTryHard(mxModel(model = "Full_mxModel", type = "RAM", manifestVars = manifest, latentVars = latent,
                                                      f1path, f2path, lpath, lcov, mcov,
                                                      mxPath(from = "one", to = manifest, values = 0, free = F),
                                                      mxData(observed = full_raw_data, type = "raw"))),
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

    ##### with CV #####
    ## lavaan
    cv_fitLavaan <- tryCatch(cfa(FullLavaanModel, sample.cov = cov(train_raw_data), sample.nobs = length(Folds$Fold1)),
                             warning = function(w){
                               print("warning: did not find proper solution")
                               return(NA)
                             },
                             error = function(e){
                               print("warning: did not find proper solution")
                               return(NA)
                             }
    )
    if(is.logical(cv_fitLavaan)){
      error = error+ 1
      next
    }

    ## OpenMx
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
      CV_covMxModel <- full_fitMxARCL_cov
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
      FIML_CV_rawMxModel <- full_fitMxARCL_FIML
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
    parameterValues$Base_Model_lavaan[,iteration] = full_fitLavaan@ParTable$est[full_fitLavaan@ParTable$label == "gamma"]
    parameterValues$Base_Model_cov[,iteration] = full_fitMxARCL_cov$A$values[which(full_fitMxARCL_cov$A$labels == "gamma")]
    parameterValues$Base_Model_FIML[,iteration] = full_fitMxARCL_FIML$A$values[which(full_fitMxARCL_FIML$A$labels == "gamma")]

    RMSE$Base_Model_lavaan[,iteration] <- computeRMSE(A_est = extractMatrices(full_fitLavaan)$A_est,
                                                      S_est = extractMatrices(full_fitLavaan)$S_est,
                                                      Apop = popAvalues,
                                                      Spop = popSvalues,
                                                      Afree = Afree,
                                                      Sfree = Sfree)
    RMSE$Base_Model_cov[,iteration] <- computeRMSE(A_est = full_fitMxARCL_cov$A$values,
                                                   S_est = full_fitMxARCL_cov$S$values,
                                                   Apop = popAvalues,
                                                   Spop = popSvalues,
                                                   Afree = Afree,
                                                   Sfree = Sfree)
    RMSE$Base_Model_FIML[,iteration] <- computeRMSE(A_est = full_fitMxARCL_FIML$A$values,
                                                    S_est = full_fitMxARCL_FIML$S$values,
                                                    Apop = popAvalues,
                                                    Spop = popSvalues,
                                                    Afree = Afree,
                                                    Sfree = Sfree)

    ###### computation without CV ######
    # regsem
    full_lavaan_reg_Model <- tryCatch(cv_regsem(model = full_fitLavaan, n.lambda = 40, pars_pen = "gamma",type = "lasso", metric = "BIC", fit.ret = c("AIC", "BIC")),

                                      error = function(e){
                                        print("warning: did not find proper solution")
                                        return(NA)
                                      }
    )
    if(!is.logical(full_lavaan_reg_Model)){
      parameterValues$regsem_BIC[,iteration] <- full_lavaan_reg_Model$final_pars["f1 -> f2"]
      RMSE$regsem_BIC[iteration] <- sqrt(mean((full_lavaan_reg_Model$final_pars - regsem_pop_Param)^2))

      AICrow <- which(full_lavaan_reg_Model$fits[,"AIC"] == min(full_lavaan_reg_Model$fits[,"AIC"], na.rm = T))
      parameterValues$regsem_AIC[,iteration] <- full_lavaan_reg_Model$parameters[AICrow, "f1 -> f2"]
      RMSE$regsem_AIC[iteration] <- sqrt(mean((full_lavaan_reg_Model$parameters[AICrow, ] - regsem_pop_Param)^2))

      print("regsem without CV successful")
    }

    # OpenMx
    ##### Covariance Based Models ####
    ## start fitting with different penalty values:

    penA <- matrix(0, ncol = ncol(full_fitMxARCL_cov$A$values), nrow = nrow(full_fitMxARCL_cov$A$values))
    penA[full_fitMxARCL_cov$A$labels == "gamma"] <- 1

    full_cov_reg_Model <- tryCatch(fitRegModels(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",
                                                fitfun = "FML",pen_on = "A",selectedA = penA,pen_start = 0,pen_end = .4,pen_stepsize = .01,fit_index = "BIC"),

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

      full_cov_reg_Model_AIC <- createRegModel(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = cov_minimum_AIC,pen_on = "A",selectedA = penA)
      full_cov_reg_Model_AIC <- mxRun(full_cov_reg_Model_AIC, silent = T)


      full_cov_reg_Model_BIC <- createRegModel(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = cov_minimum_BIC,pen_on = "A",selectedA = penA)
      full_cov_reg_Model_BIC <- mxRun(full_cov_reg_Model_BIC, silent = T)

      #### saving parameters
      parameterValues$AIC[,iteration] <- full_cov_reg_Model_AIC$BaseModel$A$values[which(full_fitMxARCL_cov$A$labels == "gamma")]
      parameterValues$BIC[,iteration] <- full_cov_reg_Model_BIC$BaseModel$A$values[which(full_fitMxARCL_cov$A$labels == "gamma")]
      RMSE$AIC[,iteration] <- computeRMSE(A_est = full_cov_reg_Model_AIC$BaseModel$A$values,
                                          S_est = full_cov_reg_Model_AIC$BaseModel$S$values,
                                          Apop = popAvalues,
                                          Spop = popSvalues,
                                          Afree = Afree,
                                          Sfree = Sfree)

      RMSE$BIC[,iteration] <- computeRMSE(A_est = full_cov_reg_Model_BIC$BaseModel$A$values,
                                          S_est = full_cov_reg_Model_BIC$BaseModel$S$values,
                                          Apop = popAvalues,
                                          Spop = popSvalues,
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
                                                 fitfun = "FIML",pen_on = "A",selectedA = penA,pen_start = 0,pen_end = .4,pen_stepsize = .01,fit_index = "BIC"),
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
      full_FIML_reg_Model_AIC <- createRegModel(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel", fitfun = "FIML",pen_value = FIML_minimum_AIC,pen_on = "A",selectedA = penA)
      full_FIML_reg_Model_AIC <- mxRun(full_FIML_reg_Model_AIC, silent = T)

      full_FIML_reg_Model_BIC <- createRegModel(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel",fitfun = "FIML",pen_value = FIML_minimum_BIC,pen_on = "A",selectedA = penA)
      full_FIML_reg_Model_BIC <- mxRun(full_FIML_reg_Model_BIC, silent = T)


      # save parameters

      parameterValues$FIML_AIC[,iteration] <- full_FIML_reg_Model_AIC$BaseModel$A$values[which(full_fitMxARCL_cov$A$labels == "gamma")]
      parameterValues$FIML_BIC[,iteration] <- full_FIML_reg_Model_BIC$BaseModel$A$values[which(full_fitMxARCL_cov$A$labels == "gamma")]

      RMSE$FIML_AIC[,iteration] <- computeRMSE(A_est = full_FIML_reg_Model_AIC$BaseModel$A$values,
                                               S_est = full_FIML_reg_Model_AIC$BaseModel$S$values,
                                               Apop = popAvalues,
                                               Spop = popSvalues,
                                               Afree = Afree,
                                               Sfree = Sfree)

      RMSE$FIML_BIC[,iteration] <- computeRMSE(A_est = full_FIML_reg_Model_BIC$BaseModel$A$values,
                                               S_est = full_FIML_reg_Model_BIC$BaseModel$S$values,
                                               Apop = popAvalues,
                                               Spop = popSvalues,
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
                                              fitfun = "FML",pen_on = "A",selectedA = penA,pen_start = 0,pen_end = .4,pen_stepsize = .01,
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
      cov_reg_Model_CV_m2LL <- createRegModel(full_fitMxARCL_cov, data_type = "cov",model_type = "mxModel",fitfun = "FML",pen_value = minimum_CV_cov_m2LL,pen_on = "A",selectedA = penA)
      cov_reg_Model_CV_m2LL <- mxRun(cov_reg_Model_CV_m2LL, silent = T)


      # save parameters

      parameterValues$CV_m2LL[,iteration] <-cov_reg_Model_CV_m2LL$BaseModel$A$values[which(CV_covMxModel$A$labels == "gamma")]
      RMSE$CV_m2LL[,iteration] <- computeRMSE(A_est = cov_reg_Model_CV_m2LL$BaseModel$A$values,
                                              S_est = cov_reg_Model_CV_m2LL$BaseModel$S$values,
                                              Apop = popAvalues,
                                              Spop = popSvalues,
                                              Afree = Afree,
                                              Sfree = Sfree)

      print("Covariance with CV successful")
    }else{
      error = error+ 1
      next
    }

    ## FIML based models ###

    FIML_CV_reg_Model <- tryCatch(fitRegModels(model = FIML_CV_rawMxModel, data_type = "raw",model_type = "mxModel",
                                               fitfun = "FIML",pen_on = "A",selectedA = penA,pen_start = 0,pen_end = .4,pen_stepsize = .01,
                                               fit_index = "BIC", CV = T, Test_Sample = FIML_test_data),


                                  error = function(e){
                                    print("warning: did not find proper solution")
                                    return(NA)
                                  }
    )

    if(!is.logical(FIML_CV_reg_Model)){

      FIML_minimum_CV_m2LL <- FIML_CV_reg_Model$`fit measures`[which(FIML_CV_reg_Model$`fit measures`[,"CV m2LL"]==min(FIML_CV_reg_Model$`fit measures`[,"CV m2LL"])),"penalty"]

      FIML_reg_Model_CV_m2LL <- createRegModel(model = full_fitMxARCL_FIML, data_type = "raw",model_type = "mxModel",fitfun = "FIML",pen_value = FIML_minimum_CV_m2LL,pen_on = "A",selectedA = penA)
      FIML_reg_Model_CV_m2LL <- mxRun(FIML_reg_Model_CV_m2LL, silent = T)

      # save parameters

      parameterValues$FIML_CV_m2LL[,iteration] <- FIML_reg_Model_CV_m2LL$BaseModel$A$values[which(CV_covMxModel$A$labels == "gamma")]
      RMSE$FIML_CV_m2LL[,iteration] <- computeRMSE(A_est = FIML_reg_Model_CV_m2LL$BaseModel$A$values,
                                                   S_est = FIML_reg_Model_CV_m2LL$BaseModel$S$values,
                                                   Apop = popAvalues,
                                                   Spop = popSvalues,
                                                   Afree = Afree,
                                                   Sfree = Sfree)


      print("FIML with CV successful")
    }else{
      error = error+ 1
      next
    }

    # with regsem
    # note: a difficulty with regsem is that the variable order in the model implied covariance matrix is often different from the one in the data set provided
    # step 1: adjust variable order in cv-sample:
    reorderedCVSample <- test_raw_data[,rownames(fitted(cv_fitLavaan)$cov)]
    reorderedCVSampleCov <- cov(reorderedCVSample)

    CV_regsem_out <- tryCatch(cv_regsem(cv_fitLavaan, n.lambda = 40, pars_pen = "gamma", metric = "chisq", fit.ret2 = "test", test.cov = reorderedCVSampleCov, test.n.obs = nrow(reorderedCVSample)),
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

      print("regsem with CV successful")

      # extract best penalty
      lambda <- CV_regsem_out$fits[which(CV_regsem_out$fits[,"chisq"] == min(CV_regsem_out$fits[,"chisq"], na.rm = T)),"lambda"]

      cv_lavaan_model <- regsem(model = full_fitLavaan, lambda = lambda, pars_pen = "gamma", type = "lasso")
      # save parameters:
      # chisqu
      parameterValues$regsem_CV_chisq[,iteration] <- as.numeric(cv_lavaan_model$coefficients["f1 -> f2"])
      RMSE$regsem_CV_chisq[iteration] <- sqrt(mean((cv_lavaan_model$out$pars - regsem_pop_Param)^2))
      print("regsem with CV successful")

    }else{
      error = error+ 1
      next
    }

    # set progress bar:
    setTkProgressBar(global_pb, iteration, label=paste( round(iteration/total_repetitions*100, 0),
                                                        "% done"))

    iteration = iteration+1


  }





  ###########################
  # overall evaluataion:

  for(i in 1:total_repetitions){
    # a12 is non-zero in true model:
    if(gamma>0){
      overall_evaluation$AIC_false_zero[i] <- sum(abs(parameterValues$AIC[,i]) < .001)
      overall_evaluation$FIML_AIC_false_zero[i] <- sum(abs(parameterValues$FIML_AIC[,i]) < .001)
      overall_evaluation$BIC_false_zero[i] <- sum(abs(parameterValues$BIC[,i]) < .001)
      overall_evaluation$FIML_BIC_false_zero[i] <- sum(abs(parameterValues$FIML_BIC[,i]) < .001)
      overall_evaluation$CV_m2LL_false_zero[i] <- sum(abs(parameterValues$CV_m2LL[,i]) < .001)
      overall_evaluation$FIML_CV_m2LL_false_zero[i] <- sum(abs(parameterValues$FIML_CV_m2LL[,i]) < .001)
      overall_evaluation$regsem_AIC_false_zero[i] <- sum(abs(parameterValues$regsem_AIC[,i]) < .001)
      overall_evaluation$regsem_BIC_false_zero[i] <- sum(abs(parameterValues$regsem_BIC[,i]) < .001)
      overall_evaluation$regsem_CV_chisq_false_zero[i] <- sum(abs(parameterValues$regsem_CV_chisq[,i]) < .001)
    }else{
      # a21 is zero in rue model

      overall_evaluation$AIC_false_nonzero[i] <-  sum(abs(parameterValues$AIC[,i]) > .001)
      overall_evaluation$FIML_AIC_false_nonzero[i] <- sum(abs(parameterValues$FIML_AIC[,i]) > .001)
      overall_evaluation$BIC_false_nonzero[i] <- sum(abs(parameterValues$BIC[,i]) > .001)
      overall_evaluation$FIML_BIC_false_nonzero[i] <- sum(abs(parameterValues$FIML_BIC[,i]) > .001)
      overall_evaluation$CV_m2LL_false_nonzero[i] <- sum(abs(parameterValues$CV_m2LL[,i]) > .001)
      overall_evaluation$FIML_CV_m2LL_false_nonzero[i] <- sum(abs(parameterValues$FIML_CV_m2LL[,i]) > .001)
      overall_evaluation$regsem_AIC_false_nonzero[i] <- sum(abs(parameterValues$regsem_AIC[,i]) > .001)
      overall_evaluation$regsem_BIC_false_nonzero[i] <- sum(abs(parameterValues$regsem_BIC[,i]) > .001)
      overall_evaluation$regsem_CV_chisq_false_nonzero[i] <- sum(abs(parameterValues$regsem_CV_chisq[,i]) > .001)
    }

  }


  sumsOverallEval <- apply(overall_evaluation,2, sum, na.rm = TRUE)
  names(sumsOverallEval) <- colnames(overall_evaluation)
  perRunOverAll <- sumsOverallEval/total_repetitions


  # note: regularization makes fit worse. The reason is that both parameters, the non-zero and the true-zero one are regularized. The non-zero one
  # gets pulled away from its true value and this impacts the over-all RMSEA more than setting a parameter that is close to zero to zero

  #save(overall_evaluation, parameterValues,total_repetitions, RMSE, gamma, error, file = "Study2_N200_1_1000.RData")

  save(gamma,overall_evaluation, parameterValues, RMSE, sampleSize, error, seed, total_repetitions, seed, file = paste("Simulation2_N", sampleSize, "_gamma",gamma,"_1_1000.RData", sep = ""))

}
