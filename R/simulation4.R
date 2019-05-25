#' Simulation study 4
#'
#'
#'
#' @param sampleSize sampleSize
#' @param seed seed
#' @param wd working directory. The results will be saved here
#' @param total_repetitions total number of repetitions
#' @author Jannik Orzek
#' @import OpenMx caret laremm regsem tcltk ctsem
#'
#' @examples
#'
#' @export
#'
#'
simulation4 <- function(sampleSize, seed, wd, total_repetitions, Scale){
  setwd(wd)
  ##### Settings #####

  numLatent = 2
  Timepoints = 5
  burning = 50
  sampleSize = sampleSize
  set.seed(seed)

  ##### Parameters #####

  Adisc <- matrix(c(.5,.2,0,.5), ncol = 2,byrow = T)

  DRIFT <- logm(Adisc)
  DRIFT_hash <- kronecker(DRIFT, diag(2)) + kronecker(diag(2),DRIFT)

  atempsol <- solve(DRIFT_hash)%*%(expm(DRIFT_hash)-diag(4))

  Q4 <- .75/atempsol[4,4]
  Q3 <- (0 - atempsol[3,4]*Q4)/atempsol[3,3]
  Q2 <- (0 - atempsol[2,4]*Q4)/atempsol[2,2]
  Q1 <- (.71 - atempsol[1,2]*Q2  - atempsol[1,3]*Q3 - atempsol[1,4]*Q4 )/atempsol[1,1]
  Qmat <- matrix(c(Q1,Q2,Q3,Q4), nrow = 2, byrow = F)
  DIFFUSION <- t(chol(Qmat))

  MANIFESTVAR=diag(0,2) # items have 0 measurement-error-variance
  LAMBDA=diag(1,2) # loadings set to 1
  TRAITVAR=matrix(c(0,0,0,0),nrow=2) # no trait variance
  CINT=matrix(c(0,0),nrow=2) # intercepts set to zero
  T0MEANS=matrix(0,ncol=1,nrow=2) # initial means set to zero
  T0VAR=diag(1,2) # initial variance set to 1

  generatingModel<-ctModel(Tpoints= Timepoints,n.latent=numLatent,n.TDpred=0,n.TIpred=0,n.manifest=2,
                           MANIFESTVAR=MANIFESTVAR, # items have 0 measurement-error-variance
                           LAMBDA=LAMBDA, # loadings set to 1
                           DRIFT=DRIFT,
                           TRAITVAR=TRAITVAR,
                           DIFFUSION=DIFFUSION, # diffusion matrix G
                           CINT=CINT, # intercepts set to zero
                           T0MEANS=T0MEANS, # initial means set to zero
                           T0VAR=T0VAR # initial variance set to 1
  )

  ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
  # although it is 0 in the population.

  analysisModel <- ctModel(Tpoints= Timepoints,n.latent=numLatent,n.TDpred=0,n.TIpred=0,n.manifest=2,
                           LAMBDA=diag(1,2),
                           MANIFESTVAR=diag(0,2),
                           CINT=matrix(c(0,0),nrow=2),
                           MANIFESTMEANS = matrix(c(0,0),nrow=2),
                           DIFFUSION=matrix(c('eta1_eta1', 'eta2_eta1','eta1_eta2','eta2_eta2'),2),
                           T0MEANS=matrix(0,ncol=1,nrow=2),
                           T0VAR=diag(1,2))

  ##### prepare result data sets #####
  total_repetitions = total_repetitions
  iteration = 1
  improper_solutions = 0
  error = 0

  overall_evaluation <- data.frame(
    "FIML_AIC_false_nonzero" = rep(NA,total_repetitions),
    "FIML_AIC_false_zero"=rep(NA,total_repetitions) ,
    "FIML_BIC_false_nonzero"=rep(NA,total_repetitions),
    "FIML_BIC_false_zero"=rep(NA,total_repetitions),
    "FIML_CV_m2LL_false_nonzero"=rep(NA,total_repetitions),
    "FIML_CV_m2LL_false_zero"=rep(NA,total_repetitions)
  )

  parameter_Table <- matrix(NA, nrow = length(cvectorize(DRIFT)), ncol = total_repetitions)

  rownames(parameter_Table) <- c("eta1 -> eta1", "eta1 -> eta2", "eta2 -> eta1", "eta2 -> eta2")

  parameterValues <- list(
    "Base_Model_FIML" = parameter_Table,
    "FIML_BIC"=parameter_Table, "FIML_AIC"=parameter_Table,
    "FIML_CV_m2LL" =parameter_Table
  )


  RMSE_table <- matrix(NA, nrow = 3, ncol = total_repetitions)
  rownames(RMSE_table) <- c("RMSE_DRIFT", "RMSE_DIFFUSION", "RMSE_comb")
  RMSE <- list("Base_Model_FIML" = RMSE_table,
               "FIML_BIC"=RMSE_table, "FIML_AIC"=RMSE_table,
               "FIML_CV_m2LL" =RMSE_table)

  ##### start simulation

  # Progress bar
  global_pb <- tkProgressBar(title = "progress bar", min = 0,
                             max = total_repetitions, width = 300) # progress - bar

  while(iteration <= total_repetitions){

    # simulate CT-SEM data:

    full_raw_data <- ctGenerate(generatingModel,n.subjects = sampleSize, burnin = burning)

    # split dataset in train and test dataset
    Folds <- createCVfolds(full_raw_data,2)

    train_raw_data <- full_raw_data[Folds$Fold1,]

    test_raw_data <- full_raw_data[Folds$Fold2,]

    if(Scale == T){
      train_raw_data[,c(1,3,5,7,9)] <- scale(train_raw_data[,c(1,3,5,7,9)])
      train_raw_data[,c(2,4,6,8,10)] <- scale(train_raw_data[,c(2,4,6,8,10)])

      test_raw_data[,c(1,3,5,7,9)] <- scale(test_raw_data[,c(1,3,5,7,9)])
      test_raw_data[,c(2,4,6,8,10)] <- scale(test_raw_data[,c(2,4,6,8,10)])

      full_raw_data[,c(1,3,5,7,9)] <- scale(full_raw_data[,c(1,3,5,7,9)])
      full_raw_data[,c(2,4,6,8,10)] <- scale(full_raw_data[,c(2,4,6,8,10)])
    }
    ##### build models with ctsem:#####

    full_fit_ctsem <- tryCatch(ctFit(full_raw_data, analysisModel),
                               warning = function(w){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               },
                               error = function(e){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               }
    )
    if(is.logical(full_fit_ctsem)){
      error = error+ 1
      next
    }

    ## with CV:
    # run Models:
    CV_train_ctsem <- tryCatch(ctFit(train_raw_data, analysisModel),
                               warning = function(w){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               },
                               error = function(e){
                                 print("warning: did not find proper solution")
                                 return(NA)
                               }
    )


    ## check the model fitting for errors #####
    if(is.logical(CV_train_ctsem)){
      error = error+ 1
      next
    }


    # save base parameters
    parameterValues$Base_Model_FIML[,iteration] <- cvectorize(full_fit_ctsem$mxobj$DRIFT$values)

    RMSE$Base_Model_FIML[,iteration] <- computeRMSE(A_est = full_fit_ctsem$mxobj$DRIFT$values,
                                                    S_est = t(chol(full_fit_ctsem$mxobj$DIFFUSION$result)),
                                                    Apop = DRIFT,
                                                    Spop = DIFFUSION,
                                                    Afree = matrix(TRUE,2,2),
                                                    Sfree = matrix(TRUE, 2,2)
    )


    ###### regularized computation without CV ######
    ## start fitting with different penalty values:

    full_FIML_reg_Model <- tryCatch(fitRegModels(model = full_fit_ctsem, model_type = "ctsem",
                                                 fitfun = "FIML",data_type = "raw",
                                                 pen_on = "DRIFT", selectedDrifts = "cross",
                                                 pen_start = 0, pen_end = 1, pen_stepsize = .01, fit_index = "BIC"),
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
      full_FIML_reg_Model_AIC <- createRegModel(model = full_fit_ctsem, model_type = "ctsem",
                                                fitfun = "FIML", data_type = "raw",
                                                pen_on = "DRIFT", selectedDrifts = "cross",
                                                pen_value = FIML_minimum_AIC, driftexpo = T, DRIFT_dt = 1)
      full_FIML_reg_Model_AIC <- mxRun(full_FIML_reg_Model_AIC, silent = T)

      full_FIML_reg_Model_BIC <- createRegModel(model = full_fit_ctsem, model_type = "ctsem",
                                                fitfun = "FIML", data_type = "raw",
                                                pen_on = "DRIFT", selectedDrifts = "cross",
                                                pen_value = FIML_minimum_BIC, driftexpo = T, DRIFT_dt = 1)
      full_FIML_reg_Model_BIC <- mxRun(full_FIML_reg_Model_BIC, silent = T)


      # save parameters

      parameterValues$FIML_AIC[,iteration] <- cvectorize(full_FIML_reg_Model_AIC$BaseModel$DRIFT$values)
      parameterValues$FIML_BIC[,iteration] <- cvectorize(full_FIML_reg_Model_BIC$BaseModel$DRIFT$values)

      RMSE$FIML_AIC[,iteration] <- computeRMSE(A_est = full_FIML_reg_Model_AIC$BaseModel$DRIFT$values,
                                               S_est = t(chol(full_FIML_reg_Model_AIC$BaseModel$DIFFUSION$result)),
                                               Apop = DRIFT,
                                               Spop = DIFFUSION,
                                               Afree = matrix(TRUE,2,2),
                                               Sfree = matrix(TRUE, 2,2))

      RMSE$FIML_BIC[,iteration] <- computeRMSE(A_est = full_FIML_reg_Model_BIC$BaseModel$DRIFT$values,
                                               S_est = t(chol(full_FIML_reg_Model_BIC$BaseModel$DIFFUSION$result)),
                                               Apop = DRIFT,
                                               Spop = DIFFUSION,
                                               Afree = matrix(TRUE,2,2),
                                               Sfree = matrix(TRUE, 2,2))

      print("FIML without CV successful")
    }else{
      error = error+ 1
      next
    }



    ##### with CV ####

    # OpenMx
    ## Note: Cross validation is performed differently from Jacobucci (2016): it is assumed that the researcher only
    # has one sample of site sample_size and this sample is split into a training and a testing sample

    # note: ctsem renames the rows and columns when fitting a model. To get
    # the right names, we fit the ctsem model with the new dataset and then extract
    # the dataset, where rows and columns have now been renamed

    fit_myModel <- ctFit(test_raw_data, analysisModel, useOptimizer = F)
    testdata <- fit_myModel$mxobj$data

    FIML_CV_reg_Model <- tryCatch(fitRegModels(model = CV_train_ctsem, model_type = "ctsem",
                                               fitfun = "FIML",data_type = "raw",
                                               pen_on = "DRIFT", selectedDrifts = "cross",
                                               pen_start = 0, pen_end = 1, pen_stepsize = .01, fit_index = "CV_m2LL",
                                               CV = TRUE, Test_Sample = testdata, driftexpo = T, DRIFT_dt = 1),


                                  error = function(e){
                                    print("warning: did not find proper solution")
                                    return(NA)
                                  }
    )

    if(!is.logical(FIML_CV_reg_Model)){

      # save parameters

      parameterValues$FIML_CV_m2LL[,iteration] <- cvectorize(FIML_CV_reg_Model$bestmodel$DRIFT$values)

      RMSE$FIML_CV_m2LL[,iteration] <- computeRMSE(A_est = FIML_CV_reg_Model$bestmodel$DRIFT$values,
                                                   S_est = t(chol(FIML_CV_reg_Model$bestmodel$DIFFUSION$result)),
                                                   Apop = DRIFT,
                                                   Spop = DIFFUSION,
                                                   Afree = matrix(TRUE,2,2),
                                                   Sfree = matrix(TRUE, 2,2))

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
    overall_evaluation$FIML_AIC_false_zero[i] <- sum(abs(parameterValues$FIML_AIC["eta2 -> eta1",i]) < .001)
    overall_evaluation$FIML_BIC_false_zero[i] <- sum(abs(parameterValues$FIML_BIC["eta2 -> eta1",i]) < .001)
    overall_evaluation$FIML_CV_m2LL_false_zero[i] <- sum(abs(parameterValues$FIML_CV_m2LL["eta2 -> eta1",i]) < .001)

    # a21 is zero in true model

    overall_evaluation$FIML_AIC_false_nonzero[i] <- sum(abs(parameterValues$FIML_AIC["eta1 -> eta2",i]) > .001)
    overall_evaluation$FIML_BIC_false_nonzero[i] <- sum(abs(parameterValues$FIML_BIC["eta1 -> eta2",i]) > .001)
    overall_evaluation$FIML_CV_m2LL_false_nonzero[i] <- sum(abs(parameterValues$FIML_CV_m2LL["eta1 -> eta2",i]) > .001)

  }


  sumsOverallEval <- apply(overall_evaluation,2, sum, na.rm = TRUE)
  names(sumsOverallEval) <- colnames(overall_evaluation)
  perRunOverAll <- sumsOverallEval/total_repetitions


  # note: regularization makes fit worse. The reason is that both parameters, the non-zero and the true-zero one are regularized. The non-zero one
  # gets pulled away from its true value and this impacts the over-all RMSEA more than setting a parameter that is close to zero to zero
  save(overall_evaluation, parameterValues, RMSE, total_repetitions, DRIFT, DIFFUSION, error, seed, file = paste("Simulation4_N", sampleSize, "_1_1000.RData", sep = ""))



}
