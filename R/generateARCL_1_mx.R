#' generateARCL_1_mx
#'
#'
#' generateARCL_1_mx creates an ARCL(1)-SEM model syntax for OpenMx with 2 latent variables and 1 manifest per latent and time point for data simulation
#' Important: the naming of the variables starts at t = 0 and increases to t = timepoints+burning-1. The data-points used in the analysis
#' are: t_burning (corresponds to the first observation in the data set because the indexing starts at 0) to t_(burning+timepoints-1).
#'
#' @param timepoints number of time points observed in the final data set
#' @param burning number of observations to burn before the initial time point
#' @param sample_size sample size
#' @param raw_data raw data set - e.g. simulated with generateARCL_1_lavaan
#'
#' @author Jannik Orzek
#' @import OpenMx
#'
#' @examples
#' library(lavaan)
#' library(OpenMx)
#' a_11 = .5
#' a_22 = .5
#' a_12 = .5
#' a_21 = 0
#' phi1_init = 1
#' phi2_init = 1
#' timepoints = 10
#' burning = 100
#' sample_size = 100
#' seed = 1234
#'
#' # generate lavaan ARCL(1) with two latent variables:
#' lavaanModel <- generateARCL_1(a_11 = a_11, a_22 = a_22,a_12 = a_12,a_21 = a_21,
#'                              phi1_init = phi1_init, phi2_init = phi2_init,
#'                              timepoints = timepoints, burning = burning,
#'                              sample_size = sample_size, seed = seed)
#'
#' myModel <- cfa(lavaanModel$AnalysisModel, sample.cov = cov(lavaanModel$raw_data[,(2*burning+1):(2*burning+2*timepoints)]), sample.nobs = lavaanModel$sample_size)
#' summary(myModel)
#'
#' mymxModel <- generateARCL_1_mx(timepoints = timepoints, burning = burning, sample_size = sample_size, raw_data = lavaanModel$raw_data)
#'
#' fitMxARCL <- mxTryHard(mymxModel)
#' summary(fitMxARCL)
#'
#' @export
#'
#'


generateARCL_1_mx <- function(timepoints, burning, sample_size, raw_data){

  latentVars = c()
  for(time in burning:(burning+timepoints-1)){
    for(lv in 1:2){
      var = paste("eta", lv,"_t", time, sep = "")
      latentVars = c(latentVars,var)
    }
  }

  manifestVars = c()
  for(time in burning:(burning+timepoints-1)){
    for(manif in 1:2){
      man = paste("y", manif,"_t", time, sep = "")
      manifestVars = c(manifestVars,man)
    }
  }

  ##### Amatrix #####
  Avalues = matrix(0, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints+2*timepoints)
  Afree = matrix(FALSE, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints+2*timepoints)
  Alabel = matrix(NA, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints+2*timepoints)
  # loadings of latent on manifest set to 1
  diag(Avalues[1:(2*timepoints), (2*timepoints+1): (2*timepoints+2*timepoints)]) = 1

  # regression coefficients
  Alabel_unit <- matrix(c("a_11", "a_12", "a_21", "a_22"), ncol = 2, byrow = T)
  #initial effects
  Alabel[2*timepoints+3,2*timepoints+1] <- paste("a_11_", burning,sep = "")
  Alabel[2*timepoints+3,2*timepoints+2] <- paste("a_12_", burning,sep = "")
  Alabel[2*timepoints+4,2*timepoints+1] <- paste("a_21_", burning,sep = "")
  Alabel[2*timepoints+4,2*timepoints+2] <- paste("a_22_", burning,sep = "")
  for(rowstart in seq(5,(2*timepoints), by = 2)){
    Alabel[(2*timepoints+rowstart):(2*timepoints+rowstart+1),(2*timepoints+rowstart-2):(2*timepoints+rowstart+1-2)] = Alabel_unit
  }

  # free
  for(rowstart in seq(3,(2*timepoints), by = 2)){
    Afree[(2*timepoints+rowstart):(2*timepoints+rowstart+1),(2*timepoints+rowstart-2):(2*timepoints+rowstart+1-2)] = TRUE
  }

  # build matrix
  Amatrix = mxMatrix(type = "Full", nrow = 2*timepoints+2*timepoints, ncol =2*timepoints+2*timepoints,
                     free = Afree, values = Avalues, labels = Alabel, name = "A", dimnames = list(c(manifestVars, latentVars), c(manifestVars, latentVars)))

  ##### Smatrix #####
  Svalues = matrix(0, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints+2*timepoints)
  Sfree = matrix(FALSE, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints+2*timepoints)
  Slabel = matrix(NA, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints+2*timepoints)

  # latent disturbances
  Slabel_unit <- matrix(c("phi1", NA, NA, "phi2"), ncol = 2, byrow = T)
  #initial disturbance
  Slabel[2*timepoints+1,2*timepoints+1] <- paste("phi1_",burning,sep = "")
  Slabel[2*timepoints+1,2*timepoints+2] <- paste("phi12_",burning,sep = "")
  Slabel[2*timepoints+2,2*timepoints+1] <- paste("phi12_",burning,sep = "")
  Slabel[2*timepoints+2,2*timepoints+2] <- paste("phi2_",burning,sep = "")

  # all other:
  for(rowstart in seq(3,(2*timepoints), by = 2)){
    Slabel[(2*timepoints+rowstart):(2*timepoints+rowstart+1),(2*timepoints+rowstart):(2*timepoints+rowstart+1)] = Slabel_unit
  }
  # last disturbance:
  Slabel[(2*timepoints+2*timepoints-1):(2*timepoints+2*timepoints),(2*timepoints+2*timepoints-1):(2*timepoints+2*timepoints)] <-
    matrix(c("phi1",paste("phi12_",burning+timepoints,sep = ""),
             paste("phi12_",burning+timepoints,sep = ""),"phi2"), ncol = 2, byrow = T)


  # S free
  # initial:
  Sfree[(2*timepoints+1):(2*timepoints+2),(2*timepoints+1):(2*timepoints+2)] = TRUE

  # all other time points:
  Sfree_unit = matrix(c(T,F,F,T), byrow = T)
  for(rowstart in seq(3,(2*timepoints), by = 2)){
    Sfree[(2*timepoints+rowstart):(2*timepoints+rowstart+1),(2*timepoints+rowstart):(2*timepoints+rowstart+1)] = Sfree_unit
  }
  # last time point:
  Sfree[(2*timepoints+2*timepoints-1):(2*timepoints+2*timepoints),(2*timepoints+2*timepoints-1):(2*timepoints+2*timepoints)] = TRUE


  # Svalues
  Svalues[(2*timepoints+1):(2*timepoints+2),(2*timepoints+1):(2*timepoints+2)] = .5
  Svalues[(2*timepoints+2*timepoints-1):(2*timepoints+2*timepoints),(2*timepoints+2*timepoints-1):(2*timepoints+2*timepoints)] = .5
  Svalues_unit = matrix(c(.5,0,0,.5), byrow = T)
  for(rowstart in seq(3,(2*timepoints), by = 2)){
    Svalues[(2*timepoints+rowstart):(2*timepoints+rowstart+1),(2*timepoints+rowstart):(2*timepoints+rowstart+1)] = Svalues_unit
  }

  # build matrix
  Smatrix = mxMatrix(type = "Full", nrow = 2*timepoints+2*timepoints, ncol =2*timepoints+2*timepoints,
                     free = Sfree, values = Svalues, labels = Slabel, name = "S", dimnames = list(c(manifestVars, latentVars), c(manifestVars, latentVars)))

  ##### FMatrix #####
  Fvalues = matrix(0, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints)
  Ffree = matrix(FALSE, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints)
  Flabel = matrix(NA, ncol = 2*timepoints+2*timepoints, nrow = 2*timepoints)

  diag(Fvalues[1:(2*timepoints), 1:(2*timepoints)]) = 1

  Fmatrix = mxMatrix(type = "Full", nrow = 2*timepoints, ncol =2*timepoints+2*timepoints, free =
                       Ffree, values = Fvalues,  name = "F",dimnames = list(manifestVars, c(manifestVars, latentVars) ))

  ##### fit function #####
  expect <- mxExpectationRAM(A="A", S="S", F="F")

  mxARCL <- mxModel(model = "ARCL",
                    mxData(observed = cov(raw_data[,(2*burning+1):(2*burning+2*timepoints)]),type = "cov", numObs =  sample_size),
                    Amatrix, Smatrix, Fmatrix,
                    expect, mxFitFunctionML())
  return(mxARCL)
}

