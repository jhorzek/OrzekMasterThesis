#' generateARCL_1_lavaan
#'
#'
#' generateARCL_1_lavaan creates an ARCL(1)-SEM with 2 latent variables and 1 manifest per latent and time point for data simulation.
#' Important: the naming of the variables starts at t = 0 and increases to t = timepoints+burning-1. The data-points used in the analysis
#' are: t_burning (corresponds to the first observation in the data set because the indexing starts at 0) to t_(burning+timepoints-1).
#'
#' @param a_11 autoregressive parameter for eta_1
#' @param a_22 autoregressive parameter for eta_2
#' @param a_12 cross-lagged parameter for eta_2 -> eta_1
#' @param a_21 cross-lagged parameter for eta_1 -> eta_2
#' @param phi1_init variance of the initial time point for eta1
#' @param phi2_init variance of the initial time point for eta2
#' @param timepoints number of time points observed in the final data set
#' @param burning number of observations to burn before the initial time point
#' @param sample_size sample size
#' @param seed seed for random data simulation
#'
#'
#'#' @author Jannik Orzek
#' @import lavaan
#'
#'
#' @examples
#' library(lavaan)
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
#' lavaanModel <- generateARCL_1_lavaan(a_11 = a_11, a_22 = a_22,a_12 = a_12,a_21 = a_21,
#'                              phi1_init = phi1_init, phi2_init = phi2_init,
#'                              timepoints = timepoints, burning = burning,
#'                              sample_size = sample_size, seed = seed)
#'
#' myModel <- cfa(lavaanModel$AnalysisModel, sample.cov = cov(lavaanModel$raw_data[,(2*burning+1):(2*burning+2*timepoints)]), sample.nobs = lavaanModel$sample_size)
#' summary(myModel)
#'
#' @export
#'
#'
generateARCL_1_lavaan <- function(a_11, a_22, a_12, a_21, phi1_init, phi2_init, timepoints, burning, sample_size, seed){
  tp <- timepoints + burning


  SimulationModel <- ''

  latentEffects <- "#latent regressions"
  for(i in (1):(tp-1)){
    eta1 <- paste("eta1_t", i, " ~ ", a_11, "*eta1_t", i-1, " + ", a_12, "*eta2_t", i-1   ,sep = "")
    eta2 <- paste("eta2_t", i, " ~ ", a_22, "*eta2_t", i-1, " + ", a_21, "*eta1_t", i-1   ,sep = "")
    latentEffects <- paste(latentEffects, eta1, eta2, sep = '\n')
  }

  manifestEffects <- "#definition of latent variables"
  for(i in 0:(tp-1)){
    eta1 <- paste("eta1_t", i, " =~ 1*y1_t", i ,sep = "")
    eta2 <- paste("eta2_t", i, " =~ 1*y2_t", i ,sep = "")
    manifestEffects <- paste(manifestEffects, eta1, eta2, sep = '\n')
  }

  latentDisturbances <- paste("#latent disturbances",
                              paste("eta1_t0 ~~ ",phi1_init,"*eta1_t0", sep = ""),
                              paste("eta2_t0 ~~ ", phi2_init,"*eta2_t0", sep =""),
                              sep = "\n")
  for(i in (1):(tp-1)){
    eta1 <- paste("eta1_t", i, " ~~ (1-1*", a_11, "^2-1*", a_12, "^2)*eta1_t", i ,sep = "")
    eta2 <- paste("eta2_t", i, " ~~ (1-1*", a_22, "^2-1*", a_21, "^2)*eta2_t", i ,sep = "")
    latentDisturbances <- paste(latentDisturbances, eta1, eta2, sep = '\n')
  }

  manifestErrors <- "#manifest errors"
  for(i in 0:(tp-1)){
    y1 <- paste("y1_t", i, " ~~ 0*y1_t", i ,sep = "")
    y2 <- paste("y2_t", i, " ~~ 0*y2_t", i ,sep = "")
    manifestErrors <- paste(manifestErrors, y1, y2, sep = '\n')
  }

  SimulationModel <- paste(latentEffects, manifestEffects, latentDisturbances, manifestErrors, sep = '\n')

  # Simulate data
  set.seed(seed)
  raw_data <- simulateData(model = SimulationModel, sample.nobs = sample_size)

  ###############

  AnalysisModel <- ''

  latentEffects <- paste("#latent regressions")
  #,
  #                       paste("eta1_t", burning+1, " ~ a_11_", burning,"*eta1_t", burning, " + a_12_", burning,"*eta2_t", burning, sep = ""),
  #                       paste("eta2_t", burning+1, " ~ a_22_", burning,"*eta2_t", burning, " + a_21_", burning,"*eta1_t", burning, sep = ""),
  #                       sep = "\n")
  for(i in (burning+1):(tp-1)){
    eta1 <- paste("eta1_t", i, " ~ a_11*eta1_t", i-1, " + a_12*eta2_t", i-1   ,sep = "")
    eta2 <- paste("eta2_t", i, " ~ a_22*eta2_t", i-1, " + a_21*eta1_t", i-1   ,sep = "")
    latentEffects <- paste(latentEffects, eta1, eta2, sep = '\n')
  }

  manifestEffects <- "#definition of latent variables"
  for(i in burning:(tp-1)){
    eta1 <- paste("eta1_t", i, " =~ 1*y1_t", i ,sep = "")
    eta2 <- paste("eta2_t", i, " =~ 1*y2_t", i ,sep = "")
    manifestEffects <- paste(manifestEffects, eta1, eta2, sep = '\n')
  }

  latentDisturbances <- paste("#latent disturbances",
                              paste("eta1_t",burning," ~~ phi1_",burning,"*eta1_t",burning, "+phi12_",burning,"*eta2_t",burning, sep = ""),
                              paste("eta2_t",burning," ~~ phi2_",burning,"*eta2_t",burning, sep = ""),
                              sep = "\n")
  for(i in (burning+1):(tp-1)){
    eta1 <- paste("eta1_t", i, " ~~ phi1*eta1_t", i ,sep = "")
    eta2 <- paste("eta2_t", i, " ~~ phi2*eta2_t", i ,sep = "")
    latentDisturbances <- paste(latentDisturbances, eta1, eta2, sep = '\n')
  }

  manifestErrors <- "#manifest errors"
  for(i in burning:(tp-1)){
    y1 <- paste("y1_t", i, " ~~ 0*y1_t", i ,sep = "")
    y2 <- paste("y2_t", i, " ~~ 0*y2_t", i ,sep = "")
    manifestErrors <- paste(manifestErrors, y1, y2, sep = '\n')
  }

  AnalysisModel <- paste(manifestEffects, latentEffects, latentDisturbances, manifestErrors, sep = '\n')


  ret = list("SimulationModel" = SimulationModel, "AnalysisModel" = AnalysisModel, "raw_data" = raw_data, "sample_size" = sample_size)
  return(ret)
}
