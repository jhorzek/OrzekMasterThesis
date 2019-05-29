#' simARCL_scaled_10
#'
#'
#' simARCL_scaled_10 creates an ARCL(1)-SEM model with OpenMx with 5 latent variables and 1 manifest per latent and time point for data simulation. It returns the simulated data set and a covariance based and a raw data based model. 10 Timepoints are created of which the initial five are omitted (burned)
#' The model has 10% non-zero cross-lagged effects
#'
#'
#' @author Jannik Orzek
#' @import OpenMx
#'
#' @export
#'
#'
simARCL_scaled_10 <- function(seed, sampleSize, autoEffect, crossEffect){
  set.seed(seed)

  x1_t1 <-rnorm(n = sampleSize, mean = 0, sd = 1)
  x2_t1 <-rnorm(n = sampleSize, mean = 0, sd = 1)
  x3_t1 <-rnorm(n = sampleSize, mean = 0, sd = 1)
  x4_t1 <-rnorm(n = sampleSize, mean = 0, sd = 1)
  x5_t1 <-rnorm(n = sampleSize, mean = 0, sd = 1)

  x1_t2 <- autoEffect*x1_t1 + crossEffect*x2_t1 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t1) + crossEffect^2*var(x2_t1)+ autoEffect*crossEffect*2*cov(x1_t1,x2_t1))))
  x2_t2 <- autoEffect*x2_t1 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t1))))
  x3_t2 <- autoEffect*x3_t1 + crossEffect*x4_t1 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t1) + crossEffect^2*var(x4_t1)+ autoEffect*crossEffect*2*cov(x3_t1,x4_t1))))
  x4_t2 <- autoEffect*x4_t1 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t1))))
  x5_t2 <- autoEffect*x5_t1 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t1))))

  x1_t3 <- autoEffect*x1_t2 + crossEffect*x2_t2 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t2) + crossEffect^2*var(x2_t2)+ autoEffect*crossEffect*2*cov(x1_t2,x2_t2))))
  x2_t3 <- autoEffect*x2_t2 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t2))))
  x3_t3 <- autoEffect*x3_t2 + crossEffect*x4_t2 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t2) + crossEffect^2*var(x4_t2)+ autoEffect*crossEffect*2*cov(x3_t2,x4_t2))))
  x4_t3 <- autoEffect*x4_t2 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t2))))
  x5_t3 <- autoEffect*x5_t2 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t2))))

  x1_t4 <- autoEffect*x1_t3 + crossEffect*x2_t3 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t3) + crossEffect^2*var(x2_t3)+ autoEffect*crossEffect*2*cov(x1_t3,x2_t3))))
  x2_t4 <- autoEffect*x2_t3 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t3))))
  x3_t4 <- autoEffect*x3_t3 + crossEffect*x4_t3 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t3) + crossEffect^2*var(x4_t3)+ autoEffect*crossEffect*2*cov(x3_t3,x4_t3))))
  x4_t4 <- autoEffect*x4_t3 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t3))))
  x5_t4 <- autoEffect*x5_t3 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t3))))

  x1_t5 <- autoEffect*x1_t4 + crossEffect*x2_t4 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t4) + crossEffect^2*var(x2_t4)+ autoEffect*crossEffect*2*cov(x1_t4,x2_t4))))
  x2_t5 <- autoEffect*x2_t4 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t4))))
  x3_t5 <- autoEffect*x3_t4 + crossEffect*x4_t4 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t4) + crossEffect^2*var(x4_t4)+ autoEffect*crossEffect*2*cov(x3_t4,x4_t4))))
  x4_t5 <- autoEffect*x4_t4 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t4))))
  x5_t5 <- autoEffect*x5_t4 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t4))))

  x1_t6 <- autoEffect*x1_t5 + crossEffect*x2_t5 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t5) + crossEffect^2*var(x2_t5)+ autoEffect*crossEffect*2*cov(x1_t5,x2_t5))))
  x2_t6 <- autoEffect*x2_t5 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t5))))
  x3_t6 <- autoEffect*x3_t5 + crossEffect*x4_t5 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t5) + crossEffect^2*var(x4_t5)+ autoEffect*crossEffect*2*cov(x3_t5,x4_t5))))
  x4_t6 <- autoEffect*x4_t5 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t5))))
  x5_t6 <- autoEffect*x5_t5 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t5))))

  x1_t7 <- autoEffect*x1_t6 + crossEffect*x2_t6 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t6) + crossEffect^2*var(x2_t6)+ autoEffect*crossEffect*2*cov(x1_t6,x2_t6))))
  x2_t7 <- autoEffect*x2_t6 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t6))))
  x3_t7 <- autoEffect*x3_t6 + crossEffect*x4_t6 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t6) + crossEffect^2*var(x4_t6)+ autoEffect*crossEffect*2*cov(x3_t6,x4_t6))))
  x4_t7 <- autoEffect*x4_t6 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t6))))
  x5_t7 <- autoEffect*x5_t6 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t6))))

  x1_t8 <- autoEffect*x1_t7 + crossEffect*x2_t7 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t7) + crossEffect^2*var(x2_t7)+ autoEffect*crossEffect*2*cov(x1_t7,x2_t7))))
  x2_t8 <- autoEffect*x2_t7 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t7))))
  x3_t8 <- autoEffect*x3_t7 + crossEffect*x4_t7 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t7) + crossEffect^2*var(x4_t7)+ autoEffect*crossEffect*2*cov(x3_t7,x4_t7))))
  x4_t8 <- autoEffect*x4_t7 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t7))))
  x5_t8 <- autoEffect*x5_t7 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t7))))

  x1_t9 <- autoEffect*x1_t8 + crossEffect*x2_t8 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t8) + crossEffect^2*var(x2_t8)+ autoEffect*crossEffect*2*cov(x1_t8,x2_t8))))
  x2_t9 <- autoEffect*x2_t8 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t8))))
  x3_t9 <- autoEffect*x3_t8 + crossEffect*x4_t8 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t8) + crossEffect^2*var(x4_t8)+ autoEffect*crossEffect*2*cov(x3_t8,x4_t8))))
  x4_t9 <- autoEffect*x4_t8 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t8))))
  x5_t9 <- autoEffect*x5_t8 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t8))))

  x1_t10 <- autoEffect*x1_t9 + crossEffect*x2_t9 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x1_t9) + crossEffect^2*var(x2_t9)+ autoEffect*crossEffect*2*cov(x1_t9,x2_t9))))
  x2_t10 <- autoEffect*x2_t9 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x2_t9))))
  x3_t10 <- autoEffect*x3_t9 + crossEffect*x4_t9 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x3_t9) + crossEffect^2*var(x4_t9)+ autoEffect*crossEffect*2*cov(x3_t9,x4_t9))))
  x4_t10 <- autoEffect*x4_t9 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x4_t9))))
  x5_t10 <- autoEffect*x5_t9 + rnorm(n = sampleSize, 0, sd = sqrt(1-(autoEffect^2*var(x5_t9))))


  SimulatedDataSet <- cbind(x1_t1, x2_t1, x3_t1, x4_t1, x5_t1,
                            x1_t2, x2_t2, x3_t2, x4_t2, x5_t2,
                            x1_t3, x2_t3, x3_t3, x4_t3, x5_t3,
                            x1_t4, x2_t4, x3_t4, x4_t4, x5_t4,
                            x1_t5, x2_t5, x3_t5, x4_t5, x5_t5,
                            x1_t6, x2_t6, x3_t6, x4_t6, x5_t6,
                            x1_t7, x2_t7, x3_t7, x4_t7, x5_t7,
                            x1_t8, x2_t8, x3_t8, x4_t8, x5_t8,
                            x1_t9, x2_t9, x3_t9, x4_t9, x5_t9,
                            x1_t10, x2_t10, x3_t10, x4_t10, x5_t10)


  Avalues <- diag(.5,nrow=5,ncol = 5)
  Avalues[1,2] <-crossEffect
  Avalues[3,4] <-crossEffect

  Afree <- matrix(TRUE, 5,5)
  Alabel <- matrix(c("a11", "a12", "a13", "a14", "a15", "a21", "a22", "a23", "a24", "a25", "a31", "a32", "a33", "a34", "a35", "a41", "a42", "a43", "a44", "a45", "a51", "a52", "a53", "a54", "a55"), nrow = 5, byrow = T)

  Svalues <- diag(1-.5^2,nrow=5,ncol = 5)
  Svalues[1,1] <- 1-.5^2 -crossEffect^2
  Svalues[3,3] <- 1-.5^2 -crossEffect^2

  Sfree <- diag(TRUE, 5,5)
  Slabel <- matrix(c("s11", NA, NA, NA, NA,
                     NA, "s22", NA, NA, NA,
                     NA, NA, "s33", NA, NA,
                     NA, NA, NA, "s44", NA,
                     NA, NA, NA, NA, "s55"), nrow = 5, byrow = T)


  AnalysisModels <- createARCLModel(numLatent = 5, Timepoints = 5, burning = 5,
                                    Avalues = Avalues, Afree = Afree, Alabel = Alabel,
                                    Svalues = Svalues, Sfree = Sfree, Slabel = Slabel,
                                    S_firstObsAllFree = T, SimulatedDataSet = SimulatedDataSet
                                    )

  ret = list("mxARCL_cov" = AnalysisModels$mxARCL_cov, "mxARCL_FIML" = AnalysisModels$mxARCL_FIML, "SimulatedDataSet" = SimulatedDataSet)

  return(ret)
}
