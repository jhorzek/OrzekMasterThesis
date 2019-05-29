#' createARCLModel
#'
#'
#' createARCLModel creates an ARCL(1)-SEM model with OpenMx with >= 1 latent variables and 1 manifest per latent and time point for data simulation. It returns the covariance based and a raw data based model
#'
#' @param numLatent number of latent varaibles per time point
#' @param Timepoints number of time points observed in the final data set
#' @param burning number of observations to burn before the initial time point
#' @param sampleSize sample size
#' @param Avalues values in ARCL matrix
#' @param Afree free parameters in ARCL matrix
#' @param Alabel labels parameters for ARCL matrix
#' @param Svalues values in S matrix
#' @param Sfree free parameters in S matrix
#' @param Slabel labels parameters for S matrix
#' @param S_firstObsAllFree should all initial observations be allowed to correlate?
#'
#' @author Jannik Orzek
#' @import OpenMx
#'
#' @examples
#'
#' library(OpenMx)
#'
#' Alabel <- matrix(c("a11", "a12", "a13", "a14", "a15", "a21", "a22", "a23", "a24", "a25", "a31", "a32", "a33", "a34", "a35", "a41", "a42", "a43", "a44", "a45", "a51", "a52", "a53", "a54", "a55"), nrow = 5, byrow = T)
#'
#' Avalues <- diag(.5,nrow=5,ncol = 5)
#' Avalues[2:5,1] <-.3
#' Avalues[3:5,2] <-.3
#' Avalues[4,5] <-.3
#' Afree <- matrix(TRUE, 5,5)
#' Alabel <- matrix(c("a11", "a12", "a13", "a14", "a15", "a21", "a22", "a23", "a24", "a25", "a31", "a32", "a33", "a34", "a35", "a41", "a42", "a43", "a44", "a45", "a51", "a52", "a53", "a54", "a55"), nrow = 5, byrow = T)
#'
#' Svalues <- diag(.75,nrow=5,ncol = 5)
#' Sfree <- diag(TRUE, 5,5)
#' Slabel <- matrix(c("s11", NA, NA, NA, NA,
#'                    NA, "s22", NA, NA, NA,
#'                    NA, NA, "s33", NA, NA,
#'                    NA, NA, NA, "s44", NA,
#'                    NA, NA, NA, NA, "s55"), nrow = 5, byrow = T)
#'
#'
#' Svalues_init <- diag(1,nrow=5,ncol = 5)
#' Sfree_init <- diag(TRUE, 5,5)
#' Slabel_init <- matrix(NA, nrow = 5, byrow = T)
#'
#' simObj <- simARCL(numLatent = 5, Timepoints = 4, burning = 10,
#'                   Avalues = Avalues, Alabel = Alabel, Afree = Afree,
#'                   Svalues = Svalues, Sfree = Sfree, Slabel = Slabel,
#'                   Svalues_init = Svalues_init, Sfree_init = Sfree_init, Slabel_init = Slabel_init,
#'                   sampleSize = 1000, S_firstObsAllFree = T)
#'
#'
#' temp <- mxRun(simObj$mxARCL_cov)
#' summary(temp)
#'
#' temp2 <- mxTryHard(simObj$mxARCL_FIML)
#' summary(temp2)
#'
#' @export
#'
#'
createARCLModel <- function(numLatent, Timepoints, burning, Avalues, Afree, Alabel, Svalues, Sfree, Slabel, S_firstObsAllFree = TRUE, SimulatedDataSet){

  simTimepoints = Timepoints + burning

  #### Analysis model
  # generate variable names:
  OBS_manifestVarNames <- c()
  OBS_latentVarNames <- c()
  for(i in (burning+1):simTimepoints){
    OBS_manifestVarNames <- c(OBS_manifestVarNames, paste("x", 1:numLatent,"_t",i, sep=""))
    OBS_latentVarNames <- c(OBS_latentVarNames, paste("eta", 1:numLatent,"_t",i, sep=""))
  }

  OBS_Amat_values <- matrix(0,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  OBS_Amat_free <- matrix(FALSE,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  OBS_Amat_label <- matrix(NA,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))

  rownames(OBS_Amat_values) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(OBS_Amat_values) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  rownames(OBS_Amat_free) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(OBS_Amat_free) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  rownames(OBS_Amat_label) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(OBS_Amat_label) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  # loadings of latent on manifest:
  OBS_Amat_values[1:length(OBS_manifestVarNames),(length(OBS_manifestVarNames)+1):(length(OBS_manifestVarNames)+length(OBS_latentVarNames))] <- diag(length(OBS_manifestVarNames))

  # regressions between latent:
  for (row in seq(which(rownames(OBS_Amat_values) == paste("eta1_t",burning+2,sep = "")), which(rownames(OBS_Amat_values) == paste("eta1_t", simTimepoints, sep = "")), by = numLatent)){

    OBS_Amat_values[row:(row+numLatent-1),(row-numLatent):(row-1)] <- Avalues
    OBS_Amat_free[row:(row+numLatent-1),(row-numLatent):(row-1)] <- Afree
    OBS_Amat_label[row:(row+numLatent-1),(row-numLatent):(row-1)] <- Alabel
  }

  # (Co)variances:
  OBS_Smat_values <- matrix(0,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  OBS_Smat_free <- matrix(FALSE,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  OBS_Smat_label <- matrix(NA,nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol =length(OBS_manifestVarNames)+length(OBS_latentVarNames))

  rownames(OBS_Smat_values) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(OBS_Smat_values) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  rownames(OBS_Smat_free) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(OBS_Smat_free) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  rownames(OBS_Smat_label) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  colnames(OBS_Smat_label) <- c(OBS_manifestVarNames, OBS_latentVarNames)

  # variances of manifest variables:
  OBS_Smat_values[1:length(OBS_manifestVarNames),1:length(OBS_manifestVarNames)] <- 0 # zero variance

  # variances of latent variables:
  # initial timepoint:
  # set (co)variances of initial observations free

  if(S_firstObsAllFree){
    row = which(rownames(OBS_Smat_values) == paste("eta1_t", burning+1 ,sep=""))
    OBS_Smat_values[row:(row+numLatent-1),row:(row+numLatent-1)] <- 0
    diag(OBS_Smat_values[row:(row+numLatent-1),row:(row+numLatent-1)]) <- 1
    OBS_Smat_free[row:(row+numLatent-1),row:(row+numLatent-1)] <- TRUE
    OBS_Slabel_first <- matrix(NA, numLatent,numLatent)
    for(row in 1:nrow(OBS_Slabel_first)){
      OBS_Slabel_first[row,] <- paste("s", row, 1:ncol(OBS_Slabel_first), "_first", sep = "")
    }
    # from https://stackoverflow.com/questions/33026183/r-make-symmetric-matrix-from-lower-diagonal
    makeSymm <- function(m) {
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
      return(m)
    }
    OBS_Slabel_first = makeSymm(OBS_Slabel_first)
    row = which(rownames(OBS_Smat_values) == paste("eta1_t", burning+1 ,sep=""))
  }else{
    # only variances free:
    row = which(rownames(OBS_Smat_values) == paste("eta1_t", burning+1 ,sep=""))
    OBS_Smat_values[row:(row+numLatent-1),row:(row+numLatent-1)] <- diag(1, numLatent,numLatent)
    OBS_Smat_free[row:(row+numLatent-1),row:(row+numLatent-1)] <- diag(TRUE, numLatent,numLatent)
    OBS_Slabel_first <- matrix(NA, numLatent,numLatent)
    temp <- paste(paste("s", 1:numLatent, "_first", sep = ""))
    diag(OBS_Slabel_first) <- temp
  }
  OBS_Smat_label[row:(row+numLatent-1),row:(row+numLatent-1)] <- OBS_Slabel_first

  # all other:
  for (row in seq(which(rownames(OBS_Smat_values) == paste("eta1_t", burning+2 ,sep="")), which(rownames(OBS_Smat_values) == paste("eta1_t", simTimepoints, sep = "")), by = numLatent)){

    OBS_Smat_values[row:(row+numLatent-1),row:(row+numLatent-1)] <- Svalues
    OBS_Smat_free[row:(row+numLatent-1),row:(row+numLatent-1)] <- Sfree
    OBS_Smat_label[row:(row+numLatent-1),row:(row+numLatent-1)] <- Slabel
  }

  # Fmatrix
  OBS_Fmat_values = matrix(0, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames), nrow = length(OBS_manifestVarNames))
  OBS_Fmat_free = matrix(FALSE, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames), nrow = length(OBS_manifestVarNames))
  OBS_Fmat_label = matrix(NA, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames), nrow = length(OBS_manifestVarNames))

  rownames(OBS_Fmat_values) <- c(OBS_manifestVarNames)
  colnames(OBS_Fmat_values) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  rownames(OBS_Fmat_free) <- c(OBS_manifestVarNames)
  colnames(OBS_Fmat_free) <- c(OBS_manifestVarNames, OBS_latentVarNames)
  rownames(OBS_Fmat_label) <- c(OBS_manifestVarNames)
  colnames(OBS_Fmat_label) <- c(OBS_manifestVarNames, OBS_latentVarNames)

  diag(OBS_Fmat_values[1:length(OBS_manifestVarNames), 1:length(OBS_manifestVarNames)]) = 1

  # M matrix
  ##### MMatrix #####
  OBS_Mmat_values = matrix(0, nrow = 1, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  OBS_Mmat_free = matrix(FALSE, nrow = 1, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames))
  OBS_Mmat_label = matrix(NA, nrow = 1, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames))

  # build mx matrices:
  OBS_A <- mxMatrix(type = "Full", nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames),
                    values = OBS_Amat_values, free = OBS_Amat_free, labels = OBS_Amat_label, name = "A")
  OBS_S <- mxMatrix(type = "Full", nrow = length(OBS_manifestVarNames)+length(OBS_latentVarNames), ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames),
                    values = OBS_Smat_values, free = OBS_Smat_free, labels = OBS_Smat_label, name = "S")
  OBS_F <- mxMatrix(type = "Full", nrow = length(OBS_manifestVarNames), ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames),
                    values = OBS_Fmat_values, free = OBS_Fmat_free, labels = OBS_Fmat_label, name = "F")
  OBS_M <- mxMatrix(type = "Full", nrow = 1, ncol = length(OBS_manifestVarNames)+length(OBS_latentVarNames),
                    values = OBS_Mmat_values, free = OBS_Mmat_free, labels = OBS_Mmat_label, name = "M",dimnames = list("NA",c(OBS_manifestVarNames, OBS_latentVarNames)))


  ##### fit function #####

  expect_cov <- mxExpectationRAM(A="A", S="S", F="F")

  expect_FIML <- mxExpectationRAM(A="A", S="S", F="F", M = "M")



    mxARCL_cov <- mxModel(model = "ARCL",
                          mxData(observed = cov(SimulatedDataSet[,OBS_manifestVarNames]),type = "cov", numObs =  sampleSize),
                          OBS_A, OBS_S, OBS_F,
                          #manifestVars = OBS_manifestVarNames,
                          #latentVars = OBS_latentVarNames,
                          expect_cov, mxFitFunctionML())
    mxARCL_FIML <- mxModel(model = "ARCL",
                           mxData(observed = SimulatedDataSet[,OBS_manifestVarNames],type = "raw"),
                           OBS_A, OBS_S, OBS_F, OBS_M,
                           manifestVars = OBS_manifestVarNames,
                           latentVars = OBS_latentVarNames,
                           expect_FIML, mxFitFunctionML())


  ret = list("mxARCL_cov" = mxARCL_cov, "mxARCL_FIML" = mxARCL_FIML)

  return(ret)
}
