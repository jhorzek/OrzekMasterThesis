#' compute the Root Mean Squared Error for parameters in A and S
#'
#' @author Jannik Orzek
#' @import OpenMx
#'
#' @examples
#'
#' @export
computeRMSE <- function(A_est,S_est, Apop, Spop, Afree, Sfree){
  RMSE_Amat <- sqrt(mean((Apop[Afree] - A_est[Afree])^2))
  RMSE_Smat <- sqrt(mean((Spop[Sfree] - S_est[Sfree])^2))
  RMSE_comb <- sqrt(mean(c((Apop[Afree] - A_est[Afree])^2,(Spop[Sfree] - S_est[Sfree])^2)))

  ret = matrix(NA, nrow= 3, ncol = 1)
  rownames(ret) <- c("RMSE_Amat", "RMSE_Smat", "RMSE_comb")
  ret["RMSE_Amat",1] <- RMSE_Amat
  ret["RMSE_Smat",1] <- RMSE_Smat
  ret["RMSE_comb",1] <- RMSE_comb
  return(ret)
}
