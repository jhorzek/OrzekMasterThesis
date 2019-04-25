#' getQdeltaT
#'
#'
#' @author Jannik Orzek
#' @import OpenMx
#'
#' @examples
#'
#' @export
#'
#'

getQdeltaT <- function(DRIFT, DIFFUSION, deltaT){
  Qmat <- DIFFUSION%*%t(DIFFUSION)
  Asharp <- kronecker(DRIFT,diag(ncol(DRIFT))) + kronecker(diag(ncol(DRIFT)), DRIFT)
  QdeltaT <- matrix(solve(Asharp) %*% (expm(Asharp*deltaT)-diag(ncol(expm(Asharp*deltaT)))) %*% cvectorize(Qmat), ncol= ncol(DRIFT), nrow = nrow(DRIFT), byrow = T)
  return(QdeltaT)
  }

