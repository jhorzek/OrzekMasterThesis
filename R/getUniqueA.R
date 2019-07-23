#' getUniqueA
#'
#' Extracts the unique elements from an mxModel A-matrix. Used for ARCL models with equality constraints
#'
#' @param ARCL_AMatrix A matrix from mxModel (full mxMatrix, not just values)
#'
#' @author Jannik Orzek
#' @import OpenMx
#'
#' @examples
#'
#' @export
#'
#'

getUniqueA <- function(ARCL_AMatrix){
  values <- ARCL_AMatrix$values[ARCL_AMatrix$free]
  labels <- ARCL_AMatrix$labels[ARCL_AMatrix$free]
  unique <- unique(labels)

  ret_values <- vector()
  ret_labels <- vector()
  for(i in 1:length(labels)){
    if(labels[i] %in% unique){
      value <- values[i]
      label <- labels[i]
      ret_values <- c(ret_values, value)
      ret_labels <- c(ret_labels, label)
      unique <- unique[!unique == labels[i]]
    }

  }
  ret <- matrix(ret_values, nrow = 1)
  colnames(ret) <- ret_labels
  return(ret)
}
