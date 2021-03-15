#' Tabulate.ext.err.func
#'
#' An assistant function used to calculate the error rate of classification.
#'
#' @param pred.obj predicted group classes
#' @param obs.grp actual group classes
#' @return error rate of predicted classification
#' @keywords classification
#' @export tabulate.ext.err.func
#' @examples
#' pred.obj = c("V", "V", "V", "V", "E", "V", "E", "E", "E", "E")
#' obs.grp = c("V", "V", "V", "V", "V", "E", "E", "E", "E", "E")
#' tabulate.ext.err.func(pred.obj, obs.grp)
#'



"tabulate.ext.err.func" <- function(pred.obj, obs.grp)
  return(1 - sum(diag(table(pred.obj, obs.grp)))/length(obs.grp))
