#' create.storage
#'
#' @title Create Storage for Output
#' @description During the implement of \code{precision.simulate} function, according to the length of \code{design.list},
#' \code{class.list}, and \code{norm.list}, to create a storage list to store the output from \code{precision.simulate}.
#' @param design.list List of design methods to be applied in \code{precision.simulate}.
#' @param class.list List of classification methods to be applied in \code{precision.simulate}
#' @param norm.list List of normalization methods to be applied in \code{precision.simulate}
#' @param validating.sets Types of validation methods.
#' @return A list whose size is determined by the input parameters.
#' @keywords storage
#' @export
#' @examples
#' create.storage(design.list = c("CC+", "CC-", "PC+", "PC-"),
#'                class.list = c("PAM", "LASSO"),
#'                norm.list = c("NN", "QN"),
#'                validating.sets = c("test", "gold.standard"))
#'
#'


"create.storage" <- function(design.list = c("CC+", "CC-", "PC+", "PC-"),
                             class.list = c("PAM", "LASSO"),
                             norm.list = c("NN", "QN"),
                             validating.sets = c("test", "gold.standard")){
  n.design <- length(design.list)
  n.norm <- length(norm.list)
  n.class <- length(class.list)
  n.valid <- length(validating.sets)

  storage <- matrix(rep(list(rep(list(NA), n.norm*n.class)), n.design*n.valid),
                    ncol = n.design, nrow = n.valid)
  colnames(storage) <- design.list
  rownames(storage) <- validating.sets
  for(i in 1:nrow(storage)){
    for(j in 1:ncol(storage)){
      names(storage[i, j][[1]]) <- paste(rep(class.list, each = n.norm), norm.list, sep = ".")
    }
  }

  return(storage)
}
