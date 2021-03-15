#' Support Vector Machine  Classifier
#'
#' Build a support vactor machine classifier using internal cross validation to choose the turning parameter,
#' with a 5-fold cross validation as default.
#'
#' @references https://cran.r-project.org/web/packages/e1071/index.html
#' @param X dataset to be trained. This dataset must have rows as probes and columns as samples.
#' @param y a vector of sample group of each sample for the dataset to be trained.
#' It must have an equal length to the number of samples in \code{X}.
#' @param kfold number of folds. By default, \code{kfold = 5}.
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @return a list of 4 elements:
#' \item{mc}{an internal misclassification error rate}
#' \item{time}{the processing time of performing internal validation with SVM}
#' \item{model}{a SVM classifier, resulted from \code{cv.fit}}
#' @export svm.intcv
#' @import e1071
#' @keywords classification
#' @examples
#' set.seed(101)
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#' biological.effect.nc <- biological.effect[!rownames(biological.effect)
#'   %in% ctrl.genes, ]
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(
#'   group.id == "E"), size = 64),
#'   sample(which(group.id == "V"), size = 64))]
#' biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#'
#' svm.int <- svm.intcv(X = biological.effect.nc.tr,
#'                      y = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)
#'

"svm.intcv" <- function(kfold = 5, X, y, seed){
  ptm <- proc.time()
  set.seed(seed)

  svm_tune = tune.svm(x = data.matrix(t(X)), y = factor(y), tunecontrol = tune.control(cross = kfold))

  time <- proc.time() - ptm
  return(list(mc = svm_tune$best.performance, time = time, model = svm_tune$best.model))
}




