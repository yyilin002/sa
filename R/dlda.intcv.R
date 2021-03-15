#' Diagonal Linear Discriminant Classifier
#'
#' Build a Diagonal Linear Discriminant classifier.
#'
#' @references https://cran.r-project.org/web/packages/HiDimDA/index.html
#' @param X dataset to be trained. This dataset must have rows as probes and columns as samples.
#' @param y a vector of sample group of each sample for the dataset to be trained.
#' It must have an equal length to the number of samples in \code{X}.
#' @param kfold placeholder with no meaning, default as NULL.
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @return a list of 4 elements:
#' \item{mc}{an internal misclassification error rate}
#' \item{time}{the processing time}
#' \item{model}{a DLDA classifier}
#' @export
#' @export dlda.intcv
#' @import HiDimDA
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
#' dlda.int <- dlda.intcv(X = biological.effect.nc.tr,
#'                      y = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                      kfold = NULL, seed = 1)
#'




"dlda.intcv" <- function(kfold = NULL, X, y, seed){
  ptm <- proc.time()
  set.seed(seed)

  dlda <- Dlda(data = t(X), grouping = factor(y), ldafun = "classification")
  pred <- predict(dlda, t(X), grpcodes = levels(factor(y)))$class
  mc <- tabulate.ext.err.func(pred, y)

  time <- proc.time() - ptm
  return(list(mc = mc, time = time, model = dlda, cfs = NULL))
}




