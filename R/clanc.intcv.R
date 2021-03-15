#' Classification to Nearest Centroids  Classifier
#'
#' Build a Classification to Nearest Centroids classifier on the objective data.
#'
#' @references Alan R. Dabney, Author Notes.(2005) ClaNC: point-and-click software for classifying microarrays to nearest centroids,
#' https://academic.oup.com/bioinformatics/article/22/1/122/219377
#' @param X dataset to be trained. This dataset must have rows as probes and columns as samples.
#' @param y a vector of sample group of each sample for the dataset to be trained.
#' It must have an equal length to the number of samples in \code{X}.
#' @param kfold number of folds. By default, \code{kfold = 5}.
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @return a list of 4 elements:
#' \item{mc}{an internal misclassification error rate}
#' \item{time}{the processing time of performing internal validation with ClaNC}
#' \item{model}{a ClaNC classifier, resulted from \code{cv.fit}}
#' @export
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
#' clanc.int <- clanc.intcv(X = biological.effect.nc.tr,
#'                      y = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                      kfold = 5, seed = 1)
#'

"clanc.intcv" <- function(kfold = 5, X, y, seed){
    ptm <- proc.time()
    set.seed(seed)

    id = ifelse(y == "E", 1, 2)
    gene_names = rownames(X)
    class_names = unique(y)
    clanc_cv = cvClanc(X, id)

    active = min(which(clanc_cv$overallErrors == min(clanc_cv$overallErrors)))
    #plot(clanc_cv$overallErrors)

    train_out = trainClanc(X, id, gene_names)
    build_out = buildClanc(X, id, class_names, train_out, active = active)

    time <- proc.time() - ptm

    return(list(mc = build_out$overallError, time = time, model = build_out))
  }





