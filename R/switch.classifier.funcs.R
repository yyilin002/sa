#' switch.classifier.funcs
#'
#' @title Switch classfication functions
#' @description Transform classfication method names into a fitting function name and a prediction function name,
#' for using in the classfication procedure
#' @param class.list A list of classification method name to be transformed. The build-in methods are "PAM" and "LASSO".
#' @param class.funcs Other classfication functions that users can create by themselves.
#' @param pred.funcs Other prediction funcstions that users can create by themselves.
#' @return A list containing \code{build.funcs}, which means the fitting function names, and \code{pred.funcs},
#' which means the prediction function names. The fitting functions is to be fitted on the training set and prediction
#' function is to be used on the test set.
#' @keywords switch names
#' @export
#' @examples
#' switch.classifier.funcs.class(class.list =  c("PAM", "LASSO"),
#'                               class.funcs = NULL, pred.funcs = NULL)
#'
#'

"switch.classifier.funcs" <- function(class.list = c("PAM", "LASSO"),
                                      class.funcs = NULL,
                                      pred.funcs = NULL){
  temp <- tolower(class.list) %in% c("pam", "lasso")
  stopifnot(length(class.list[!temp]) == length(class.funcs))
  stopifnot(length(class.list[!temp]) == length(pred.funcs))

  new.class.funcs <- tolower(class.list)
  new.class.funcs <- gsub("pam", "pam.intcv",
                          gsub("lasso", "lasso.intcv", tolower(new.class.funcs)))
  new.class.funcs <- c(new.class.funcs[temp], class.funcs)

  new.pred.funcs <- tolower(class.list)
  new.pred.funcs <- gsub("pam", "pam.predict",
                         gsub("lasso", "lasso.predict", tolower(new.pred.funcs)))
  new.pred.funcs <- c(new.pred.funcs[temp], pred.funcs)

  return(list(build.funcs = new.class.funcs,
              pred.funcs = new.pred.funcs))
}
