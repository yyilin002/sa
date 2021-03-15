#' switch.classifier.funcs.class
#'
#' @title Switch classfication functions
#' @description Transform classfication method names into a fitting function name and a prediction function name,
#' for using in the classfication procedure
#' @param class.list A list of classification method name to be transformed. The build-in methods are "PAM", "LASSO",
#' "ClaNC", "ranFor", "SVM", "kNN" and "DLDA".
#' @param class.funcs Other classfication functions that users can create by themselves.
#' @param pred.funcs Other prediction funcstions that users can create by themselves.
#' @return A list containing \code{build.funcs}, which means the fitting function names, and \code{pred.funcs},
#' which means the prediction function names. The fitting functions is to be fitted on the training set and prediction
#' function is to be used on the test set.
#' @keywords switch names
#' @export
#' @examples
#' switch.classifier.funcs.class(class.list = c("PAM", "LASSO", "ClaNC", "ranFor", "SVM"), class.funcs = NULL, pred.funcs = NULL)
#'



"switch.classifier.funcs.class" <- function(class.list = c("PAM", "LASSO", "ClaNC", "ranFor", "SVM", "kNN", "DLDA", "custom"),
                                             class.funcs = NULL,
                                             pred.funcs = NULL){
  temp <- tolower(class.list) %in% c("pam", "lasso", "clanc", "ranfor", "svm", "knn", "dlda", "custom")
  stopifnot(length(class.list[!temp]) == length(class.funcs))
  stopifnot(length(class.list[!temp]) == length(pred.funcs))

  new.class.funcs <- tolower(class.list)
  new.class.funcs <- gsub("pam", "pam.intcv",
                          gsub("lasso", "lasso.intcv",
                               gsub("clanc", "clanc.intcv",
                                    gsub("ranfor", "ranfor.intcv",
                                         gsub("svm", "svm.intcv",
                                              gsub("knn", "knn.intcv",
                                                   gsub("dlda", "dlda.intcv",
                                                        gsub("custom", "custom.intcv",
                                                             tolower(new.class.funcs)))))))))

  new.class.funcs <- c(new.class.funcs[temp], class.funcs)

  new.pred.funcs <- tolower(class.list)
  new.pred.funcs <- gsub("pam", "pam.predict",
                         gsub("lasso", "lasso.predict",
                              gsub("clanc", "clanc.predict",
                                   gsub("ranfor", "ranfor.predict",
                                        gsub("svm", "svm.predict",
                                             gsub("knn", "knn.predict",
                                                  gsub("dlda", "dlda.predict",
                                                       gsub("custom", "custom.predict",
                                                            tolower(new.pred.funcs)))))))))

  new.pred.funcs <- c(new.pred.funcs[temp], pred.funcs)


  return(list(build.funcs = new.class.funcs,
              pred.funcs = new.pred.funcs))
}
