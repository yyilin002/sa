#' switch.norm.funcs
#'
#' @title Switch Normalization Functions
#' @description Switch the normalization method name into function name, for running in the normalization procedure
#' @param norm.list Switch all the build-in normalization methods into function names, including "NN", "MN", "QN", and "VSN".
#' @param norm.funcs New functions that user can create by themselves.
#' @return A list that transforms the normalization method name into function name.
#' @keywords switch names
#' @export
#' @examples
#' switch.norm.funcs(norm.list = c("NN", "MN", "QN"), norm.funcs = NULL)
#'


"switch.norm.funcs" <- function(norm.list = c("NN", "MN", "QN", "VSN"),
                                 norm.funcs = NULL){
  temp <- tolower(norm.list) %in% c("nn", "mn", "qn", "vsn")
  stopifnot(length(unique(norm.list[!temp])) == length(norm.funcs))

  new.norm.funcs <- tolower(norm.list)
  new.norm.funcs <- gsub("nn", "nn.norm",
                         gsub("mn", "med.norm",
                              gsub("qn", "quant.norm",
                                   gsub("vsn", "vs.norm",
                                        tolower(new.norm.funcs)))))
  new.norm.funcs <- c(new.norm.funcs[temp], norm.funcs)
  return(new.norm.funcs)
}
