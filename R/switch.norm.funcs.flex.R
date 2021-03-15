#' switch.norm.funcs.flex
#'
#' @title Switch normalization funcsions in a flexible way
#' @description Transform the normalization method names into functions names, for running in the normalization procedure
#' @param norm.list Switch all the build-in normalization methods into function names, including "NN", "MN", "QN", and "VSN".
#' @param norm.funcs New functions that user can create by themselves.
#' @return A list that transforms the normalization method name into function name.
#' @keywords switch names
#' @export
#' @examples
#' switch.norm.funcs.flex(norm.list = c("NN", "MN", "QN"), norm.funcs = NULL)
#'


"switch.norm.funcs.flex" <- function(norm.list = c("NN", "QN"),
                                norm.funcs = NULL){

  temp <- tolower(norm.list) %in% c("nn", "mn", "qn", "vsn")
  stopifnot(length(unique(norm.list[!temp])) == length(norm.funcs))

  new.norm.funcs <- tolower(norm.list)

  switch_x <- sapply(new.norm.funcs, function(x) switch(x,
                                                        "nn" = "nn.norm",
                                                        "mn" = "med.norm",
                                                        "qn" = "quant.norm",
                                                        "vsn" = "vs.norm"))

  switch_x <- as.character(switch_x)
  switch_x[which(switch_x == "NULL")] <-
    paste0(tolower(norm.list)[which(switch_x == "NULL")], ".norm")
  new.norm.funcs <- switch_x

  stopifnot(new.norm.funcs %in% c("nn.norm", "med.norm", "quant.norm", "vs.norm",
                                  norm.funcs))

  return(new.norm.funcs)
}
