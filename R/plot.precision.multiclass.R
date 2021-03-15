#' plot.precision.multiclass
#'
#' Plot the average external error rate of different classifications from the output of \code{precision.simulate.multiclass}
#'
#' @param data output variable of \code{precision.simulate.multiclass}
#' @param mytitle plot title
#' @param class.order the order of classifications
#' @param ylim the length of y axis in the plot
#' @param save.name plot name to be saved if not NULL
#'
#' @export plot.precision.multiclass
#' @keywords result.manage
#' @examples
#' \dontrun{
#' set.seed(101)
#' biological.effect <- estimate.biological.effect(uhdata = uhdata.pl)
#' handling.effect <- estimate.handling.effect(uhdata = uhdata.pl,
#'                              nuhdata = nuhdata.pl)
#'
#' ctrl.genes <- unique(rownames(uhdata.pl))[grep("NC", unique(rownames(uhdata.pl)))]
#'
#' biological.effect.nc <- biological.effect[!rownames(biological.effect) %in% ctrl.genes, ]
#' handling.effect.nc <- handling.effect[!rownames(handling.effect) %in% ctrl.genes, ]
#'
#' group.id <- substr(colnames(biological.effect.nc), 7, 7)
#'
#' # randomly split biological effect data into training and test set with
#' # equal number of endometrial and ovarian samples
#' biological.effect.train.ind <- colnames(biological.effect.nc)[c(sample(which(group.id == "E"), size = 64),
#'                                           sample(which(group.id == "V"), size = 64))]
#' biological.effect.test.ind <- colnames(biological.effect.nc)[!colnames(biological.effect.nc) %in% biological.effect.train.ind]
#' biological.effect.train.test.split =
#'   list("tr" = biological.effect.train.ind,
#'        "te" = biological.effect.test.ind)
#'
#' # non-randomly split handling effect data into training and test set
#' handling.effect.train.test.split =
#'   list("tr" = c(1:64, 129:192),
#'        "te" = 65:128)
#'
#' biological.effect.nc.tr <- biological.effect.nc[, biological.effect.train.ind]
#' biological.effect.nc.te <- biological.effect.nc[, biological.effect.test.ind]
#' handling.effect.nc.tr <- handling.effect.nc[, c(1:64, 129:192)]
#' handling.effect.nc.te <- handling.effect.nc[, 65:128]
#'
#' # Simulation
#' precision.multiclass.results = precision.simulate.multiclass(seed = 0, N = 3,
#'                                             biological.effect.tr = biological.effect.nc.tr,
#'                                             biological.effect.te = biological.effect.nc.te,
#'                                             handling.effect.tr = handling.effect.nc.tr,
#'                                             handling.effect.te = handling.effect.nc.te,
#'                                             group.id.tr = substr(colnames(biological.effect.nc.tr), 7, 7),
#'                                             group.id.te = substr(colnames(biological.effect.nc.te), 7, 7),
#'                                             train.design.met = "BLK",
#'                                             test.design.met = "STR",
#'                                             train.norm.met = "MN",
#'                                             test.norm.met = "fMN",
#'                                             class.list = c("SVM", "kNN", "LASSO"),
#'                                             train.batch.id = list(1:40, 41:64, (129:152)-64, (153:192)-64),
#'                                             test.batch.id = list((65:80)-64,(81:114)-64,(115:128)-64))
#'
#' # Plot
#' plot.precision.multiclass(data = precision.multiclass.results,
#'                           mytitle = "Average external error rates",
#'                           class.order = c("SVM", "kNN", "ClaNC"),
#'                           ylim = c(0,0.5),
#'                           save.name = "myimage")
#' }
#'






plot.precision.multiclass <- function(data,
                                 mytitle = "Average external error rates",
                                 class.order = NULL,
                                 ylim = c(0,0.5),
                                 save.name = NULL){

  errors <- c()
  for(i in 1:length(data$error_store)){
    errors[i] <- mean(data$error_store[[i]][[2]])
  }
  plot <- barplot(errors, names.arg = class.order, col = "blue", main = mytitle, ylim = ylim)

  if(!is.null(save.name)){
    dev.copy(png, paste0(save.name, ".png"))
    dev.off()
  }
}
