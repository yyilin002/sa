#' @title precision simulation with classification
#'
#' @description  Classification analysis applied on different normalization and study design.
#' Perform the simulation study in Qin et al. (see reference).
#'
#' @references Qin LX, Huang HC, Begg CB. Cautionary note on cross validation in molecular classification. Journal of Clinical Oncology. 2016
#'
#' @details The classification anlaysis of simulation study consists of the following main steps:
#'
#' First, The generation of training and test sets is the same as \code{precision.simultate}.
#' \code{precision.simulate.class} requires the training and test sets for both estimated biological effects and estimated handling effects.
#' The effects can be simulated as follows (using \code{estimate.biological.effect} and \code{estimate.handling.effect}).
#' The uniformly-handled dataset are used to approximate the biological effect for each sample,
#' and the difference between the two arrays (one from the uniformly-handled dataset and
#' the other from the nonuniformly-handled dataset, subtracting the former from the latter)
#' for the same sample are used to approximate the handling effect for each array in the nonuniformly-handled dataset.
#'
#' The samples are randomly split into a training set and a test set, balanced by tumor type (in Qin et al., training-to-test ratio is 2:1).
#' The arrays were then non-randomly split to a training set and a test set (in Qin et al., training set n = 128 -- the first 64 and last 64 arrays
#' in the order of array processing; test set n = 64 -- the middle 64 arrays).
#' This setup allows different pairings of arrays and samples by various different training-and-test-set splits.
#' Furthermore, biological signal strength and confounding level of the handling effects can be modified
#' (using \code{reduce.signal} and \code{amplify.handling.effect}).
#'
#' Second, apply "virtual re-hybridization" methods (using \code{rehybridize}) on the training and test sets.
#' There are 6 different methods to choose, besides doing no hybridization.
#' And it is also allowed to produce different methods on training set and test set.
#' This is specified in \code{train.design.met} and \code{test.design.met}.
#'
#' Third, apply normalization on the training and test sets. There are three normalization methods to choose on the training set,
#' median normalization, quantile normalization and variance stabilizing normalization. With these three methods,
#' test set can be processed the same method, or corresponding frozen normalization and pool normalization.
#' Besides, both sets can choose to do no normalization. The normalization methods are specified in \code{train.norm.met}
#' and \code{test.norm.met}. Data preprocessing and batch effects can be adjusted specified with \code{icombat}, \code{isva} and \code{iruv}.
#'
#' Fourth, choose a classfication method, specified with \code{class.met}. This will fit a classification on the training set,
#' with a 5-fold cross validation and predict on the test set. The internal and external validation misclassificaiton error estimation,
#' and the adjusted rand index on the test set will be both included in the output.
#'
#' For a given split of samples to training set versus test set,
#' \code{N} datasets will be simulated and analyzed for each array-assignment scheme.
#' For user-defined normalization method or classification method, please refer to the vignette.
#'
#' @param seed an integer used to initialize a pseudorandom number generator.
#' @param N number of simulation runs.
#' @param biological.effect.tr the training set of the estimated biological effects. This dataset must have rows as probes and columns as samples.
#' @param biological.effect.te the test set of the estimated biological effects. This dataset must have rows as probes and columns as samples.
#' It must have the same number of probes and the same probe names as the training set of the estimated biological effects.
#' @param handling.effect.tr the training set of the estimated handling effects. This dataset must have rows as probes and columns as samples.
#' It must have the same dimensions and the same probe names as the training set of the estimated biological effects.
#' @param handling.effect.te the test set of the estimated handling effects. This dataset must have rows as probes, columns as samples.
#' It must have the same dimensions and the same probe names as the training set of the estimated handling effects.
#' @param group.id.tr a vector of sample-group labels for each sample of the training set of the estimated biological effects.
#' It must be a 2-level non-numeric factor vector.
#' @param group.id.te a vector of sample-group labels for each sample of the test set of the estimated biological effects.
#' It must be a 2-level non-numeric factor vector.
#' @param train.design.met a string for study design to be applied on the training set.
#' The built-in designs are "NONE", CC+", "CC-", "PC+", "PC-", "BLK", and "STR" for "No Rehybridization", "Complete Confounding 1",
#' "Complete Confounding 2", "Partial Confounding 1", "Partial Confounding 2", "Blocking", and "Stratification" in Qin et al.
#' @param test.design.met a string for study design to be applied on the test set.
#' The built-in designs are "NONE", CC+", "CC-", "PC+", "PC-", "BLK", and "STR" for "No Rehybridization", "Complete Confounding 1",
#' "Complete Confounding 2", "Partial Confounding 1", "Partial Confounding 2", "Blocking", and "Stratification" in Qin et al.
#' @param train.norm.met a string for normalization method to be applied on the training set.
#' The build-in available normalization methods are "NN", "QN", "MN", "VSN" for "No Normalization", "Quantile Normalization",
#' "Median Normalization", and "Variance Stabilizing Normalization".
#' User can provide a list of normalization methods given the functions are supplied (also see norm.funcs).
#' @param test.norm.met a string for normalization method to be applied on the test set.
#' The build-in available normalization methods are "NN", "MN", "QN", "fMN", "fQN", "pMN", "pQN", "fVSN", for "No Normalization",
#' "Median Normalization", "Quantile Normalization", "Frozen Median Normalization", "Frozen Quantile Normalization",
#' "Pool Median Normalization", "Pool Quantile Normalization", "Frozen Vairance Stability Normalization".
#' User can provide a list of normalization methods given the functions are supplied (also see norm.funcs).
#' @param class.met a string for classfication method to be fitted on the training set and predict on test set.
#' The build-in available classfication methods are "PAM", "LASSO", "ClaNC", "ranFor", "SVM", "kNN" and "DLDA",
#' for "Prediction Analysis for Microarrays", "Least Absolute Shrinkage and Selection Operator", "Classification to Nearest Centroids",
#' "Random Forest", "Support Vector Machine", "K-Nearest Neighbors" and "Diagonal Linear Discriminant".
#' User can provide a list of classification methods given the correponding model-building and
#' predicting functions are supplied (also see \code{class.funcs} and \code{pred.funcs}). You can also use your own classification method
#' by setting \code{class.met} as "custom". The format to create \code{custom.intcv} and \code{custom.predict} please refers to other
#' classification methods in this package.
#' @param batch.id a list of array indices grouped by batches when data were profiled.
#' The length of the list must be equal to the number of batches in the data;
#' the number of array indices must be the same as the number of samples.
#' This is required if stratification study design is specified in \code{design.list}; otherwise \code{batch.id = NULL}.
#' @param icombat an indicator for combat adjustment. By default, \code{icombat = FALSE} for no ComBat adjustment.
#' @param isva an indicator for sva adjustment. By default, \code{isva = FALSE} for no sva adjustment.
#' @param iruv an indicator for RUV-4 adjustment. By default, \code{iruv = FALSE} for no RUV-4 adjustment.
#' @param biological.effect.tr.ctrl the training set of the negative-control probe biological effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same number of samples and the same sample names as \code{biological.effect.tr}.
#' @param handling.effect.tr.ctrl the training set of the negative-control probe handling effect data if \code{iruv = TRUE}.
#' This dataset must have rows as probes and columns as samples.
#' It also must have the same dimensions and the same probe names as \code{biological.effect.tr.ctrl}.
#' @param norm.funcs a list of strings for names of user-defined normalization method functions, in the order of \code{norm.list},
#' excluding any built-in normalization methods.
#' @param class.funcs a list of strings for names of user-defined classification model-building functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#' @param pred.funcs a list of strings for names of user-defined classification predicting functions, in the order of \code{class.list},
#' excluding any built-in classification methods.
#'
#' @return simulation study results -- a list of array-to-sample assignments, fitted models,
#' and misclassification error rates across simulation runs:
#' \item{assign_store}{array-to-sample assignments for the study design, classified by "Train" and "Test"}
#' \item{model_store}{models for the combination of study designs, normalization methods, and classification methods}
#' \item{error_store}{internal and external misclassification error rates for the combination of study designs,
#' normalization methods, and classification methods, classfied by "Train" and "Test"}
#' \item{ari_store}{adjusted rand index of the prediction on test data}
#'
#' @keywords simulation
#'
#' @import mclust
#' @export precision.simulate.class
#'
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
#'
#' redhalf.biological.effect.nc <- reduce.signal(biological.effect = biological.effect.nc, group.id = substr(colnames(biological.effect.nc), 7, 7),group.id.level = c("E", "V"),reduce.multiplier = 1/2)
#'
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
#' precision.result = precision.simulate.class(seed = 0, N = 3,
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
#'                                             class.met = "LASSO",
#'                                             train.batch.id = list(1:40, 41:64, (129:152)-64, (153:192)-64),
#'                                             test.batch.id = list((65:80)-64,(81:114)-64,(115:128)-64))
#'
#'
#'
#' }
#'



"precision.simulate.class" <- function(seed, N,
                                       biological.effect.tr, biological.effect.te,
                                       handling.effect.tr, handling.effect.te,
                                       group.id.tr, group.id.te,
                                       train.design.met = "NONE",
                                       test.design.met = "NONE",
                                       train.norm.met = "NN",
                                       test.norm.met = "NN",
                                       class.met = "LASSO",
                                       train.batch.id = NULL, test.batch.id = NULL,
                                       icombat = FALSE, isva = FALSE, iruv = FALSE,
                                       biological.effect.tr.ctrl = NULL,
                                       handling.effect.tr.ctrl = NULL,
                                       norm.funcs = NULL,
                                       class.funcs = NULL,
                                       pred.funcs = NULL){


  stopifnot(train.design.met %in% c("NONE", "CC+", "CC-", "PC+", "PC-", "BLK", "STR"))
  stopifnot(test.design.met %in% c("NONE", "CC+", "CC-", "PC+", "PC-", "BLK", "STR"))
  stopifnot(train.norm.met %in% c("NN", "MN", "QN", "VSN"))
  stopifnot(test.norm.met %in% c("NN", "MN", "QN", "fMN", "fQN", "pMN", "pQN", "VSN", "fVSN"))
  stopifnot(class.met %in% c("PAM", "LASSO", "ClaNC", "ranFor", "SVM", "kNN", "DLDA", "custom"))




  ### Initialize storage

  assign_store <- list()

  model_store <- list()

  error_store <- list()

  ari_store <- list() # adjusted random index


  for (k in 1:N) {# each of the N simulation

    cat(k, "round seed used:", seed + k, "\n")

    # Rehybridize on training set

    if(train.design.met == "NONE"){
      train.design.met2 <- tolower(train.design.met)
      eval(parse(text = paste0(train.design.met2, ".tr <- biological.effect.tr")))
    } else{
      if(train.design.met == "CC+") train.cc1.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "complete", rev.order = FALSE)
      if(train.design.met == "CC-") train.cc2.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "complete", rev.order = TRUE)
      if(train.design.met == "PC+") train.pc1.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "partial", rev.order = FALSE)
      if(train.design.met == "PC-") train.pc2.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.tr), degree = "partial", rev.order = TRUE)
      if(train.design.met == "BLK") train.blk.ind <- blocking.design(seed = seed + k, num.array = ncol(handling.effect.tr))
      if(train.design.met == "STR") train.str.ind <- stratification.design(seed = seed + k, num.array = ncol(handling.effect.tr), batch.id = train.batch.id)

      train.design.met2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", train.design.met)))

      param <- ifelse(icombat, ", icombat = TRUE",
                      ifelse(isva, ", isva = TRUE",
                             ifelse(iruv, ", iruv = TRUE,
                                        biological.effect.ctrl = biological.effect.tr.ctrl,
                                        handling.effect.ctrl = handling.effect.tr.ctrl", "")))

      eval(parse(text = paste0(train.design.met2, ".tr <- rehybridize(biological.effect = biological.effect.tr,
                                   handling.effect = handling.effect.tr, group.id = group.id.tr,
                                   array.to.sample.assign = train.", train.design.met2, ".ind", param, ")")))


      eval(parse(text = paste0("assign_store[['train']][[", k, "]] <- list(train.", train.design.met2, ".ind)")))


    }





    # Rehybridize on test set

    if(test.design.met == "NONE"){
      test.design.met2 <- tolower(test.design.met)
      eval(parse(text = paste0(test.design.met2, ".te <- biological.effect.te")))
    } else{

      if(test.design.met == "CC+") test.cc1.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.te), degree = "complete", rev.order = FALSE)
      if(test.design.met == "CC-") test.cc2.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.te), degree = "complete", rev.order = TRUE)
      if(test.design.met == "PC+") test.pc1.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.te), degree = "partial", rev.order = FALSE)
      if(test.design.met == "PC-") test.pc2.ind <- confounding.design(seed = seed + k, num.array = ncol(handling.effect.te), degree = "partial", rev.order = TRUE)
      if(test.design.met == "BLK") test.blk.ind <- blocking.design(seed = seed + k, num.array = ncol(handling.effect.te))
      if(test.design.met == "STR") test.str.ind <- stratification.design(seed = seed + k, num.array = ncol(handling.effect.te), batch.id = test.batch.id)

      test.design.met2 <- tolower(gsub("[-]", "2", gsub("[+]", "1", test.design.met)))

      param <- ifelse(icombat, ", icombat = TRUE",
                      ifelse(isva, ", isva = TRUE",
                             ifelse(iruv, ", iruv = TRUE,
                                        biological.effect.ctrl = biological.effect.te.ctrl,
                                        handling.effect.ctrl = handling.effect.te.ctrl", "")))

      eval(parse(text = paste0(test.design.met2, ".te <- rehybridize(biological.effect = biological.effect.te,
                                   handling.effect = handling.effect.te, group.id = group.id.te,
                                   array.to.sample.assign = test.", test.design.met2, ".ind", param, ")")))


      eval(parse(text = paste0("assign_store[['test']][", k, "] <- list(test.", test.design.met2, ".ind)")))

    }




    # # Unfinished
    # if(isva){ # isva has special output from rehybridize()
    #
    #   eval(parse(text = paste0("fsvaobj <- fsva(", dd2, ".tr$trainData, ",
    #                            dd2, ".tr$trainMod, ",
    #                            dd2, ".tr$trainSV, biological.effect.te)")))
    #
    #   eval(parse(text = paste0(dd2, ".tr <- fsvaobj$db # train.sva")))
    #   eval(parse(text = "biological.effect.te2 <- fsvaobj$new # test.fsva"))
    # } else {
    #   biological.effect.te2 <- biological.effect.te
    # }


    train.norm.list = c("MN", "QN", "VSN")
    new.norm.funcs <- switch.norm.funcs(norm.list = train.norm.list, norm.funcs = norm.funcs)


    # Normalization on training set
    if(train.norm.met == "NN"){

      # normalization function
      train.norm.met2 <- tolower(train.norm.met)

      # summarize
      eval(parse(text = paste0(train.design.met2, ".tr.", train.norm.met2, ".fin <- med.sum.pbset(", train.design.met2, ".tr)")))

    }
    else{

      # normalization function
      train.norm.met2 <- tolower(train.norm.met)
      train.norm.func <- new.norm.funcs[train.norm.list == train.norm.met]

      # normalize on training set
      eval(parse(text = paste0("temp <- ", train.norm.func,"(", train.design.met2, ".tr, ", test.design.met2, ".te)")))
      eval(parse(text = paste0(train.design.met2, ".tr.", train.norm.met2, " <- temp$train.", train.norm.met2)))

      #summarize
      eval(parse(text = paste0(train.design.met2, ".tr.", train.norm.met2, ".fin <- med.sum.pbset(", train.design.met2, ".tr.", train.norm.met2, ")")))

    }




    # MN and QN on test set
    if(test.norm.met %in% c("MN", "QN", "VSN")){

      # norm function
      test.norm.met2 <- tolower(test.norm.met)
      test.norm.func <- new.norm.funcs[train.norm.list == test.norm.met]

      # normalize on test set
      eval(parse(text = paste0("temp1 <- ", test.norm.func,"(", test.design.met2, ".te)")))   # Without dd2.tr
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, " <- temp1$train.", test.norm.met2)))

      # summarize
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, ".fin <- med.sum.pbset(gs.", test.design.met2, ".", test.norm.met2, ")")))

    }


    # fMN and fQN on test set
    if(test.norm.met %in% c("fMN", "fQN")){

      # normorlization function
      test.norm.met2 <- tolower(test.norm.met)
      test.norm.met.sub <- substr(test.norm.met, 2, 3)
      test.norm.met.sub2 <- tolower(test.norm.met.sub)
      test.norm.func <- new.norm.funcs[train.norm.list == test.norm.met.sub]


      # normalize on test set
      eval(parse(text = paste0("temp1 <- ", test.norm.func,"(", train.design.met2, ".tr, ", test.design.met2, ".te)")))   # With dd2.tr
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, " <- temp1$test.", test.norm.met2)))

      # summarize
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, ".fin <- med.sum.pbset(gs.", test.design.met2, ".", test.norm.met2, ")")))

    }

    # fVSN on test set
    if(test.norm.met %in% c("fVSN")){

      # normorlization function
      test.norm.met2 <- tolower(test.norm.met)
      test.norm.met.sub <- substr(test.norm.met, 2, 4)
      test.norm.met.sub2 <- tolower(test.norm.met.sub)
      test.norm.func <- new.norm.funcs[train.norm.list == test.norm.met.sub]


      # normalize on test set
      eval(parse(text = paste0("temp1 <- ", test.norm.func,"(", train.design.met2, ".tr, ", test.design.met2, ".te)")))   # With dd2.tr
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, " <- temp1$test.", test.norm.met2)))

      # summarize
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, ".fin <- med.sum.pbset(gs.", test.design.met2, ".", test.norm.met2, ")")))

    }


    # pMN and pQN on test set
    if(test.norm.met %in% c("pMN", "pQN")){

      # normorlization function
      test.norm.met2 <- tolower(test.norm.met)
      test.norm.met.sub <- substr(test.norm.met, 2, 3)
      test.norm.met.sub2 <- tolower(test.norm.met.sub)
      test.norm.func <- new.norm.funcs[train.norm.list == test.norm.met.sub]

      # combine training set and test set
      eval(parse(text = paste0("combine <- cbind(", train.design.met2, ".tr, ", test.design.met2, ".te)")))

      # normalize on combined set
      eval(parse(text = paste0("temp1 <- ", test.norm.func, "(combine)")))   # Without dd2.tr
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, " <- temp1$train.", test.norm.met.sub2, "[,129:192]"))) # Choose the test part

      # summarize
      eval(parse(text = paste0("gs.", test.design.met2, ".", test.norm.met2, ".fin <- med.sum.pbset(gs.", test.design.met2, ".", test.norm.met2, ")")))

    }



    # No normalization on test set
    if(test.norm.met == "NN"){

      # normalization function
      test.norm.met2 <- tolower(test.norm.met)

      # summarize
      eval(parse(text = paste0("gs.",test.design.met2, ".", test.norm.met2, ".fin <- med.sum.pbset(", test.design.met2, ".te)"))) # non normalization

    }


    class.list = c("PAM", "LASSO", "ClaNC", "ranFor", "SVM", "kNN", "DLDA", "custom")
    new.class.funcs <- switch.classifier.funcs.class(class.list = class.list,
                                                     class.funcs = class.funcs,
                                                     pred.funcs = pred.funcs)

    # Do classification

    cat(paste0("- build ", class.met, " (Train: ", train.design.met, " & ", train.norm.met, "; Test: ", test.design.met, " & ", test.norm.met, ")\n"))
    class.func <- new.class.funcs$build.funcs[class.list == class.met]
    pred.func <- new.class.funcs$pred.funcs[class.list == class.met]

    # build
    eval(parse(text = paste0(train.design.met2, ".", tolower(class.met), ".int.", train.norm.met2,
                             " <- ", class.func, "(kfold = 5, X = ", train.design.met2, ".tr.",
                             train.norm.met2, ".fin, y = group.id.tr, seed = ", seed + k, ")")))


    # store model
    eval(parse(text = paste0("model_store[", k, "] <- list(", train.design.met2, ".",
                             tolower(class.met), ".int.", train.norm.met2, ")")))
    # store feature
    # can extract later

    # predict
    eval(parse(text = paste0(test.design.met2, ".", tolower(class.met), ".pred.", test.norm.met2,
                             " <- ", pred.func, "(", train.design.met2, ".", tolower(class.met),
                             ".int.", train.norm.met2, ", gs.", test.design.met2, '.',
                             test.norm.met2, ".fin , group.id.te)")))


    # store coefficients
    # can extract later


    # store errors
    eval(parse(text = paste0("error_store[['internal']][", k, "]<- ", train.design.met2, ".",
                             tolower(class.met), ".int.", train.norm.met2, "$mc")))
    eval(parse(text = paste0("error_store[['external']][", k, "] <- ", test.design.met2, ".",
                             tolower(class.met), ".pred.", test.norm.met2, "$mc")))


    # store adjusted random index
    eval(parse(text = paste0("ari_store[", k, "]<- adjustedRandIndex(group.id.te, ",
                             test.design.met2, ".", tolower(class.met), ".pred.", test.norm.met2, "$pred)")))



  } # simulation loop


  return(list(assign_store = assign_store,
              model_store = model_store,
              error_store = error_store,
              ari_store = ari_store))

}


