# All the necessary functions for PLS methods comparison for Feature Selection 

# Note: use the document outline pane in RStudio to 
# navigate through this script


# GENERATE DATASETS #####

#++++++++++++++++++++++++++++++++++++++++++++
# draws random ids for a 2 class-problem
#++++++++++++++++++++++++++++++++++++++++++++

idsimul <- function(class, size, nsimul) {
  # class: vector of class, must be in [0,1]
  # size: size of the total dataset
  # nsimul: number of simulations
  
  class <- as.factor(class)
  class <- relevel(class , ref = "0")
  tab <- table(class)
  prop <- tab/sum(tab)
  
  i0 <- which(class == levels(class)[1])
  i1 <- which(class == levels(class)[2])
  
  ind0 <- matrix(data = NA, ncol = nsimul, nrow = round(prop[1]*size))
  ind1 <- matrix(data = NA, ncol = nsimul, nrow = round(prop[2]*size))
  
  for (i in 1:nsimul) {
    ind0[,i] <- sample(i0, size = round(prop[1]*size))
    ind1[,i] <- sample(i1, size = round(prop[2]*size))   
  }
  
  idsimul <- list(ind0=ind0,ind1=ind1)
  
  return(idsimul)
  
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Generates nsimul subdatasets with randomised observations
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# x : spectral matrix
# y : response vector
# size : number of obs

genData <- function(x, y, nsimul, size, seed = 1){ 
  sizetot <- size # size of the total dataset
  
  pander("size of the total dataset:")
  pander(sizetot)
  
  # size: taille d'ech par class
  #### balanced classes
  # ia <- which(y==levels(as.factor(y))[1])
  # ib <- which(y==levels(as.factor(y))[2])
  
  set.seed(seed)
  idsimul_res <- idsimul(class = y, size = size, nsimul = nsimul)
  
  idsimul_res <- do.call(rbind, idsimul_res)
  
  x_boot <- mapply(function(id) x[idsimul_res[,id],], c(1:nsimul), SIMPLIFY = FALSE)
  y_boot <- mapply(function(id) y[idsimul_res[,id]], c(1:nsimul), SIMPLIFY = FALSE)
  
  # randomize the data
  id_random <- mapply(function(id) sample(1:sizetot,sizetot, replace = FALSE), 
                      c(1:nsimul) ,SIMPLIFY = TRUE)
  
  for (i in 1:nsimul){
    x_boot[[i]] <- x_boot[[i]][id_random[,i],]
    y_boot[[i]] <- y_boot[[i]][id_random[,i]]
  }
  
  
  names(x_boot) <- names(y_boot) <- paste0("Sim", 1:nsimul)
  
  n <- dim(x_boot[[1]])[1]
  m <- dim(x_boot[[1]])[2]
  
  prop <- sapply(y_boot, function(y) table(y)/size)
  
  return(list(xx=x_boot,yy=y_boot,n=n,m=m, prop  = prop))
  
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Creates K folds For K-fold CV (K=10)  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Keeps the same proportion of classes in the folds.
# Creates K train and test sets
# x : spectral matrix
# y : response vector
# k : number of folds

createKFolds <- function(x, y, k = 10) {
  
  df <- data.frame(rowname = rownames(x), Class = y)
  n <- dim(x)[1]
  
  # indication of in which fold are the samples
  folds <- caret::createFolds(df$Class, k = 10, list = FALSE)
  
  df$fold <- folds
  # mean value of the response per fold
  Prop1 <- ddply(df, "fold", summarise, prop=mean(Class))
  
  dftrain <- dftest <- vector("list", k)
  for (j in 1:k) {
    id <- which(folds == j)
    
    dftrain[[j]]$y <- y[-id]
    dftrain[[j]]$x <- x[-id,]
    names(dftrain[[j]]$y) <- rownames(dftrain[[j]]$x)
    
    dftest[[j]]$y <- y[id]
    dftest[[j]]$x <- x[id,]
    
  }
  
  print(paste0("nLevels:", paste0(sapply(dftest, function(x)
    nlevels(as.factor(x$y))), collapse=" ")))
  
  print(paste0("meanTestLevels", 
               paste0(sapply(dftest, function(x) mean(x$y)),collapse=" ")))
  
  return(list(dftrain=dftrain, dftest=dftest, proportion=Prop1))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# creates bootstrap samples a number of times
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# returns train and test samples
# x : spectral matrix
# y : response vector
# prop: proportion of the total size of the resampled datasets
# times: number of bootstrapped samples

createBoot <- function(y,x, times, prop) {

  size <- round(prop*length(y))
  
  id1 <- which(y==levels(as.factor(y))[1])
  size1 <- length(id1)
  
  id2 <- which(y==levels(as.factor(y))[2])
  size2 <- length(id2)
  
  mat <- matrix(data=NA, ncol=times, nrow=size)
  dftrain <- dftest <- vector(mode = "list", times)
  for (i in 1:times) {
    id1sampl <- sample(id1,size=ceiling(size/2))
    id2sampl <- sample(id2,size=ceiling(size/2))
    ids <- sample(c(id1sampl,id2sampl),replace = FALSE,
                  size = length(c(id1sampl,id2sampl)))
    mat[,i] <- ids
    
    dftrain[[i]]$y <- y[mat[,i]]
    dftrain[[i]]$x <- x[mat[,i],]
    
    dftest[[i]]$y <- y[-mat[,i]]
    dftest[[i]]$x <- x[-mat[,i],]
    
  }
  
  
  print(paste0("nLevels:", paste0(sapply(dftest, function(x)
    nlevels(as.factor(x$y))), collapse=" ")))
  print(paste0("meanTestLevels", 
               paste0(sapply(dftest, function(x) mean(x$y)),collapse=" ")))
  
  Prop1 <- sapply(dftest, function(x) mean(x$y, na.rm = TRUE))
  names(Prop1) <- paste0("boot", 1:times)
  
  return(list(dftrain=dftrain,  dftest=dftest, proportion=Prop1))
  
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# creates 1 bootstrap sample
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# returns train and test samples
# x : spectral matrix
# y : response vector
# prop: proportion of the total size of the resampled datasets


create1Boot <- function(y,x,  prop) {
  # prop: proportion equal to the total size of the resampled datasets
  
  size <- round(prop*length(y))
  
  id1 <- which(y==levels(as.factor(y))[1])
  size1 <- length(id1)
  
  id2 <- which(y==levels(as.factor(y))[2])
  size2 <- length(id2)
  
  mat <- matrix(data=NA, ncol=1, nrow=size)
  dftrain <- dftest <- list()
  id1sampl <- sample(id1,size=ceiling(size/2))
  id2sampl <- sample(id2,size=ceiling(size/2))
  ids <- sample(c(id1sampl,id2sampl),replace = FALSE,
                size = length(c(id1sampl,id2sampl)))
  mat <- ids
  
  dftrain$y <- y[mat]
  dftrain$x <- x[mat,]
  
  dftest$y <- y[-mat]
  dftest$x <- x[-mat,]
  
  print(paste0("nLevels: ",  nlevels(as.factor(dftest$y)), collapse=" "))
  print(paste0("meanTestLevels: ", 
               paste0(mean(dftest$y),collapse=" ")))
  
  Prop1 <- mean(dftest$y, na.rm = TRUE)
  return(list(dftrain=dftrain,  dftest=dftest, proportion=Prop1))
  
}

# ROC #####

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# distance to the optimal situation in ROC curves
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# x and y : coordinates

ROC_dist_fun <- function(x, y) {
  which.min(sqrt((0-x)^2 + (1-y)^2))
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Computes the AUC
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# function from https://blog.revolutionanalytics.com/2016/11/calculating-auc.html
# TPR: true positive rate
# TPR: false positive rate

AUC_fun <- function(TPR, FPR){
  # inputs already sorted, best scores first 
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ROC indices for 1 threshold, and vectors of true and predicted responses 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# t: threshold
# y_true: true response
# y_pred: predicted response

ROC_index <- function(t, y_true, y_pred){
  # t: threshold
  # assuming to have y = 0 (control) or 1 (case)
  
  y_pred_bin <- rep(0, length(y_true))
  y_pred_bin[as.numeric(as.character(y_pred)) >= t] <- 1
  
  CM <- matrix(NA, ncol=2, nrow=2, 
               dimnames = list(c("pred1","pred0"), c("true1", "true0")))
  CM[1,1] <- sum(y_pred_bin==1 & y_true==1)
  CM[1,2] <- sum(y_pred_bin==1 & y_true==0)
  CM[2,1] <- sum(y_pred_bin==0 & y_true==1)
  CM[2,2] <- sum(y_pred_bin==0 & y_true==0)
  
  # CM <- table(y_pred_bin, y_true)
  TP <- CM[1,1]
  FP <- CM[1,2]
  FN <- CM[2,1]
  TN <- CM[2,2]
  
  sensit <- ifelse(is.na(TP/(TP+FN)),0, TP/(TP+FN))
  specif <- ifelse(is.na(TN/(TN+FP)),0, TN/(TN+FP))    
  PPV <- ifelse(is.na(TP/(TP+FP)),0, TP/(TP+FP))   
  NPV <- ifelse(is.na(TN/(TN+FN)),0, TN/(TN+FN))  
  
  AvPredAb <- (PPV + NPV)/2
  AvConfClas <- (sensit + specif)/2
  
  MCC <- (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)))
  MCC <- ifelse(is.na(MCC),0, 
                (TP*TN - FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  F1 <- (2*TP)/(2*TP+FP+FN)
  F1 <- ifelse(is.na(MCC),0, (2*TP)/(2*TP+FP+FN))
    
  list(y_pred_bin=y_pred_bin, sensit=sensit, specif=specif, PPV=PPV, NPV=NPV, 
       AvPredAb=AvPredAb,  AvConfClas=AvConfClas, MCC=MCC, F1=F1)
}

# ROC_index(t = 0.5, y_true = c(1,1,1,1,1,1), y_pred = c(1,1,1,0,1,1))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ROC indices with all possible thresholds
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X: 1 redicted response vector
# y_true: 1 true response vector

ROC_indices <- function(X, y_true){
  # Function with all the possible thresholds, 
  # one vector y_true and one vector of y_pred
  
  y_pred <- X
  
  thresh <- seq(from = min(y_pred), to = max(y_pred), 
                by = diff(range(y_pred))/100)
  
  res <- sapply(thresh, FUN = function(v) 
    ROC_index(t=v, y_pred=y_pred, y_true = y_true),
    simplify = TRUE, USE.NAMES = TRUE)
  
  sensit_res <- unlist(res["sensit",])
  specif_res <- unlist(res["specif",])
  AvPredAb_res <- unlist(res["AvPredAb",])
  AvConfClas_res <- unlist(res["AvConfClas",])
  PPV_res <- unlist(res["PPV",])
  F1_res <- unlist(res["F1",])
  
  id <- ROC_dist_fun(x= (1-specif_res), y = sensit_res)
  
  Opt_thres <- thresh[id]
  
  
  # at best threshold
  sensit <- sensit_res[id]
  specif <- specif_res[id]
  AvPredAb <- AvPredAb_res[id]
  AvConfClas <- AvConfClas_res[id]
  PPV <- PPV_res[id]
  F1 <- F1_res[id]
  
  return(list(sensit_res=sensit_res, specif_res=specif_res, 
              AvPredAb_res=AvPredAb_res, PPV_res=PPV_res,F1_res=F1_res,
              AvConfClas_res=AvConfClas_res, Opt_thres=Opt_thres,thresh=thresh,
              Opt_thres_res = list(sensit=sensit,specif=specif,
                                   AvPredAb =AvPredAb,AvConfClas= AvConfClas, PPV=PPV, 
                                   F1=F1)))
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Function with all the possible thresholds, one vector y_true and one list of y_pred folds
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X: predicted response vector
# y_true: list of true response vectors 

ROC_indices_list <- function(X, y_true){
  # Function with all the possible thresholds, 
  # one vector y_true and one list of y_pred folds

  y_pred <- X
  sensit_res <- specif_res <- res <- AvPredAb_res <- 
    AvConfClas_res <- PPV_res <- F1_res <- vector(mode = "list")
  Opt_thres <- id <- c()
  
  for (i in 1:length(y_pred)){
    
    thresh <- seq(from = min(y_pred[[i]]), to = max(y_pred[[i]]), 
                  by = diff(range(y_pred[[i]]))/100)
    
    res[[i]] <- sapply(thresh, FUN = function(v) 
      ROC_index(t=v, y_pred=y_pred[[i]], y_true = y_true[[i]]),
      simplify = TRUE, USE.NAMES = TRUE)
    
    sensit_res[[i]] <- unlist(res[[i]]["sensit",])
    specif_res[[i]] <- unlist(res[[i]]["specif",])
    AvPredAb_res[[i]] <- unlist(res[[i]]["AvPredAb",])
    AvConfClas_res[[i]] <- unlist(res[[i]]["AvConfClas",])
    PPV_res[[i]] <- unlist(res[[i]]["PPV",])
    F1_res[[i]] <- unlist(res[[i]]["F1",])
    
    id[i] <- ROC_dist_fun(x= (1-specif_res[[i]]), y = sensit_res[[i]] )
    
    Opt_thres[i] <- thresh[id[i]]
    
  }
  
  sensit <- specif <- AvPredAb <- AvConfClas <- PPV <- F1 <- c()
  for (i in 1:length(y_pred)){
    sensit[i] <- sensit_res[[i]][id[i]]
    specif[i] <- specif_res[[i]][id[i]]
    AvPredAb[i] <- AvPredAb_res[[i]][id[i]]
    AvConfClas[i] <- AvConfClas_res[[i]][id[i]]
    PPV[i] <- PPV_res[[i]][id[i]]
    F1[i] <- F1_res[[i]][id[i]]
  }
  
  return(list(sensit_res=sensit_res, specif_res=specif_res, 
              AvPredAb_res=AvPredAb_res,PPV_res=PPV_res, 
              AvConfClas_res=AvConfClas_res,F1_res=F1_res, Opt_thres=Opt_thres,
              Opt_thres_res = list(sensit=sensit,specif=specif,
                                   AvPredAb =AvPredAb,AvConfClas= AvConfClas,
                                   PPV=PPV,F1=F1)))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Computes the RMSE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X: predicted response
# y_true: true response

RMSE_fun <- function(X, y_true){
  # X:y_pred
  y_pred <- X
  RMSEP <- sqrt(sum((y_true - y_pred)^2)/length(y_pred))
  RMSEP
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Constr_ROC_feat_rank
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X = n biomarkers selected (threshold) 
# ranked_B = matrix with ranked feature indices; columns = simulations
# Bexp = true biomarkers

Constr_ROC_feat_rank <- function(X, ranked_B, Bexp) {
  nbiom <- X
  Ntot <- dim(ranked_B)[1]
  sel_feat <- matrix(ranked_B[1:nbiom,], nrow = nbiom, byrow = FALSE)
  notsel_feat <- ranked_B[min(nbiom+1,Ntot):Ntot,]
  
  sensit <- apply(sel_feat, 2, function(a) sum(is.element(a,Bexp))/length(Bexp))
  fdr <- apply(sel_feat, 2, function(a) sum(!is.element(a,Bexp))/nbiom)
  precision <-  apply(sel_feat, 2, function(a) sum(is.element(a,Bexp))/nbiom)
  # specif <- apply(notsel_feat, 2, function(a) sum(!is.element(a,Bexp))/(Ntot-x))
  ppv <- 1 - fdr
  F1 <- 2*((ppv*sensit)/(ppv+sensit))
  
  sdsensit <- sd(sensit)
  # sdspecif <- sd(specif)
  sdfdr <- sd(fdr)
  sdppv <- sd(ppv)
  sdF1 <- sd(F1) 
  sdprecision <- sd(precision)
  
  msensit <- mean(sensit) # moyenne de sensitivite des simulations pour un nb et une methode
  # mspecif <- mean(specif)
  mfdr <- mean(fdr) # moyenne de fdr des simulations pour un nb et une methode
  mppv <- mean(ppv)
  mF1 <- mean(F1)
  mprecision <- mean(precision)
  
  
  return(list(msensit = msensit,mfdr = mfdr, mppv = mppv, mF1 = mF1,
              sdsensit = sdsensit, sdfdr=sdfdr,
              mprecision=mprecision,sdprecision=sdprecision, 
              sdF1=sdF1,sdppv=sdppv))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Constr_ROC_feat_sel
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ranked_B = matrix with ranked feature indices; columns = simulations
# Bexp = true biomarkers

Constr_ROC_feat_sel <- function(ranked_B, Bexp) {

  
  nbiom <- sapply(ranked_B, length)
  # Ntot <- dim(ranked_B)[1]
  
  sel_feat <- ranked_B
  
  # notsel_feat <- ranked_B[min(nbiom+1,Ntot):Ntot,]
  
  sensit <- sapply(sel_feat, function(a) sum(is.element(a,Bexp))/length(Bexp))
  
  fdr <- mapply(function(sel_feat, nbiom) sum(!is.element(sel_feat,Bexp))/nbiom, 
                sel_feat=sel_feat, nbiom=nbiom)
  
  precision <-  mapply(function(sel_feat, nbiom) sum(is.element(sel_feat,Bexp))/nbiom, 
                       sel_feat=sel_feat, nbiom=nbiom)
  
  ppv <- 1 - fdr
  F1 <- 2*((ppv*sensit)/(ppv+sensit))
  
  # specif <- apply(notsel_feat, 2, function(a) sum(!is.element(a,Bexp))/(Ntot-x))
  
  sdsensit <- sd(sensit)
  # sdspecif <- sd(specif)
  sdfdr <- sd(fdr)
  sdppv <- sd(ppv)
  sdF1 <- sd(F1) 
  sdprecision <- sd(precision)
  
  msensit <- mean(sensit) # moyenne de sensitivite des simulations pour un nb et une methode
  # mspecif <- mean(specif)
  mfdr <- mean(fdr) # moyenne de fdr des simulations pour un nb et une methode
  mppv <- mean(ppv)
  mF1 <- mean(F1)
  mprecision <- mean(precision)

  
  return(list(msensit = msensit,mfdr = mfdr,mppv = mppv, mF1 = mF1,
              sdsensit = sdsensit, sdfdr=sdfdr,
              mprecision=mprecision,sdprecision=sdprecision, 
              sdF1=sdF1,sdppv=sdppv))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Prediction accuracy
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# res_outerloop: list with truetest and predtest to be compared together

pred_accuracy_FUN <- function(res_outerloop){
  # Average the prediction accuracies over boot samples
  res_ROC_indices <- mapply(ROC_indices, y_true = res_outerloop$truetest, 
                            X = res_outerloop$predtest, SIMPLIFY = FALSE)
  Opt_thres_res <- sapply(res_ROC_indices, function(a) a[["Opt_thres_res"]], 
                          simplify = TRUE)
  
  Opt_thres_res <- matrix(unlist(Opt_thres_res), dimnames = dimnames(Opt_thres_res), 
                          ncol = dim(Opt_thres_res)[2],nrow = dim(Opt_thres_res)[1], byrow = FALSE)
  
  mean_res_ROC_indices <- apply(Opt_thres_res,1,mean)
  
  specif_res <- sapply(res_ROC_indices, function(x) x[["specif_res"]])
  mspecif <- rowMeans(specif_res)
  sensit_res <- sapply(res_ROC_indices, function(x) x[["sensit_res"]])
  msensit <- rowMeans(sensit_res)
  
  
  AvConfClas_res <- sapply(res_ROC_indices, function(x) x[["AvConfClas_res"]])
  mAvConfClas <- rowMeans(AvConfClas_res)
  AvPredAb_res <- sapply(res_ROC_indices, function(x) x[["AvPredAb_res"]])
  mAvPredAb <- rowMeans(AvPredAb_res)
  F1_res <- sapply(res_ROC_indices, function(x) x[["F1_res"]])
  mF1 <- rowMeans(F1_res)
  
  return(list(mean_res_ROC_indices = mean_res_ROC_indices, mspecif = mspecif, 
              msensit = msensit, mAvConfClas = mAvConfClas, mAvPredAb = mAvPredAb, 
              mF1 = mF1))
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Extracts selected features and a binary matrix coding for the selected features 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# type = c("thresh")
# rank_B =rank_B_byeta[[1]]
# ranked_B= ranked_B_sel_allEta[[j]][[1]]
# thresh = median(abs(rank_B))
# nbiom = 10

# FS stability ##################

ExtrFeat <- function(rank_B, ranked_B, type = c("thresh", "nbiom"), 
                     thresh = median(abs(rank_B)), nbiom = 10){
  # Extracts selected features and a binary matrix coding for the selected features 
  type = match.arg(type)
  
  if (type=="thresh"){
    # threshold
    # thresh <- median(abs(rank_B_all$coefficients))
    FeatSel <- apply(rank_B, 1, function(a) which(abs(a)> thresh))
  } else{
    # nbiom
    FeatSel <- ranked_B[1:nbiom,]
    if(class(FeatSel)=="integer"){
      FeatSel  <- t(as.matrix(FeatSel))
    }
    FeatSel <- apply(FeatSel, 2, list)
    FeatSel <- lapply(FeatSel,unlist)
  }
  
  binMat_feat<- vector(mode="list")
  binMat_feat <- matrix(0, ncol = dim(rank_B)[2], nrow = dim(rank_B)[1],
                        dimnames=list(NULL, 1:dim(rank_B)[2]))
  for (i in 1:dim(rank_B)[1]){
    binMat_feat[i,FeatSel[[i]]] <- 1
  }
  
  
  
  return(list(FeatSel=FeatSel, binMat_feat=binMat_feat))
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ASMstab
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From Lustgarten, J. L., Gopalakrishnan, V., & Visweswaran, S. (2009). 
# Measuring stability of feature selection in biomedical datasets. 
# In AMIA annual symposium proceedings (Vol. 2009, p. 406). 
# American Medical Informatics Association.

# FeatSets: list with selected features
# n: n features

ASMstab <- function(FeatSets, n){
  c <- length(FeatSets)
  sim <- matrix(NA, ncol = (c-1), nrow = (c-1))
  for (i in 1:(c-1)){
    for (j in 2:c){ 
      
      Si <- FeatSets[[i]]
      Sj <- FeatSets[[j]]
      
      ki <- length(Si)
      kj <- length(Sj)
      
      r <- sum(FeatSets[[i]] %in% FeatSets[[j]])
      
      sim[i,(j-1)] <- (r-((ki*kj)/n))/(min(ki,kj) - max(0, ki+kj-n))
      
    }
  }
  
  return(sum(sim)/(c*(c-1)))
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# getStability
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From : Nogueira, S., Sechidis, K., & Brown, G. (2017). 
# On the Stability of Feature Selection Algorithms. Journal 
# of Machine Learning Research, 18, 174-1.

getStability <- function(X,alpha=0.05) {
  ## the input X is a binary matrix of size M*d where:
  ## M is the number of bootstrap replicates
  ## d is the total number of features
  ## alpha is the level of significance (e.g. if alpha=0.05, we will get 95% confidence intervals)
  ## it's an optional argument and is set to 5% by default
  ### first we compute the stability
  
  M<-nrow(X)
  d<-ncol(X)
  hatPF<-colMeans(X)
  kbar<-sum(hatPF)
  v_rand=(kbar/d)*(1-kbar/d)
  stability<-1-(M/(M-1))*mean(hatPF*(1-hatPF))/v_rand ## this is the stability estimate
  
  ## then we compute the variance of the estimate
  ki<-rowSums(X)
  phi_i<-rep(0,M)
  for(i in 1:M){ 
    phi_i[i]<-(1/v_rand)*((1/d)*sum(X[i,]*hatPF)-(ki[i]*kbar)/d^2-(stability/2)*((2*kbar*ki[i])/d^2-ki[i]/d-kbar/d+1))
  }
  phi_bar=mean(phi_i)
  var_stab=(4/M^2)*sum((phi_i-phi_bar)^2) ## this is the variance of the stability estimate
  
  ## then we calculate lower and upper limits of the confidence intervals
  z<-qnorm(1-alpha/2) # this is the standard normal cumulative inverse at a level 1-alpha/2
  upper<-stability+z*sqrt(var_stab) ## the upper bound of the (1-alpha) confidence interval
  lower<-stability-z*sqrt(var_stab) ## the lower bound of the (1-alpha) confidence interval
  
  return(list(stability=stability,variance=var_stab,lower=lower,upper=upper))
  
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VIP
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# W = loading.weights
# Q = Yloadings
# TT = scores
# opt.comp: number of components


VIP_fun <- function(W, Q, TT,opt.comp){
  ### from the plsVarSel package
  # Arguments:
  # W = loading.weights
  # Q = Yloadings
  # TT = scores
  W <- as.matrix(W)
  TT <- as.matrix(TT)
  p <- dim(W)[1]
  Q2 <- as.numeric(Q) * as.numeric(Q)
  Q2TT <- Q2[1:opt.comp] * diag(crossprod(TT))[1:opt.comp]
  WW <- W * W/apply(W, 2, function(x) sum(x * x))
  VIP <- sqrt(p * apply(sweep(WW[, 1:opt.comp, drop = FALSE],
                              2, Q2TT, "*"), 1, sum)/sum(Q2TT))
  VIP
}

# PLS ###########

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop for PLS 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# trainIn: list of folds with with x and y for the training set
# testIn: list of folds with with x and y for the test set
# nLVmax: number of maximum components for PLS
# mc.cores: number of cores for mclapply

innerloop_PLS <- function(trainIn, testIn, nLVmax = 15, mc.cores = 2){
  
  # PLS-DA model 
  
  pls_fun <- function(j){
      predictions <- vector(mode="list")
      
      df <- as.data.frame(cbind(y = trainIn[[j]]$y, trainIn[[j]]$x))
      res_mvr <- mvr(y ~ . -1, ncomp = nLVmax, data = df ,
                     validation = "none", method = "simpls") # method = SIMPLS
      
      predict(object = res_mvr, 
              newdata = testIn[[j]]$x)
      
  }
    
  predtest <- mclapply(1:length(trainIn), pls_fun, mc.cores = mc.cores)
    
  
  predtestList <- vector("list")
  for (i in 1:nLVmax){
    predtestList[[i]] <- sapply(predtest, function(x) x[,1,i], simplify = FALSE)
  }
  
  y_true <-  sapply(testIn, function(x) x$y, simplify=FALSE)
  
  
  return(list(predtestList = predtestList, y_true = y_true))
  
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop evaluation for PLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# predtestList : list returned from innerloop_PLS
# y_true : list of true response
# mc.core : number of cores for mclapply

innerloop_PLS_eval <- function(predtestList, y_true,mc.cores=2){ 
  
  nLVmax <- length(predtestList)
  
  res_ROC_indices <- mclapply(predtestList, ROC_indices_list, 
                              y_true = y_true, mc.cores=mc.cores)
  res_ROC_indices <- simplify2array(res_ROC_indices)
  
  AvPredAb <- mclapply(res_ROC_indices["Opt_thres_res",], 
                       function(x) x$AvPredAb,
                       mc.cores=mc.cores)
  AvPredAb <- simplify2array(AvPredAb)
  
  ConfClas <- mclapply(res_ROC_indices["Opt_thres_res",], function(x) x$AvConfClas,mc.cores=mc.cores)
  ConfClas <- simplify2array(ConfClas)
  
  mean_AvPredAb <- apply(AvPredAb,2,mean)
  plot(mean_AvPredAb, type="b", main = "mean AvPredAb")
  
  mean_ConfClas <- apply(ConfClas,2,mean)
  plot(mean_ConfClas, type="b", main = "mean ConfClas")
  
  sensit_res <- res_ROC_indices["sensit_res",]
  specif_res <- res_ROC_indices["specif_res",]
  
  res_AUC_fun <- matrix(NA, ncol=length(sensit_res[[1]]), nrow=nLVmax)
  for (i in 1:nLVmax){
    for (j in 1:length(sensit_res[[1]])){
      res_AUC_fun[i,j] <- abs(AUC_fun(TPR = sensit_res[[i]][[j]], 
                                      FPR = (1-specif_res[[i]][[j]])))
    }
  }
  
  mean_AUC <- apply(res_AUC_fun,1,mean)
  plot(mean_AUC, type="b", main = "mean AUC")
  
  
  #  RMSEP  
  res_RMSEP <- sapply(predtestList, function(x) mapply(RMSE_fun, X=x, y_true=y_true))
  mean_RMSEP <- apply(res_RMSEP,2,mean)
  plot(mean_RMSEP, type="b", main = "mean RMSEP")
  
  # final criterion for the optimal number of nLV +++++
  # nLV_opt <- min(15,which.min(mean_RMSEP))
  nLV_opt <- min(15,which.max(mean_AUC))
  
  return(nLV_opt)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Outer loop for PLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all_subsamples : list with elements (also lists): 
  # trainOut, 
  # testOut, 
  # trainIn,
  # testIn
  # trainOutX
  # trainOuty (cf. Comp_methods_main.rmd)
# nLVmax : number of LV max in PLS
# mc.cores : number of cores for mclapply

outerloop_PLS <- function(all_subsamples,  nLVmax = 15, mc.cores=2){
  
  trainOut <- all_subsamples$trainOut
  testOut <- all_subsamples$testOut
  trainIn <- all_subsamples$trainIn
  testIn <- all_subsamples$testIn
  # trainOutX <- all_subsamples$trainOutX
  # trainOuty <- all_subsamples$trainOuty
  
  # PLS-DA model +++++
  
  truetest <- sapply(testOut, function(x) x$y, simplify = FALSE) # true value of test sets
  
  predtest <- HP <- vector(mode="list", length = length(trainIn))
  rank_B_all <- vip1 <- vector(mode="list")
  coefficients <- vip <- matrix(NA, ncol = dim(trainOut[[1]]$x)[2], 
                                nrow = length(trainIn))
  
  for (j in 1:length(trainIn)){
    
    #### innerloop results
    HP[[j]] <- innerloop_PLS(trainIn = trainIn[[j]], testIn =  testIn[[j]], 
                             nLVmax = nLVmax, mc.cores=mc.cores)
    
    HP[[j]] <- innerloop_PLS_eval(predtestList = HP[[j]]$predtestList, 
                                  y_true = HP[[j]]$y_true,mc.cores=mc.cores)
    
    #### outerloop results
    
    # PLS +++++++
    df <- as.data.frame(cbind(y = trainOut[[j]]$y, trainOut[[j]]$x))
    res_mvr <- pls::mvr(y ~ . -1, ncomp = HP[[j]], data = df,
                        validation = "none", method = "simpls") # method = SIMPLS
    names(HP) <- rep("nLV",length(HP))
    
    coefficients[j,] <- res_mvr$coefficients[,,HP[[j]]]
    res_mvr$loading.weights <- res_mvr$projection # missing loading.weights when simpls...
    
    vip[j,] <- plsVarSel::VIP(res_mvr, opt.comp = HP[[j]])
    
    predtest[[j]] <- as.vector(predict(object = res_mvr, newdata = testOut[[j]]$x,
                                       ncomp = HP[[j]]))
    
    vip1[[j]] <- which(vip[j,]>=1)
  }  
  
  rank_B_all$PLS_coef <- coefficients
  rank_B_all$PLS_vip <- vip
  rank_B_all$PLS_vip1 <- vip
  
  # Sorted ranked indexes for all the simulations for VIP and coef
  ranked_B_all <- sapply(rank_B_all, 
                         function(x) apply(x,1, function(a) order(abs(a), 
                                                                  decreasing=TRUE)),
                         simplify=FALSE, USE.NAMES = TRUE)
  # Sorted ranked indexes for all the simulations for VIP1
  ranked_B_all$PLS_vip1 <- vip1
    
  return(list(predtest = predtest, truetest = truetest, HP = HP, 
              rank_B_all = rank_B_all, ranked_B_all = ranked_B_all, 
              trainOut=trainOut))
}

# RIDGE ###########

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Ridge regression and evaluation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# trainIn: list of folds with with x and y for the training set
# testIn: list of folds with with x and y for the test set

innerloop_RIDGE_eval <- function(trainIn, testIn){
  
  X = rbind(trainIn[[1]]$x, testIn[[1]]$x)
  Y = c(trainIn[[1]]$y, testIn[[1]]$y)
  
  allXtest <- sapply(testIn, function(x) x$x, simplify = FALSE)
  
  # sample ID by fold
  index <- sapply(allXtest, function(x) rownames(x), simplify = FALSE)
  
  all_ids <- unique(unlist(index))
  foldnum <- c()
  for (i in 1:10){
    foldnum[index[[i]]] <- i
  }
  
  # innerloop_RIDGE ===
  lambda2 <- optL2(response = Y, penalized = X, lambda1 = 0, minlambda2=1e-13,
                   model = "logistic", trace = FALSE, fold = foldnum)$lambda
  
  return(lambda2)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# outer loop for Ridge
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all_subsamples : list with elements (also lists): 
  # trainOut, 
  # testOut, 
  # trainIn,
  # testIn
  # trainOutX
  # trainOuty (cf. Comp_methods_main.rmd)

outerloop_RIDGE <- function(all_subsamples){
  
  trainOut <- all_subsamples$trainOut
  testOut <- all_subsamples$testOut
  trainIn <- all_subsamples$trainIn
  testIn <- all_subsamples$testIn
  
  
  truetest <- sapply(testOut, function(a) as.vector(a$y), simplify = FALSE) # true value of test sets
  
  predtest <- lambda2_opt <- vector(mode = "list", length = length(trainIn))
  rank_B_all <- vector(mode = "list")
  coefficients <-  matrix(NA, ncol = dim(trainOut[[1]]$x)[2], 
                          nrow = length(trainIn))
  
  for (j in 1:length(trainIn)){
    
    #### innerloop results
    res_eval <- innerloop_RIDGE_eval(trainIn=trainIn[[j]], 
                                     testIn=testIn[[j]])
    lambda2_opt[[j]] <- res_eval
    
    #### outerloop results
    
    # Regression logistique ridge avec lambda optimisé
    fit.ridge <- penalized(response = trainOut[[j]]$y, penalized = trainOut[[j]]$x,
                           model = "logistic",lambda1 = 0, 
                           lambda2 = lambda2_opt[[j]] , trace = FALSE)
    
    coefficients[j,] <- attr(fit.ridge, "penalized")
    
    predtest[[j]] <- predict(fit.ridge, penalized = testOut[[j]]$x)
    
  }  
  
  
  rank_B_all$RIDGE_coef <- coefficients
  
  # Sorted ranked indexes for all the simulations
  ranked_B_all <- sapply(rank_B_all, 
                         function(x) apply(x,1, function(a) order(abs(a), 
                                                                  decreasing=TRUE)),
                         simplify=FALSE, USE.NAMES = TRUE)
  
  return(list(predtest = predtest, truetest = truetest, lambda2_opt = lambda2_opt, 
              rank_B_all = rank_B_all, ranked_B_all = ranked_B_all, trainOut = trainOut))
}

# T-test ########
 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# outerloop for t-test
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all_subsamples : list with elements (also lists): 
  # trainOut, 
  # testOut, 
  # trainIn,
  # testIn
  # trainOutX
  # trainOuty (cf. Comp_methods_main.rmd)

outerloop_ttest <- function(all_subsamples){
  
  yy <- all_subsamples$trainOuty
  # yy <- mapply(function(a,b) c(a$y, b$y), a = all_subsamples$trainOut, 
  #              b = all_subsamples$testOut, SIMPLIFY = FALSE)
  
  xx <- all_subsamples$trainOutX
  # xx <- mapply(function(a,b) rbind(a$x, b$x), a = all_subsamples$trainOut, 
  #              b = all_subsamples$testOut, SIMPLIFY = FALSE)
  
  # loop for all bootstraped samples
  rank_B_all <- ranked_B_all <- vector(mode = "list")
  selection <- vector(mode = "list", length = 2)
  names(selection) <- c("ranked_B_all","rank_B_all")
  selection$rank_B_all <- selection$ranked_B_all <- vector(mode = "list")
  
  # ranked_B_all$tstat <- vector(mode = "list", length = length(xx))
  ranked_B_all$tstat <- matrix(NA, nrow = dim(xx[[1]])[2], ncol = length(xx))
  rank_B_all$tstat <- selection$rank_B_all$tstat <-  matrix(NA, ncol = dim(xx[[1]])[2], nrow = length(xx))
  selection$ranked_B_all$tstat <- vector(mode = "list")
  
  for(j in 1:length(xx)){
    tstat <- pval <- ratio <- c()
    for (i in 1:dim(xx[[1]])[2]){
      xx0 <- xx[[j]][yy[[j]]==0,i]
      xx1 <- xx[[j]][yy[[j]]==1,i]
      
      res <- t.test(x=xx1, y = xx0,
                    alternative = c("two.sided"),
                    var.equal = FALSE,
                    conf.level = 0.95)
      tstat[i] <- res$statistic
      pval[i] <- res$p.value
      
      ratio[i] <- mean(xx1)/mean(xx0)
    }
    
    # plot(ratio, type="l", main ="ratio")
    # plot(tstat, type="l", main ="tstat")
    
    # fdr corrected p-value
    pval[which(is.na(pval))] <- 1
    
    cor_pvalue <- p.adjust(p = pval, method = "fdr")
    
    tstat[which(is.na(tstat))] <- 0
    
    selection$rank_B_all$tstat[j,] <- tstat
    # fdr of 10%
    selection$ranked_B_all$tstat[[j]] <- which(cor_pvalue<=0.05)
    
    rank_B_all$tstat[j,] <- tstat
    ranked_B_all$tstat[,j] <- order(abs(rank_B_all$tstat[j,]), decreasing = TRUE)
    
    # volcano plot  
    # plot(log2(ratio),-log10(cor_pvalue), main = "volcano plot \n fdr pvalue correction", pch = 16)
    # abline(h =-log10(0.05), col="green") 
    # abline(v=c(-1,1), lty= 2) #  modulations de + ou - 2x
    # text(x = log2(ratio)[index], y = -log10(cor_pvalue)[index],
    #      labels = index, pos = c(2,3,4), cex = 0.7)
    
  }
  
  res_outerloop_ttest <- vector(mode = "list")
  res_outerloop_ttest$rank_B_all <- rank_B_all
  res_outerloop_ttest$ranked_B_all <- ranked_B_all
  res_outerloop_ttest$selection <- selection
  
  return(res_outerloop_ttest)
}

# SPLS ######

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop for SPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# trainIn: list of folds with with x and y for the training set
# testIn: list of folds with with x and y for the test set
# eta: sparsity parameter for SPLS
# K: number of LV for SPLS
# mc.cores: number of cores for mclapply

innerloop_SPLS <- function(trainIn, testIn, eta, K, mc.cores=2){
  
  # grid of Meta-Parameters
  grid_HP <- expand.grid(list(K=K, eta=eta))
  
  # SPLS FUN ++++
  spls_fun <- function(train, test, eta, K) {
    res_spls <- spls(x = train$x, y = train$y, 
                     K = K, eta = eta, kappa=0.5, 
                     select="pls2", fit="simpls",
                     scale.x=FALSE, scale.y=FALSE, eps=1e-4, 
                     maxstep=100, trace=FALSE)
    predictions <- predict(object = res_spls, newx = test$x)
    predictions
  }  
  
  # predictions from the inner loop ++++++
  inloop <- function(j){
    mapply(spls_fun, eta = grid_HP$eta, 
           K = grid_HP$K,MoreArgs = list(train =
                                           trainIn[[j]], test = testIn[[j]]))
  }

  predtest <- mclapply(X=c(1:length(trainIn)), FUN=inloop,
                       mc.cores=mc.cores)
 
  predtestList <- mclapply(X=c(1:dim(grid_HP)[1]), FUN=function(i)
    lapply(X = predtest, 
           FUN = function(x)
           {as.matrix(x[,i])}),
    mc.cores=mc.cores)
  
  
  y_true <-  sapply(testIn, function(x) x$y, simplify=FALSE)
  
  
  return(list(predtestList = predtestList, y_true = y_true, 
              grid_HP = grid_HP))
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop evaluation for SPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# res_innerloop_SPLS : list returned from innerloop_SPLS
# mc.core : number of cores for mclapply


innerloop_SPLS_eval <- function(res_innerloop_SPLS, mc.cores=2){ 
  
  # needed objects from res_innerloop_SPLS
  predtestList <- res_innerloop_SPLS$predtestList
  y_true <- res_innerloop_SPLS$y_true
  grid_HP <- res_innerloop_SPLS$grid_HP
  
  
  HP_num <- length(predtestList)
  
  res_ROC_indices <- mclapply(predtestList, ROC_indices_list, 
                              y_true = y_true,mc.cores=mc.cores)
  res_ROC_indices <- simplify2array(res_ROC_indices)
  
  AvPredAb <- mclapply(res_ROC_indices["Opt_thres_res",], 
                       function(x) x$AvPredAb,
                       mc.cores=mc.cores)
  AvPredAb <- simplify2array(AvPredAb)
  
  ConfClas <- mclapply(res_ROC_indices["Opt_thres_res",], function(x) x$AvConfClas,mc.cores=mc.cores)
  ConfClas <- simplify2array(ConfClas)
  
  mean_AvPredAb <- apply(AvPredAb,2,mean)
  plot(mean_AvPredAb, type="b", main = "mean AvPredAb")
  
  mean_ConfClas <- apply(ConfClas,2,mean)
  plot(mean_ConfClas, type="b", main = "mean ConfClas")
  
  sensit_res <- res_ROC_indices["sensit_res",]
  specif_res <- res_ROC_indices["specif_res",]
  
  
  res_AUC_fun <- matrix(NA, ncol=length(sensit_res[[1]]), nrow=HP_num)
  for (i in 1:HP_num){
    for (j in 1:length(sensit_res[[1]])){
      res_AUC_fun[i,j] <- abs(AUC_fun(TPR = sensit_res[[i]][[j]], 
                                      FPR = (1-specif_res[[i]][[j]])))
    }
  }
  
  mean_AUC <- apply(res_AUC_fun,1,mean)
  plot(mean_AUC, type="b", main = "mean AUC")
  
  
  #  RMSEP  +++
  res_RMSEP <- mclapply(predtestList, function(x) 
    mapply(RMSE_fun, X=x, y_true=y_true),
    mc.cores=mc.cores)
  res_RMSEP <- simplify2array(res_RMSEP)
  
  mean_RMSEP <- apply(res_RMSEP,2,mean)
  plot(mean_RMSEP, type="b", main = "mean RMSEP")
  
  # final criterion for the optimal number of K and eta ++++
  df <- data.frame(grid_HP, mean_AUC)
  Kbyeta_opt <- by(data = df$mean_AUC, INDICES = df$eta,  FUN = which.max)
  Kbyeta_opt <- as.matrix(Kbyeta_opt)
  colnames(Kbyeta_opt) <- "ncomp"
  
  # HP_opt_index <- which.min(mean_RMSEP)
  HP_opt_index <- which.max(mean_AUC)
   
  HP_opt <- grid_HP[HP_opt_index,]
  grid_HP <- cbind(grid_HP, mean_AUC=mean_AUC)
  
  return(list(HP_opt = HP_opt, Kbyeta_opt = Kbyeta_opt, grid_HP=grid_HP))
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Outerloop for SPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all_subsamples : list with elements (also lists): 
  # trainOut, 
  # testOut, 
  # trainIn,
  # testIn
  # trainOutX
  # trainOuty (cf. Comp_methods_main.rmd)
# eta: range of sparsity parameter for SPLS
# K: range of number of LV for SPLS
# mc.cores : number of cores for mclapply

outerloop_SPLS_allEta <- function(all_subsamples, 
                                  eta =  c(0.1,0.9), 
                                  K =  c(1:15),mc.cores=2){
  
  trainOut <- all_subsamples$trainOut
  testOut <- all_subsamples$testOut
  trainIn <- all_subsamples$trainIn
  testIn <- all_subsamples$testIn
  
  
  truetest <- sapply(testOut, function(a) as.vector(a$y), simplify = FALSE) # true value of test sets
  
  predtest <- Kbyeta_opt <- HP_opt <- index_opt <- coefficients <- 
    vector(mode = "list", length = length(trainIn))
  
  rank_B_all <- ranked_B_all <- vector(mode = "list")
  
  #### innerloop results
  inloop <- function(j){
    # apply inner loop
    res_innerloop_SPLS <- innerloop_SPLS(trainIn = trainIn[[j]], 
                                         testIn =  testIn[[j]], 
                                         eta = eta, K = K,
                                         mc.cores=mc.cores)
    
    # evaluate performances of the inner loop
    res_eval <- innerloop_SPLS_eval(res_innerloop_SPLS,
                                    mc.cores=mc.cores)

    index_opt <- which(rownames(res_eval$Kbyeta_opt) ==
                         as.character(res_eval$HP_opt$eta))

    return(list(HP_opt = res_eval$HP_opt,Kbyeta_opt=res_eval$Kbyeta_opt,
                index_opt = index_opt, grid_HP=res_eval$grid_HP))
  }
  
  res_inloop <- sapply(X=1:length(trainIn), FUN = inloop)
  HP_opt <- res_inloop["HP_opt",]
  Kbyeta_opt <- res_inloop["Kbyeta_opt",]
  index_opt <- res_inloop["index_opt",]
  grid_HP <- res_inloop["grid_HP",]$grid_HP
  
  #### outerloop results
  
  outloop <- function(i,j){
    eta_val <- as.numeric(rownames(Kbyeta_opt[[j]])[i]) 
    K_val <-  Kbyeta_opt[[j]][i,]
    res_spls <- spls(x = trainOut[[j]]$x, y = trainOut[[j]]$y, 
                     K = K_val, eta = eta_val, 
                     kappa=0.5,select="pls2", fit="simpls",
                     scale.x=FALSE, scale.y=FALSE, eps=1e-4, 
                     maxstep=100, trace=FALSE)
    res_spls
  }
  
  
  res_outloop <- sapply(1:length(trainIn), function(j) 
    sapply(1:dim(Kbyeta_opt[[j]])[1], function(i) 
      outloop(j = j, i = i), simplify = FALSE),
    simplify = FALSE)
  
  
  coefficients <- sapply(res_outloop, 
                         function(x) sapply(1:dim(Kbyeta_opt[[1]])[1], FUN = function(y) 
                           predict(x[[y]] , type = "coefficient")), 
                         simplify = FALSE)
  coefficients <- sapply(coefficients, FUN = function(x) t(x), simplify=FALSE)
  
  
  
  # best meta parameters +++++++
  best <- vector(mode = "list", length = 4)
  names(best) <- c("rank_B_all","ranked_B_all", "predtest", 
                   "truetest")
  best$rank_B_all <- best$ranked_B_all <- 
                  vector(mode = "list", length = 1)
  best$predtest <- best$truetest <- vector(mode = "list")
  names(best$rank_B_all) <- names(best$ranked_B_all) <-  "SPLS_coef"
  
  best$rank_B_all$SPLS_coef <- t(mapply(function(i, coef) coef[i,], 
                                        coef = coefficients, i = index_opt))
  
  best$ranked_B_all$SPLS_coef <- sapply(1:dim(best$rank_B_all$SPLS_coef)[1],
                                        function(i) 
                                          which(abs(best$rank_B_all$SPLS_coef[i,])> 0), 
                                        simplify = FALSE)
  
  
  testOutx <- sapply(testOut, function(x) x$x, simplify = FALSE)
  predtest <- mapply(FUN=function(x,y)  predict(object = x, newx = y), 
         x = res_outloop, y = testOutx, SIMPLIFY = FALSE)
  
  # use best MPs for final regression
  best$predtest <- mapply(function(i, pred) as.vector(pred[[i]]), 
                pred = predtest, i = index_opt, SIMPLIFY = FALSE)
  
  best$truetest <- truetest
  
  # all etas +++++++
  
  rank_B_all$SPLS_coef_eta <- coefficients
  
  fun_ranked <- function(i,j){
    which(abs(rank_B_all$SPLS_coef_eta[[j]][i,]) > 0)
  }
  
  ranked_B_all$SPLS_coef_eta <- sapply(1:dim(rank_B_all$SPLS_coef_eta[[1]])[1], 
                                       function(x) sapply(1:length(trainIn), 
                                                          function(y) fun_ranked(j=y, i=x),
                                                          simplify = FALSE), simplify = FALSE)
  names(ranked_B_all$SPLS_coef_eta) <- paste0("eta",rownames(Kbyeta_opt[[1]]))
  
  return(list(predtest = predtest, truetest = truetest, 
              Kbyeta_opt = Kbyeta_opt, HP_opt = HP_opt,
              rank_B_all = rank_B_all, ranked_B_all = ranked_B_all,
              trainOut = trainOut, best = best, grid_HP=grid_HP))
}


# OPLS ######
 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop for OPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# trainIn: list of folds with with x and y for the training set
# testIn: list of folds with with x and y for the test set
# nOcmax: Number of maximum components for OPLS
# mc.cores: number of cores for mclapply

innerloop_OPLS <- function(trainIn, testIn, nOcmax = 15, mc.cores=2){
  
  # OPLS-DA model +++++
  
  opls_fun <- function(train, test, nOc) {
    ropls <- OPLSDA(x=train$x, y=train$y, 
                    impT = FALSE, impG = FALSE, no = nOc)
    
    predtest <- OPLSDA_pred(ropls = ropls, x.new = test$x)
    predtest
  }
  
  # predictions ++++
    pred_fun <- function(j){
      mapply(opls_fun, nOc = c(1:nOcmax),
             MoreArgs = list(train = trainIn[[j]], test = testIn[[j]]),
             SIMPLIFY = FALSE)
    }
    
    predtest <- mclapply(1:length(trainIn), pred_fun, mc.cores = mc.cores)
 
  # plot(predtest[[1]][[15]])
  # abline(h=0.5)
  
  predtestList <- vector("list")
  for (i in 1:nOcmax){
    predtestList[[i]] <- sapply(predtest, function(x) x[[i]], simplify = FALSE)
  }
  predtestList[[1]]
  y_true <-  sapply(testIn, function(x) x$y, simplify=FALSE)
  
  
  return(list(predtestList = predtestList, y_true = y_true))
  
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop evaluation for OPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# predtestList : list returned from innerloop_OPLS
# y_true : list of true response
# mc.core : number of cores for mclapply

innerloop_OPLS_eval <- function(predtestList, y_true, mc.cores=2){ 
  
  nOcmax <- length(predtestList)
 
  res_ROC_indices <- mclapply(predtestList, ROC_indices_list, 
                              y_true = y_true,mc.cores=mc.cores)
  res_ROC_indices <- simplify2array(res_ROC_indices)
  
  AvPredAb <- mclapply(res_ROC_indices["Opt_thres_res",], 
                       function(x) x$AvPredAb,
                       mc.cores=mc.cores)
  AvPredAb <- simplify2array(AvPredAb)
  
  ConfClas <- mclapply(res_ROC_indices["Opt_thres_res",], function(x) x$AvConfClas,mc.cores=mc.cores)
  ConfClas <- simplify2array(ConfClas)
  
  mean_AvPredAb <- apply(AvPredAb,2,mean)
  plot(mean_AvPredAb, type="b", main = "mean AvPredAb")
  
  mean_ConfClas <- apply(ConfClas,2,mean)
  plot(mean_ConfClas, type="b", main = "mean ConfClas")
  
  sensit_res <- res_ROC_indices["sensit_res",]
  specif_res <- res_ROC_indices["specif_res",]
  
  
  res_AUC_fun <- matrix(NA, ncol=length(sensit_res[[1]]), nrow=nOcmax)
  for (i in 1:nOcmax){
    for (j in 1:length(sensit_res[[1]])){
      res_AUC_fun[i,j] <- abs(AUC_fun(TPR = sensit_res[[i]][[j]], 
                                      FPR = (1-specif_res[[i]][[j]])))
    }
  }
  
  mean_AUC <- apply(res_AUC_fun,1,mean)
  plot(mean_AUC, type="b", main = "mean AUC")
  
  
  #  RMSEP  ++++
  res_RMSEP <- sapply(predtestList, function(x) mapply(RMSE_fun, X=x, y_true=y_true))
  mean_RMSEP <- apply(res_RMSEP,2,mean)
  plot(mean_RMSEP, type="b", main = "mean RMSEP")
  
  # final criterion for the optimal number of no ++++
  # nOc_opt <- min(15,which.min(mean_RMSEP))
  nOc_opt <- min(15,which.max(mean_AUC))
  
  return(nOc_opt)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Outerloop for OPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all_subsamples : list with elements (also lists): 
  # trainOut, 
  # testOut, 
  # trainIn,
  # testIn
  # trainOutX
  # trainOuty (cf. Comp_methods_main.rmd)
# nOcmax : Number of maximum components for OPLS
# mc.cores : number of cores for mclapply

outerloop_OPLS <- function(all_subsamples,  nOcmax = 15, mc.cores=2){
  
  trainOut <- all_subsamples$trainOut
  testOut <- all_subsamples$testOut
  trainIn <- all_subsamples$trainIn
  testIn <- all_subsamples$testIn
  
  # OPLS-DA model -++++
  
  
  truetest <- sapply(testOut, function(x) x$y, simplify = FALSE) # true value of test sets
  
  rank_B_all <- vector(mode="list")
 
  #### innerloop results
  inloop <- function(j){
    nOcopt<- innerloop_OPLS(trainIn=trainIn[[j]], testIn=testIn[[j]], nOcmax = nOcmax)
    
    nOcopt <- innerloop_OPLS_eval(predtestList = nOcopt$predtestList, 
                                       y_true = nOcopt$y_true)
    nOcopt
  }
  
  nOcopt <- sapply(1:length(trainIn), inloop, simplify = FALSE)
  
  #### outerloop results
  outloop <- function(j){
    # OPLS ++++
    OPLSDA(x=trainOut[[j]]$x, y=trainOut[[j]]$y, impT = FALSE, impG = FALSE, 
                    no = nOcopt[[j]])
  }
  
  ropls <- sapply(1:length(trainIn), outloop, simplify = FALSE)
  
  coefficients <- t(sapply(ropls, function(x) x$b))
  
  vip <- t(sapply(ropls, function(x) VIP_fun(W = x$W, Q = x$C , TT = x$Tp, 
                                             opt.comp = 1)))
  
  predtest <- mapply(function(res_opls, testout) {as.vector(OPLSDA_pred(ropls = res_opls, 
                                                            x.new = testout$x))},
         res_opls=ropls, testout=testOut, SIMPLIFY = FALSE)

  vip1 <- vector(mode = "list")
  for (i in 1:nrow(vip)){
    vip1[[i]] <-which(vip[i,]>=1)
  }


  rank_B_all$OPLS_coef <- coefficients
  rank_B_all$OPLS_VIP <- vip
  rank_B_all$OPLS_VIP1 <- vip
  
  # Sorted ranked indexes for all the simulations for coef and VIP
  ranked_B_all <- sapply(rank_B_all, 
                         function(x) apply(x,1, 
                                   function(a) order(abs(a), decreasing=TRUE)),
                         simplify=FALSE, USE.NAMES = TRUE)
  
  # Sorted ranked indexes for all the simulations for VIP1
  ranked_B_all$OPLS_VIP1 <- vip1
  
  return(list(predtest = predtest, truetest = truetest, nOcopt = nOcopt, 
              rank_B_all = rank_B_all, ranked_B_all = ranked_B_all, 
              trainOut=trainOut))
}



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Manual OPLS:
opls_fun <- function(train, test, no){
  OPLSDA_res <- MBXUCL::OPLSDA(train$x, train$y, impT = FALSE, 
                               impG = FALSE, no = no)
  
  Xopls1 <- OPLSDA_res$Xopls
  Wortho <- OPLSDA_res$Wortho
  Portho <- OPLSDA_res$Portho
  
  y_pred <- OPLSDA_pred(ropls = OPLSDA_res, x.new = test$x)
  
  res_ROC_indices <- ROC_indices(X = y_pred, y_true = train$y)
  
  AUC <- abs(AUC_fun(TPR = res_ROC_indices$sensit_res, 
                     FPR = (1-res_ROC_indices$specif_res)))
  
  rmse <- RMSE_fun(X = y_pred, y_true = test$y)
  return(list(rmse = rmse, Xopls1 = Xopls1,Wortho=Wortho,Portho=Portho, AUC=AUC))
}

spls_fun <- function(train, test, no, eta, K, Xopls1, Wortho,Portho) {
  # SPLS on teh deflated matrix
  res_spls <- spls::spls(x=Xopls1, y=train$y, K = K, eta = eta,
                         kappa=0.5,select="pls2", 
                         fit="simpls", scale.x=FALSE, scale.y=FALSE, 
                         eps=1e-4,maxstep=100, trace=FALSE)
  ##+++++++++
  # In order to apply b directly on a new X, b becomes bcorr
  b <- as.vector(predict(res_spls, type = c("coefficient")))
  bcorr <- (diag(1, dim(Xopls1)[2], dim(Xopls1)[2]) - Wortho %*% solve(t(Portho) %*% Wortho) %*% t(Portho)) %*% b 
  
  x.new <- test$x - colMeans(train$x)
  predictions <- x.new %*% bcorr
  plot(predictions, type="p", col = (test$y+1))
  ##+++++++++
  # predictions <- predict(object = res_spls, newx = test$x)
  return(predictions)
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Inner loop for LSOPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# trainIn: list of folds with with x and y for the training set
# testIn: list of folds with with x and y for the test set
# nOmax: Number of maximum components for LSOPLS
# eta: sparsity parameter for SPLS
# K: number of LV for SPLS
# mc.cores: number of cores for mclapply

innerloop_LsOPLS <- function(trainIn, testIn, nOmax, eta, 
                             K, mc.cores=2){
  
  y_true <- sapply(testIn, function(x) x$y, simplify=FALSE)
  
  
  
  # inloop_opls +++++++
  inloop_opls <- function(j){
    sapply(1:nOmax, 
           function(i) opls_fun(train = trainIn[[j]], 
                                test = testIn[[j]], no = i))
  }
  
  res_opls <- mclapply(X=c(1:length(trainIn)), FUN=inloop_opls,
                       mc.cores=mc.cores)
  
  ### RMSE
  rmse <- matrix(unlist(sapply(res_opls,function(x) x["rmse",])), nrow =nOmax ,
                 byrow = FALSE)
  mean_rmse <- apply(rmse, 1, mean)
  # id <- which.min(mean_rmse)
  ### AUC
  AUC <- matrix(unlist(sapply(res_opls,function(x) x["AUC",])), nrow =nOmax ,
                byrow = FALSE)
  mean_AUC <- apply(AUC, 1, mean)
  id <- which.max(mean_AUC)
  
  
  
  Xopls1 <- sapply(res_opls, function(x,y) x["Xopls1",id])
  Portho <- lapply(res_opls, function(x,y) x["Portho",][[id]])
  Wortho <- lapply(res_opls, function(x,y) x["Wortho",][[id]])
  
  grid_HP <- expand.grid(list(K=K, eta=eta, nOc=id))
  
  # predictions +++++++
  
  inloop_spls <- function(j){
    mapply(spls_fun, eta = grid_HP$eta, 
           K = grid_HP$K, no = grid_HP$nOc, 
           MoreArgs = list(train = trainIn[[j]], test = testIn[[j]], 
                           Xopls1=Xopls1[[j]], Wortho=Wortho[[j]],Portho=Portho[[j]]))
  }
  
  predtest <- lapply(X=c(1:length(trainIn)), FUN=inloop_spls)
  
  
  predtestList <- mclapply(X=c(1:dim(grid_HP)[1]), FUN=function(i)
    lapply(X = predtest,  FUN = function(x)
           {as.matrix(x[,i])}),
    mc.cores=mc.cores)
  
  return(list(predtestList = predtestList, y_true = y_true, 
              grid_HP = grid_HP))
  
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# inner loop evaluation for LSOPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# innerloop_LsOPLS : list returned from innerloop_LsOPLS
# mc.core : number of cores for mclapply

innerloop_LsOPLS_eval <- function(innerloop_LsOPLS, mc.cores=2){ 
  
  # needed objects from res_innerloop_SPLS
  predtestList <- innerloop_LsOPLS$predtestList
  y_true <- innerloop_LsOPLS$y_true
  grid_HP <- innerloop_LsOPLS$grid_HP
  
  
  HP_num <- length(predtestList)
  
  res_ROC_indices <- mclapply(predtestList, ROC_indices_list, 
                              y_true = y_true,mc.cores=mc.cores)
  res_ROC_indices <- simplify2array(res_ROC_indices)
  
  sensit_res <- res_ROC_indices["sensit_res",]
  specif_res <- res_ROC_indices["specif_res",]
  
  res_AUC_fun <- matrix(NA, ncol=length(sensit_res[[1]]), nrow=HP_num)
  for (i in 1:HP_num){
    for (j in 1:length(sensit_res[[1]])){
      res_AUC_fun[i,j] <- abs(AUC_fun(TPR = sensit_res[[i]][[j]], 
                                      FPR = (1-specif_res[[i]][[j]])))
    }
  }
  
  mean_AUC <- apply(res_AUC_fun,1,mean)
  plot(mean_AUC, type="b", main = "mean AUC")
  
  
  #  RMSEP  +++++++
  res_RMSEP <- mclapply(predtestList, function(x) 
    mapply(RMSE_fun, X=x, y_true=y_true),
    mc.cores=mc.cores)
  res_RMSEP <- simplify2array(res_RMSEP)
  
  mean_RMSEP <- apply(res_RMSEP,2,mean)
  plot(mean_RMSEP, type="b", main = "mean RMSEP")
  
  # final criterion for the optimal number of K and eta +++++++
  df <- data.frame(grid_HP, mean_AUC)
  byeta_opt <- by(data = df, INDICES = list(df$eta),  
                  FUN = function(x) x[which.max(x$mean_AUC),])
  byeta_opt <- do.call(rbind, byeta_opt)
  byeta_opt <- as.matrix(byeta_opt)
  
  # HP_opt_index <- which.min(mean_RMSEP)
  HP_opt_index <- which.max(mean_AUC)
  
  HP_opt <- grid_HP[HP_opt_index,]
  
  return(list(HP_opt = HP_opt, byeta_opt = byeta_opt))
  
  
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Outerloop for LSOPLS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# all_subsamples : list with elements (also lists): 
  # trainOut, 
  # testOut, 
  # trainIn,
  # testIn
  # trainOutX
  # trainOuty (cf. Comp_methods_main.rmd)
# nOmax : Number of maximum components for OPLS
# eta: range of sparsity parameter for PLS
# K: range of number of LV for SPLS
# mc.cores : number of cores for mclapply

outerloop_LsOPLS_allEta <- function(all_subsamples, 
                                    eta =  c(0.1,0.9),nOmax = 3, 
                                    K =  c(1:15), mc.cores=2){
  
  trainOut <- all_subsamples$trainOut
  testOut <- all_subsamples$testOut
  trainIn <- all_subsamples$trainIn
  testIn <- all_subsamples$testIn
  
  truetest <- sapply(testOut, function(a) as.vector(a$y), simplify = FALSE) # true value of test sets
  
  predtest <- Kbyeta_opt <- HP_opt <- index_opt <- coefficients <- 
    vector(mode = "list", length = length(trainIn))
  
  rank_B_all <- ranked_B_all <- vector(mode = "list")
  
  #### innerloop results
  inloop <- function(j){
    
    res_innerloop_LsOPLS <- innerloop_LsOPLS(trainIn = trainIn[[j]], 
                                             testIn =  testIn[[j]], 
                                             nOmax = nOmax, 
                                             eta =  eta,
                                             K =  K,
                                             mc.cores = mc.cores)
    
    res_eval <- innerloop_LsOPLS_eval(res_innerloop_LsOPLS,
                                      mc.cores = mc.cores)
    
    index_opt <- which(rownames(res_eval$byeta_opt) ==
                         as.character(res_eval$HP_opt[,"eta"]))
    
    return(list(HP_opt = res_eval$HP_opt, byeta_opt=res_eval$byeta_opt,
                index_opt = index_opt))
  }
  
  ptm <- proc.time()
  res_inloop <- sapply(X=1:length(trainIn), FUN = inloop)
  proc.time() - ptm
  
  HP_opt <- res_inloop["HP_opt",]
  byeta_opt <- res_inloop["byeta_opt",]
  index_opt <- res_inloop["index_opt",]
  
  
  #### outerloop results
  outloop <- function(i,j){
    nOc_val <- byeta_opt[[j]][i,]["nOc"]
    eta_val <- byeta_opt[[j]][i,]["eta"]
    K_val <-  byeta_opt[[j]][i,]["K"]
    
    #Manual OPLS:
    #Deflate the X matrix:
    OPLSDA_res <- MBXUCL::OPLSDA(x = trainOut[[j]]$x, y = trainOut[[j]]$y, 
                                 impT = FALSE, impG = FALSE, 
                                 no = nOc_val)
    Xopls1 <- OPLSDA_res$Xopls
    Wortho <- OPLSDA_res$Wortho
    Portho <- OPLSDA_res$Portho
    
    # SPLS on the deflated matrix
    res_spls <- spls::spls(x=Xopls1, y=trainOut[[j]]$y, K = K_val, 
                           eta = eta_val, kappa=0.5, 
                           select="pls2",  fit="simpls", scale.x=FALSE,
                           scale.y=FALSE, eps=1e-4, 
                           maxstep=100, trace=FALSE)
    return(list(res_spls=res_spls,Xopls1=Xopls1,Wortho=Wortho, Portho=Portho))
  }
  
  
  
  res_outloop <- sapply(1:length(trainIn), function(j) 
    mclapply(1:dim(byeta_opt[[j]])[1], 
             function(i) outloop(j=j, i=i),
             mc.cores=mc.cores), simplify = FALSE)
  
  res_outloop_spls <-  lapply(1:length(trainIn), 
                              function(j) lapply(1:nrow(byeta_opt[[j]]),
                                                 function(i) res_outloop[[j]][[i]]$res_spls))
  
  Xopls1 <- lapply(1:length(trainIn), 
                   function(j) lapply(1:nrow(byeta_opt[[j]]),
                                      function(i) res_outloop[[j]][[i]]$Xopls1))
  Wortho  <- lapply(1:length(trainIn), 
                    function(j) lapply(1:nrow(byeta_opt[[j]]),
                                       function(i) res_outloop[[j]][[i]]$Wortho))
  Portho <- lapply(1:length(trainIn), 
                   function(j) lapply(1:nrow(byeta_opt[[j]]),
                                      function(i) res_outloop[[j]][[i]]$Portho))
 
  # uncorrected coefficients to be applied to X filtered
  buncorr <- lapply(1:length(trainIn), function(j) sapply(res_outloop_spls[[j]], 
                    function(x) as.vector(predict(x, type = c("coefficient")))))
  
  
  # corrected coefficients to be applied to X raw
  bcorr <- vector(mode="list", length = length(trainIn))
  for (j in 1:length(trainIn)){
    CO <- c()
    for (i in 1:ncol(buncorr[[j]])){
      co <- (diag(1, dim(Xopls1[[j]][[i]])[2], dim(Xopls1[[j]][[i]])[2]) - 
               Wortho[[j]][[i]] %*% solve(t(Portho[[j]][[i]]) %*% Wortho[[j]][[i]]) %*% 
               t(Portho[[j]][[i]])) %*% buncorr[[j]][,i]
      CO <- cbind(CO, co)
    }
    bcorr[[j]] <- t(CO)
  }

  
  # best meta parameters +++++++
  best <- vector(mode = "list", length = 2)
  names(best) <- c("rank_B_all","ranked_B_all")
  best$rank_B_all <- best$ranked_B_all <-  vector(mode = "list", length = 1)
  names(best$rank_B_all) <- names(best$ranked_B_all) <-  "LsOPLS_coef"
  
  best$rank_B_all$LsOPLS_coef <- t(mapply(function(i, coef) coef[,i], 
                                          coef = buncorr, i = index_opt))
  
  best$ranked_B_all$LsOPLS_coef <- sapply(1:dim(best$rank_B_all$LsOPLS_coef)[1],
                                          function(i)
                                            which(abs(buncorr[[i]][,index_opt[[i]]])> 0),
                                          simplify = FALSE)
  
  testOutx <- sapply(testOut, function(x) x$x, simplify = FALSE)
  
  ##
  
  pred_lsopls <- function(xtrain, xtest, bcorr){
    x.new <- xtest - colMeans(xtrain)
    predtest <- x.new %*% bcorr
    predtest
  }
  
  predtest <-  vector(mode="list")
  for (j in 1:length(trainIn)){
    PRED <- vector(mode="list")
    for (i in 1:ncol(buncorr[[j]])){
      PRED[[i]] <- pred_lsopls(trainOutX[[j]], xtest=testOutx[[j]], 
                               bcorr=bcorr[[j]][i,])
    }
    predtest[[j]] <- PRED
  }
  
  best$predtest <- mapply(function(i, pred) as.vector(pred[[i]]),
                          pred = predtest, i = index_opt, SIMPLIFY = FALSE)
  
  # best$predtest <- predtest
  best$truetest <- truetest
  
  
  # all etas +++++++
  buncorr <- sapply(buncorr,t, simplify = FALSE)
  rank_B_all$LsOPLS_coef_eta <- buncorr
  
  ##
  fun_ranked <- function(i,j){
    which(abs(rank_B_all$LsOPLS_coef_eta[[j]][i,]) > 0)
  }
  
  ranked_B_all$LsOPLS_coef_eta <- sapply(1:dim(rank_B_all$LsOPLS_coef_eta[[1]])[1], 
                                         function(x) as.vector(sapply(1:length(trainIn), 
                                                            function(y) fun_ranked(j=y, i=x),
                                                            simplify = FALSE)), simplify = FALSE)
  
  names(ranked_B_all$LsOPLS_coef_eta) <- paste0("eta",rownames(byeta_opt[[1]]))
  
  return(list(predtest = predtest, truetest = truetest, byeta_opt = byeta_opt, 
              rank_B_all = rank_B_all, ranked_B_all = ranked_B_all, 
              trainOut = trainOut, best=best))
}

