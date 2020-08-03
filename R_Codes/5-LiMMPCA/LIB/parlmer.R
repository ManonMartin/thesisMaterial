#' Parallel mixed model
#'
#' Parallel mixed modelling using the lme4::lmer function
#'
#' @param design design data frame with rownames as ID
#' @param outcomes response matrix (eg spectral matrix in NMR) with rownames as ID
#' @param form mixed models formula
#' @param REML if TRUE, use REML to fit the models
#'
#' @return a list with the mkMerMod output (cf. ?lme4:lmer) for each response vector.
#'
#' @details the following steps are applied: (1) lme4::lFormula, (2) lme4::mkLmerDevfun ,
#' (3) lme4::optimizeLmer, (4) lme4::mkMerMod.
#'
#' @import lme4
#'
#' @export



parlmer <- function(design, outcomes, form, REML){

  # -- checks
  if (!is.logical(REML)) {
    stop("REML is not logical")
  }

  if (!is.data.frame(design)) {
    stop("design is not a data frame")
  }

  if (!is.matrix(outcomes) | !is.numeric(outcomes)) {
    stop("outcomes is not a matrix and/or not numeric ")
  }

  # check formula formatting
  formule <- as.formula(form)

  # Checking correspondance between formula names and design names
  allvars_formule <- all.vars(formule)
  matchesVarNames <- allvars_formule %in% names(design)

  # FIXME
  # if(!all(matchesVarNames, na.rm = FALSE)){
  #   stop("Some of the variable names, present in the formula argument, do not correspond to one of the column names of the design argument. Please adapt either one of both arguments.")
  # }

  # checking that the design variables are factors
  checkfactors <- sapply(design[,allvars_formule],is.factor)
  if (!all(checkfactors)) {
    stop("one or several variables in the design is not a factor.")
  }

  # ensure that "design" and "outcomes" are identically ordered by rownames
  if (!isTRUE(all.equal(sort(rownames(design)),sort(rownames(outcomes))))){
    stop("not a perfect match between the rownames of design and outcomes")
  }



  # design <- design[sort(rownames(design)),]
  # outcomes <- outcomes[sort(rownames(outcomes)),]
  m <- ncol(outcomes)
  respnames <- colnames(outcomes)


  total.data <- cbind(design, outcomes)


  # -- lFormula
  # "lFormula takes the arguments that would normally be passed
  # to lmer, checking for errors and processing the formula and
  # data input to create a list of objects required to fit a mixed model."


  fmla <- sapply(paste0(colnames(outcomes), form), as.formula)

  # lFormula to get the names of fixed and random effects
  # lFormula_res_lapply  <- lapply(fmla, lFormula, data = total.data, REML = REML)

  # res <- lFormula(formula=paste(colnames(outcomes)[1],form), data = total.data)
  #
  # randNames <- names(res$reTrms$flist)
  fixNames <- attr(terms.formula(as.formula(form)),"term.labels")[!grepl(" | ",
                   attr(terms.formula(as.formula(form)),"term.labels"))]

####

  if (length(fixNames)>0){ # if fixed effects are present
    contrasts.arg.Values <- list()
    length(contrasts.arg.Values) <- length(fixNames[!grepl(":", fixNames)])# keeps only the main fixed effects in case of interactions
    names(contrasts.arg.Values) <- fixNames[!grepl(":", fixNames)]

    for(iList in 1:length(contrasts.arg.Values)) contrasts.arg.Values[[iList]] <- "contr.sum"
    lFormula_res_lapply  <- lapply(fmla, lFormula, data = total.data, REML = REML,
                                   contrasts = contrasts.arg.Values)
   }else{
     lFormula_res_lapply  <- lapply(fmla, lFormula, data = total.data, REML = REML)
   }

  # FIXME : tester tous les cas de formules possibles: interaction, fixNames hierarchique ...
  # random and fixed effects names in the design
  ranNames <- names(lFormula_res_lapply[[1]]$reTrms$flist)
  indexpunct <- gregexpr("[[:punct:]]", ranNames)
  indexpunct <- sapply(indexpunct, function(x) x[length(x)])+1
  # indexpunct[indexpunct==-1] <- nchar(ranNames[indexpunct==-1])+1
  indexpunct[indexpunct==0] <- 1
  ranNames <- substr(ranNames,indexpunct,nchar(ranNames))

  # A model frame containing the variables needed to create an lmerResp or glmResp instance.
  fr_list <- lapply(lFormula_res_lapply, function(x) x[["fr"]])

  # information on random effects structure
  reTrms_list <- lapply(lFormula_res_lapply, function(x) x[["reTrms"]])
  # as.matrix(t(reTrms_list[[1]]$Zt)) # design matrix for random effects

  RanModMatlist <- reTrms_list[[1]]$Ztlist
  RanModMatlist <- lapply(RanModMatlist, function(x) t(as.matrix(x)))

  # FixModMatlist <-

  # TO BE FIXED : ncol names of random model matrix???

  # colnam <- reTrms_list[[1]]$Lind
  # NamesRand <- names(reTrms_list[[1]]$flist)
  # for (i in 1:max(reTrms_list[[1]]$Lind)){
  #   colnam[colnam==i] <- NamesRand[i]
  # }
  #
  # names(reTrms_list[[1]]$Ztlist)
  #
  # rownames(reTrms_list[[1]]$Zt) <- paste0(colnam,
  #        rownames(reTrms_list[[1]]$Zt))


  # fixed-effects design matrix
  # Note: to be designed manually if necessary
  X_list <- lapply(lFormula_res_lapply, function(x) x[["X"]])
  # X_list[[1]]

  # lFormula_res <- lFormula(formula = fmla1, data = total.data, REML = REML)

  # create a list with the model frames fr for all the response vectors
  # fr <- vector(mode = "list")
  # for (i in 1:m) {
  #   fr[[i]] <- lFormula_res$fr
  #   fr[[i]][,1] <- scores_spectra[,i]
  #   colnames(fr[[i]])[1] <- colnames(outcomes)[i]
  # }


  #---- mkLmerDevfun ---- creates a deviance function
  mkLmerDevfun_mapply_res <- mapply(mkLmerDevfun, fr = fr_list,
                                      X = X_list,
                                      reTrms = reTrms_list,
                                    MoreArg = list(REML = REML), SIMPLIFY = FALSE)


  # mkLmerDevfun_mapply_res <- list()
  # for (i in 1:m) {
  #   mkLmerDevfun_mapply_res[[i]] <- mkLmerDevfun(fr = lFormula_res_lapply[[i]]$fr,
  #                                                X = lFormula_res_lapply[[i]]$X,
  #                                                reTrms = lFormula_res_lapply[[i]]$reTrms,
  #                                                REML = lFormula_res_lapply[[i]]$REML)
  #
  # }

  #---- optimizeLmer ---- takes a deviance function and optimizes over theta
  optimizeLmer_sapply_res <- lapply(mkLmerDevfun_mapply_res, optimizeLmer)


  #---- mkMerMod ---- takes the environment of a deviance function,
  # the results of an optimization, a list of random-effect terms,
  # a model frame, and a model all and produces a [g]lmerMod object

  rho_list <- lapply(mkLmerDevfun_mapply_res, environment)
  res <- mapply(mkMerMod, rho = rho_list,
                opt = optimizeLmer_sapply_res,
                fr = fr_list, reTrms = reTrms_list)

  names(res) <- respnames

  # FixedModMatlist

  fixed_form_term.labels <- attr(terms(formula(res[[1]], fixed.only = TRUE)),"term.labels")

  if (length(fixed_form_term.labels)>0){
    fix_lev_Names <- colnames(model.matrix(res$PC1,type = c("fixed")))[!colnames(model.matrix(res$PC1,type = c("fixed")))%in%"(Intercept)"]
    sign <- c()
    for (i in 1:length(fixed_form_term.labels[!grepl(":", fix_lev_Names)])){
      sign <- rbind(sign, as.integer(grepl(fixed_form_term.labels[i], fix_lev_Names)))
    }
    sign <- apply(sign,2,paste, collapse = "")
    names(sign) <- fix_lev_Names
    uniqSign <- unique(sign)


    FixedModMatlist <- vector(mode = "list", length = length(fixed_form_term.labels))
    names(FixedModMatlist) <- fixed_form_term.labels

    mat <- as.matrix(model.matrix(res[[1]],type = c("fixed")))
    for (i in 1:length(fixed_form_term.labels)){
      id <- names(sign)[sign==uniqSign[i]]
      FixedModMatlist[[i]] <- matrix(mat[,id], dimnames = list(NULL, id), ncol = length(id))
    }
    if ("(Intercept)"%in%colnames(model.matrix(res[[1]],type = c("fixed")))){
      FixedModMatlist <- append(list(matrix(model.matrix(res[[1]],type = c("fixed"))[,"(Intercept)"], ncol = 1,
                                            dimnames = list(NULL, "(Intercept)"))), FixedModMatlist)
      # FixedModMatlist$"(Intercept)" = matrix(model.matrix(res[[1]],type = c("fixed"))[,"(Intercept)"], ncol = 1,
      #                                        dimnames = list(NULL, "(Intercept)"))
      names(FixedModMatlist)[1] <- "(Intercept)"
    }
  } else{
      fixNames <- NULL
      if ("(Intercept)"%in%colnames(model.matrix(res[[1]],type = c("fixed")))){
        FixedModMatlist <- vector(mode = "list", length = 1)
        FixedModMatlist[[1]] <- matrix(model.matrix(res[[1]],type = c("fixed"))[,"(Intercept)"], ncol = 1,
               dimnames = list(NULL, "(Intercept)"))
        names(FixedModMatlist)[1] <- "(Intercept)"

      } else{
        FixedModMatlist <- NULL
      }

  }


  # append(list(FixedModMatlist$"(Intercept)"), FixedModMatlist)
  return(list(merMod_obj = res, RanModMatlist = RanModMatlist, FixedModMatlist=FixedModMatlist, ranNames= ranNames,
              fixNames=fixNames, arg = list(design=design, outcomes=outcomes,
                                           form=form, REML=REML)))

}

