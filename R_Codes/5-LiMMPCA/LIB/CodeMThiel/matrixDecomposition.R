##############################################################################
# Project: ASCA+
# 
# Description:
#
# R-script calling the function to decompose a matrix in ASCA+ such as in Thiel et al. (2017)
#
# Author: mthiel & rvanoirbeek
#
# Maintainer: mthiel
#
# Date: 25/04/17
#
# Changes: Resolving the code if there is no factor variable
#
# Version: 1.1
#
###############################################################################

matrixDecomposition <- function(formula, outcomes, design){
  
  if(missing(outcomes)) outcomes <- environment(formula)
  
  #Checking formula
    formulaChar <- as.character(formula)
  if(length(formulaChar) == 3){
    formulaDesignMatrix <- as.formula(paste(formulaChar[1], formulaChar[3]))
  } else if (length(formulaChar) == 2){
    formulaDesignMatrix <- formula
  } else {
    stop("Please put the formula argument in its right form")
  }
  
  #Checking correspondance between formula names and design names
    varNames <- all.vars(formulaDesignMatrix)
  matchesVarNames <- varNames %in% names(design) 
  if(!all(matchesVarNames, na.rm = FALSE)){
    stop("Some of the variable names, present in the formula argument, do not correspond to one of the column names of the design argument. Please adapt either one of both arguments.")
  }  
  
  #Checking which variables are factors
  factorsDesign <- names(Filter(is.factor, design))
  varNamesFactors <- intersect(factorsDesign,varNames)  
  
  #Creating model matrix
  #If factors are present, a list is created to specify which variables are considered as factors in model.matrix
  if(length(varNamesFactors)!=0)
  {
  contrasts.arg.Values <- list()
  length(contrasts.arg.Values) <- length(varNamesFactors)
  names(contrasts.arg.Values) <- varNamesFactors
  for(iList in 1:length(contrasts.arg.Values)) contrasts.arg.Values[[iList]] <- "contr.sum"
  modelMatrix <- (model.matrix(formulaDesignMatrix, contrasts.arg = contrasts.arg.Values, data = design))
  }
  
  #If factors are not present
  if(length(varNamesFactors)==0)
  {
    modelMatrix <- (model.matrix(formulaDesignMatrix, data = design))
  }

  #Creating a list containing effect matrices and model matrices by effect
  dummyVarNames <- colnames(modelMatrix)
  presencePolynomialEffects <- str_detect(dummyVarNames, '\\^[0-9]')
  covariateEffectsNames <- character(length = length(dummyVarNames))
  covariateEffectsNames[presencePolynomialEffects] <- dummyVarNames[presencePolynomialEffects]
  covariateEffectsNames[!presencePolynomialEffects] <- gsub('[0-9]', '', dummyVarNames[!presencePolynomialEffects])
  covariateEffectsNames[covariateEffectsNames == '(Intercept)'] <- 'Intercept'
  covariateEffectsNamesUnique <- unique(covariateEffectsNames)
  nEffect <- length(covariateEffectsNamesUnique)
 
  #effects matrices
  effectMatrices <- list()
  length(effectMatrices) <- nEffect
  names(effectMatrices) <- covariateEffectsNamesUnique
  
  #model matrices by effect
  modelMatrixByEffect <- list()
  length(modelMatrixByEffect) <- nEffect
  names(modelMatrixByEffect) <- covariateEffectsNamesUnique
  
  #GLM decomposition calculated by using glm.fit and alply on outcomes
  resGLM <- alply(outcomes, 2, function(xx) glm.fit(modelMatrix, xx))
  parameters <- t(laply(resGLM, function(xx) xx$coefficients)) 
  predictedValues <- t(laply(resGLM, function(xx) xx$fitted.values)) 
  residuals <- t(laply(resGLM, function(xx) xx$residuals)) 
 
  #List with type 3 residuals
  Type3Residuals <- list()
  length(Type3Residuals) <- nEffect
  names(Type3Residuals) <- covariateEffectsNamesUnique
  
  #List with Frobenius norms
  matrixVolume <- list()
  length(matrixVolume) <- nEffect + 2
  names(matrixVolume) <- c(covariateEffectsNamesUnique, "outcomes", "residuals")
  
  #List with variation percentages
  variationPercentages <- list()
  length(variationPercentages) <- nEffect
  names(variationPercentages) <- c(covariateEffectsNamesUnique[covariateEffectsNamesUnique != 'Intercept'], 'residuals')
  
  
  for(iEffect in 1:nEffect){
    selection <- which(covariateEffectsNames == covariateEffectsNamesUnique[iEffect])
    selectionComplement <- which(covariateEffectsNames != covariateEffectsNamesUnique[iEffect])
    #Effect matrices
    effectMatrices[[iEffect]] <- t(aaply(parameters, 2, function(xx) as.matrix(modelMatrix[, selection])%*%xx[selection]))
    #Model matrices by effect
    modelMatrixByEffect[[iEffect]] <- as.matrix(modelMatrix[, selection])
    
    matrixVolume[[iEffect]] <- (norm(effectMatrices[[iEffect]], "F"))^2
    if(covariateEffectsNamesUnique[iEffect] != 'Intercept'){
      resGLMComplement <- alply(outcomes, 2, function(xx) glm.fit(modelMatrix[, selectionComplement], xx))
      Type3Residuals[[iEffect]] <- t(laply(resGLMComplement, function(xx) xx$residuals))
    }
  }
  
  matrixVolume[[nEffect + 1]] <- (norm(outcomes, "F"))^2
  matrixVolume[[nEffect + 2]] <- (norm(residuals, "F"))^2
  
  denominatorSSType3 <- norm(outcomes - effectMatrices[['Intercept']], "F")^2
  numeratorFullModelSSType3 <- norm(residuals, "F")^2
  variationPercentages[1:nEffect-1] <- llply(Type3Residuals[-1], function(xx) 100*(norm(xx, "F")^2 - numeratorFullModelSSType3)/denominatorSSType3)  
  #Variation percentages
  variationPercentages[[nEffect]] <- 100*numeratorFullModelSSType3/denominatorSSType3
  
  return(list(effectMatrices = effectMatrices, residuals = residuals, predictedValues=predictedValues, parameters=parameters, variationPercentages = variationPercentages, Type3Residuals=Type3Residuals, modelMatrixByEffect = modelMatrixByEffect, modelMatrix=modelMatrix))
  
}
