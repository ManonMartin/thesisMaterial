#' Permutation test without any nested effects. 
#' 
#' By means of a formula, the P-value of all possible model effects are calculated. 
#' @param formula Formula object in the same style as used for the well-known glm function. It is important that the variable names used for the formula argument stem from the column names of the 'design' argument of this function. 
#' @param outcomes Numerical matrix consisting of the outcome values (columns) for each different experiment (rows).
#' @param design Data frame in which variables (columns) that express the experimental conditions are linked to the actual experiments (rows). Note that it is of capital importance that the column names are well specified, as they will be used to address the variables of the experimental conditions throughout the functions of this R package.
#' @param nPermutations Numerical vector of length 1 indicating how many permutations are performed to calculate the P-value. The default value corresponds to 1000.
#' @return A numerical matrix with the P-values of all the possible model effects is returned.

meanCentering <- function(inputMatrix) {
  VectorOfOnes = rep(1, nrow(inputMatrix))
  inputMatrixColumnMeans = VectorOfOnes %*% t(colMeans(inputMatrix))
  return(inputMatrix - inputMatrixColumnMeans)
}

decompositionPCA <- function(inputMatrix) {
  
  # column centering of inputMatrix
  centeredMatrix <- meanCentering(inputMatrix)
  
  # SVD of inputMatrix
  svdCenteredMatrix = svd(centeredMatrix) # Compute the singular-value decomposition 
  scoreMatrix = svdCenteredMatrix$u %*% diag(svdCenteredMatrix$d) # scores
  normalizedScoreMatrix = svdCenteredMatrix$u # normalised scores
  
  loadingMatrix = svdCenteredMatrix$v # loadings
  singularValues = svdCenteredMatrix$d # singular values
  
  # Variance explained
  explainedVariance = singularValues^2 / (length(singularValues) - 1)
  totalExplainedVariance = sum(explainedVariance)
  relativeExplainedVariance = explainedVariance/totalExplainedVariance
  
  explainedVarianceRelative = 100 * round(relativeExplainedVariance, digits = 3) # variance
  explainedVarianceRelativeCumulative = cumsum(explainedVarianceRelative) # cumulative variance
  
  #return(list(pcs = scoreMatrix, pcu = normalizedScoreMatrix, pcv = loadingMatrix, pcd = singularValues, var = explainedVarianceRelative, cumvar = explainedVarianceRelativeCumulative))
  return(list(scores = scoreMatrix, normalizedScores = normalizedScoreMatrix, loadings = loadingMatrix, singularValues = singularValues, Variances = explainedVarianceRelative, VariancesCumulative = explainedVarianceRelativeCumulative))
}

permutationTest <- function(formula, outcomes, design, nPermutations = 1000){
  
  #Checking if the formula is valid
  
  formulaChar <- as.character(formula)
  if(length(formulaChar) == 3){
    formulaDesignMatrix <- as.formula(paste(formulaChar[1], formulaChar[3]))
  } else if (length(formulaChar) == 2){
    formulaDesignMatrix <- formula
  } else {
    stop("Please put the formula argument in its right form")
  }
  
  varNames <- all.vars(formulaDesignMatrix)
  matchesVarNames <- varNames %in% names(design)
  if(!all(matchesVarNames, na.rm = FALSE)){
    stop("Some of the variable names, present in the formula argument, do not correspond to one of the column names of the design argument. Please adapt either one of both arguments.")
  }  
  
  resultDecompEffectMatrix <- matrixDecomposition(formula, outcomes, design)
  namesDecomposition <- names(resultDecompEffectMatrix$effectMatrices)
  selectionEffectMatrices <- which(namesDecomposition != "Intercept")
  nEffectMatrices <- length(selectionEffectMatrices)
  matrixVolumeOriginal <- matrix(NA, nrow = nEffectMatrices, ncol = 1)
  
  if(nEffectMatrices!= 0){
    for(iEffect in 1:nEffectMatrices){
      selectedMatrix <- resultDecompEffectMatrix$effectMatrices[[selectionEffectMatrices[iEffect]]]
      selectedMatrixDecomposed <- decompositionPCA(selectedMatrix)$scores[, 1:2]
      matrixVolumeOriginal[iEffect,] <- norm(selectedMatrixDecomposed, "F")^2 
    }
  } else {
    stop("No effect matrices were constructed after the decomposition. Please verify that the ensemble of the inputparameters is correct.")
  }
  
  matrixVolumePermuted <- matrix(NA, nrow = nEffectMatrices, ncol = nPermutations)
  varNamesEffects <- names(resultDecompEffectMatrix$effectMatrices)[selectionEffectMatrices]
  rownames(matrixVolumePermuted) <- varNamesEffects
  rownames(matrixVolumeOriginal) <- varNamesEffects
  selectedMatrixDecomposed <- matrix(NA, nrow = dim(outcomes)[1], ncol = dim(outcomes)[2])
  
  for(iPerm in 1:nPermutations){
    outcomesPermuted <- outcomes[sample(nrow(outcomes)),]
    resultDecompEffectMatrix <- matrixDecomposition(formula, outcomesPermuted, design)
    for(iEffect in 1:nEffectMatrices){
      selectedMatrix <- resultDecompEffectMatrix$effectMatrices[[selectionEffectMatrices[iEffect]]]
      selectedMatrixDecomposed <- decompositionPCA(selectedMatrix)$scores[, 1:2]
      matrixVolumePermuted[iEffect, iPerm] <- norm(selectedMatrixDecomposed, "F")^2 
    }
  }
  pValues <- matrix(NA, nrow = nEffectMatrices, ncol = 1)
  rownames(pValues) <- varNamesEffects
  for(iEffect in 1:nEffectMatrices){
    pValues[iEffect, ] <- sum(matrixVolumeOriginal[iEffect] <= matrixVolumePermuted[iEffect, ])/nPermutations
  }
  
  list(pValues = pValues, matrixVolumePermuted = matrixVolumePermuted, matrixVolumeOriginal = matrixVolumeOriginal)
}