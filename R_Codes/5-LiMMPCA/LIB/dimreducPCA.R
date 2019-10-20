#' Dimension reduction by PCA
#'
#' Dimension reduction to replace the data matrix by the m first PCA scores and loadings.
#'
#'
#' @param data numeric matrix
#' @param pcvar percentage of variance to decide on the number of scores and loadings to keep
#'
#' @return a list with the matrices of m PCA scores and loadings.
#'
#'
#' @export

dimreducPCA <- function(data, pcvar = 99) {

  if (!is.matrix(data) | !is.numeric(data)) {
    stop("data is not a matrix and/or not numeric")
  }

  if (!is.numeric(pcvar) | length(pcvar)!=1) {
    stop("pcvar is not numeric and/or has a length < 1")
  }

  if( pcvar > 100 | pcvar < 1) {
    stop("pcvar should be between 1 and 100")
  }

  # centering the variables
  data <- data - matrix(apply(data, 2, mean),
                                                nrow = dim(data)[1],
                                                ncol = dim(data)[2],
                                                byrow = TRUE)
  # apply svd
  res.pca_svd <- svd(data)

  # PCA scores
  res.pca_scores <- res.pca_svd$u %*% diag(res.pca_svd$d)

  # PCA loadings
  res.pca_loadings <- res.pca_svd$v
  res.pca_loadings <- as.matrix(res.pca_loadings)
  # res.pca_loadings <- t(res.pca_loadings)

  # measure the cumulative variance
  res.pca_singularval <- res.pca_svd$d
  names(res.pca_singularval) <- paste0("PC", 1:length(res.pca_singularval))
  res.pca_vars <- res.pca_singularval^2/(nrow(data) - 1)
  res.pca_totalvar <- sum(res.pca_vars)
  res.pca_relvars <- res.pca_vars/res.pca_totalvar
  res.pca_variances <- 100 * res.pca_relvars
  res.pca_cumvariances <- cumsum(res.pca_variances)

  # keep m scores with cum var >= 99%
  m <- which(res.pca_cumvariances >= pcvar)[1]
  res.pca_scores <- res.pca_scores[,1:m]
  res.pca_loadings <- res.pca_loadings[,1:m]

  colnames(res.pca_scores) <- paste0("PC", 1:m)
  rownames(res.pca_scores) <- rownames(data)
  rownames(res.pca_loadings) <- colnames(data)
  colnames(res.pca_loadings) <- paste0("PC", 1:m)

  return(list(pca_scores = res.pca_scores, pca_loadings = res.pca_loadings))
}
