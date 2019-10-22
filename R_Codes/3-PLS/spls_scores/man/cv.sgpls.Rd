\name{cv.sgpls}
\alias{cv.sgpls}
\title{Compute and plot the cross-validated error for SGPLS classification}
\description{
 Draw heatmap of
 v-fold cross-validated misclassification rates
 and return optimal eta (thresholding parameter)
 and K (number of hidden components).
}
\usage{
cv.sgpls( x, y, fold=10, K, eta, scale.x=TRUE, plot.it=TRUE,
        br=TRUE, ftype='iden', n.core=8 )
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector of class indices. }
  \item{fold}{ Number of cross-validation folds. Default is 10-folds. }
  \item{K}{ Number of hidden components. }
  \item{eta}{ Thresholding parameter. \code{eta} should be between 0 and 1. }
  \item{scale.x}{ Scale predictors by dividing each predictor variable
                by its sample standard deviation? }
  \item{plot.it}{ Draw the heatmap of cross-validated misclassification rates? }
  \item{br}{ Apply Firth's bias reduction procedure? }
  \item{ftype}{ Type of Firth's bias reduction procedure.
     Alternatives are \code{"iden"} (the approximated version)
     or \code{"hat"} (the original version).
     Default is \code{"iden"}.}
  \item{n.core}{ Number of CPUs to be used when parallel computing is utilized. }
}
\details{
Parallel computing can be utilized for faster computation.
Users can change the number of CPUs to be used
by changing the argument \code{n.core}.
}
\value{
Invisibly returns a list with components:
  \item{err.mat}{ Matrix of cross-validated misclassification rates.
                Rows correspond to \code{eta} and
                columns correspond to number of components (\code{K}). }
  \item{eta.opt}{ Optimal \code{eta}. }
  \item{K.opt}{ Optimal \code{K}. }
}
\references{
    Chung D and Keles S (2010), 
    "Sparse partial least squares classification for high dimensional data",
    \emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.
}
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{print.sgpls}}, \code{\link{predict.sgpls}},
and \code{\link{coef.sgpls}}. }
\examples{
    data(prostate)
    set.seed(1)

    # misclassification rate plot. eta is searched between 0.1 and 0.9 and
    # number of hidden components is searched between 1 and 5
    \dontrun{
        cv <- cv.sgpls(prostate$x, prostate$y, K = c(1:5), eta = seq(0.1,0.9,0.1),
            scale.x=FALSE, fold=5)
    }
    
    \dontshow{
        cv <- list(); cv$eta.opt <- 0.6; cv$K.opt <- 3; 
    }

    (sgpls(prostate$x, prostate$y, eta=cv$eta.opt, K=cv$K.opt, scale.x=FALSE))
}
\keyword{models}
\keyword{multivariate}
