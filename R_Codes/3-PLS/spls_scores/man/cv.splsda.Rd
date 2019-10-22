\name{cv.splsda}
\alias{cv.splsda}
\title{Compute and plot cross-validated error for SPLSDA classification}
\description{
 Draw heatmap of
 v-fold cross-validated misclassification rates
 and return optimal eta (thresholding parameter)
 and K (number of hidden components).
}
\usage{
cv.splsda( x, y, fold=10, K, eta, kappa=0.5,
        classifier=c('lda','logistic'), scale.x=TRUE, plot.it=TRUE, n.core=8 )
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector of class indices. }
  \item{fold}{ Number of cross-validation folds. Default is 10-folds. }
  \item{K}{ Number of hidden components. }
  \item{eta}{ Thresholding parameter. \code{eta} should be between 0 and 1. }
  \item{kappa}{ Parameter to control the effect of
     the concavity of the objective function
     and the closeness of original and surrogate direction vectors.
     \code{kappa} is relevant only for multicategory classification.
     \code{kappa} should be between 0 and 0.5. Default is 0.5. }
  \item{classifier}{ Classifier used in the second step of SPLSDA.
     Alternatives are \code{"logistic"} or \code{"lda"}.
     Default is \code{"lda"}.}
  \item{scale.x}{ Scale predictors by dividing each predictor variable
                by its sample standard deviation? }
  \item{plot.it}{ Draw the heatmap of the cross-validated misclassification rates? }
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
\emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.}
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{print.splsda}}, \code{\link{predict.splsda}},
and \code{\link{coef.splsda}}. }
\examples{
data(prostate)
set.seed(1)
# misclassification rate plot. eta is searched between 0.1 and 0.9 and
# number of hidden components is searched between 1 and 5
\dontrun{ cv <- cv.splsda( prostate$x, prostate$y, K = c(1:5), eta = seq(0.1,0.9,0.1),
         scale.x=FALSE, fold=5 )}
\dontshow{ cv <- list(); cv$eta.opt <- 0.8; cv$K.opt <- 3; }
(splsda( prostate$x, prostate$y, eta=cv$eta.opt, K=cv$K.opt, scale.x=FALSE ))
}
\keyword{models}
\keyword{multivariate}
