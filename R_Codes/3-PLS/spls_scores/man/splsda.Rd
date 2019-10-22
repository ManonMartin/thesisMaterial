\name{splsda}
\alias{splsda}
\title{Fit SPLSDA classification models}
\description{
 Fit a SPLSDA classification model.
}
\usage{
splsda( x, y, K, eta, kappa=0.5,
    classifier=c('lda','logistic'), scale.x=TRUE, ... )
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector of class indices. }
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
  \item{...}{ Other parameters to be passed through to \code{spls}.}
}
\details{
The SPLSDA method is described in detail in Chung and Keles (2010).
SPLSDA provides a two-stage approach for PLS-based classification with variable selection,
by directly imposing sparsity on the dimension reduction step of PLS
using sparse partial least squares (SPLS) proposed in Chun and Keles (2010).
\code{y} is assumed to have numerical values, 0, 1, ..., G,
where G is the number of classes subtracted by one.
The option \code{classifier} refers to the classifier used in the second step of SPLSDA
and \code{splsda} utilizes algorithms offered by \pkg{MASS} and \pkg{nnet} packages
for this purpose.
If \code{classifier="logistic"}, then either logistic regression or multinomial regression is used.
Linear discriminant analysis (LDA) is used if \code{classifier="lda"}.
\code{splsda} also utilizes algorithms offered by the \pkg{pls} package for fitting \code{spls}.
The user should install \pkg{pls}, \pkg{MASS} and \pkg{nnet} packages before using \code{splsda} functions.
}
\value{A \code{splsda} object is returned.
print, predict, coef methods use this object.}
\references{
Chung D and Keles S (2010), 
"Sparse partial least squares classification for high dimensional data",
\emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.

Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection",
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25. }
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{print.splsda}}, \code{\link{predict.splsda}}, and \code{\link{coef.splsda}}. }
\examples{
data(prostate)
# SPLSDA with eta=0.8 & 3 hidden components
f <- splsda( prostate$x, prostate$y, K=3, eta=0.8, scale.x=FALSE )
print(f)
# Print out coefficients
coef.f <- coef(f)
coef.f[ coef.f!=0, ]
}
\keyword{models}
\keyword{multivariate}
