\name{predict.splsda}
\alias{predict.splsda}
\alias{coef.splsda}
\title{
Make predictions or extract coefficients from a fitted SPLSDA model
}
\description{
Make predictions or extract coefficients from a fitted SPLSDA object.
}
\usage{
\method{predict}{splsda}( object, newx, type = c("fit","coefficient"),
    fit.type = c("class","response"), ... )
\method{coef}{splsda}( object, ... )
}
\arguments{
    \item{object}{ A fitted SPLSDA object.}
    \item{newx}{
    If \code{type="fit"}, then \code{newx} should be the predictor matrix of test dataset.
    If newx is omitted, then prediction of training dataset is returned.
    If \code{type="coefficient"}, then \code{newx} can be omitted.
    }
    \item{type}{
    If \code{type="fit"}, fitted values are returned.
    If \code{type="coefficient"},
    coefficient estimates of SPLSDA fits are returned.
    }
    \item{fit.type}{
    If \code{fit.type="class"}, fitted classes are returned.
    If \code{fit.type="response"}, fitted probabilities are returned.
    Relevant only when \code{type="fit"}.
    }
    \item{...}{ Any arguments for \code{predict.splsda}
                should work for \code{coef.splsda}. }
}
\details{
Users can input either only selected variables or all variables for \code{newx}.
}
\value{Matrix of coefficient estimates if \code{type="coefficient"}.
Matrix of predicted responses if \code{type="fit"}
(responses will be predicted classes if \code{fit.type="class"}
or predicted probabilities if \code{fit.type="response"}).}
\references{
Chung D and Keles S (2010), 
"Sparse partial least squares classification for high dimensional data",
\emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.}
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{print.splsda}}. }
\examples{
data(prostate)
# SPLSDA with eta=0.8 & 3 hidden components
f <- splsda( prostate$x, prostate$y, K=3, eta=0.8, scale.x=FALSE )
# Print out coefficients
coef.f <- coef(f)
coef.f[ coef.f!=0, ]
# Prediction on the training dataset
(pred.f <- predict( f, type="fit" ))
}
\keyword{models}
\keyword{multivariate}
\keyword{methods}
