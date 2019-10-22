\name{predict.spls}
\alias{predict.spls}
\alias{coef.spls}
\title{
Make predictions or extract coefficients from a fitted SPLS model
}
\description{
Make predictions or extract coefficients from a fitted SPLS object.
}
\usage{
\method{predict}{spls}( object, newx, type = c("fit","coefficient"), ... )
\method{coef}{spls}( object, ... )
}
\arguments{
    \item{object}{ A fitted SPLS object.}
    \item{newx}{
    If \code{type="fit"}, then \code{newx} should be the predictor matrix of test dataset.
    If newx is omitted, then prediction of training dataset is returned.
    If \code{type="coefficient"}, then \code{newx} can be omitted.
    }
    \item{type}{
    If \code{type="fit"}, fitted values are returned.
    If \code{type="coefficient"},
    coefficient estimates of SPLS fits are returned.
    }
    \item{...}{ Any arguments for \code{predict.spls}
                should work for \code{coef.spls}. }
}
\details{
Users can input either only selected variables or all variables for \code{newx}.
}
\value{Matrix of coefficient estimates if \code{type="coefficient"}.
Matrix of predicted responses if \code{type="fit"}.}
\references{
Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection",
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25.}
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ \code{\link{plot.spls}} and \code{\link{print.spls}}. }
\examples{
data(yeast)
# SPLS with eta=0.7 & 8 latent components
f <- spls( yeast$x, yeast$y, K=8, eta=0.7 )
# Coefficient estimates of the SPLS fit
coef.f <- coef(f)
coef.f[1:5,]
# Prediction on the training dataset
pred.f <- predict( f, type="fit" )
pred.f[1:5,]
}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
