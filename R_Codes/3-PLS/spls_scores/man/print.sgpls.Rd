\name{print.sgpls}
\alias{print.sgpls}
\title{Print function for a SGPLS object}
\description{
 Print out SGPLS fit, the number and the list of selected predictors.
}
\usage{
\method{print}{sgpls}( x, ... )
}
\arguments{
  \item{x}{ A fitted SGPLS object.}
  \item{...}{ Additonal arguments for generic \code{print}.}
}
\value{NULL.}
\references{
Chung D and Keles S (2010), 
"Sparse partial least squares classification for high dimensional data",
\emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.}
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{predict.sgpls}} and \code{\link{coef.sgpls}}. }
\examples{
data(prostate)
# SGPLS with eta=0.55 & 3 hidden components
f <- sgpls( prostate$x, prostate$y, K=3, eta=0.55, scale.x=FALSE )
print(f)
}
\keyword{models}
\keyword{multivariate}
\keyword{methods}
