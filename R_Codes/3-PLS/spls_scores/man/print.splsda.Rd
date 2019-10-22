\name{print.splsda}
\alias{print.splsda}
\title{Print function for a SPLSDA object}
\description{
 Print out SPLSDA fits, the number and the list of selected predictors.
}
\usage{
\method{print}{splsda}( x, ... )
}
\arguments{
  \item{x}{ A fitted SPLSDA object.}
  \item{...}{ Additonal arguments for generic \code{print}.}
}
\value{NULL.}
\references{
Chung D and Keles S (2010), 
"Sparse partial least squares classification for high dimensional data",
\emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.}
\author{ Dongjun Chung and Sunduz Keles. }
\seealso{ \code{\link{predict.splsda}} and \code{\link{coef.splsda}}. }
\examples{
data(prostate)
# SPLSDA with eta=0.8 & 3 hidden components
f <- splsda( prostate$x, prostate$y, K=3, eta=0.8, scale.x=FALSE )
print(f)
}
\keyword{models}
\keyword{multivariate}
\keyword{methods}
