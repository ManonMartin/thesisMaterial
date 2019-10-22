\name{print.spls}
\alias{print.spls}
\title{Print function for a SPLS object}
\description{
 Print out SPLS fit, the number and the list of selected predictors.
}
\usage{
\method{print}{spls}( x, ... )
}
\arguments{
  \item{x}{ A fitted SPLS object.}
  \item{...}{ Additonal arguments for generic \code{print}.}
}
\value{NULL.}
\references{
Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection,"
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25.}
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ \code{\link{plot.spls}}, \code{\link{predict.spls}},
and \code{\link{coef.spls}}. }
\examples{
data(yeast)
# SPLS with eta=0.7 & 8 hidden components
f <- spls( yeast$x, yeast$y, K=8, eta=0.7 )
print(f)
}
\keyword{regression}
\keyword{multivariate}
\keyword{methods}
