\name{plot.spls}
\alias{plot.spls}
\title{Plot the coefficient path of SPLS regression}
\description{
 Provide the coefficient path plot of SPLS regression as a function of
 the number of hidden components (K) when eta is fixed.
}
\usage{
\method{plot}{spls}( x, yvar=c(1:ncol(x$y)), ... )
}
\arguments{
  \item{x}{ A fitted SPLS object. }
  \item{yvar}{ Index vector of responses to be plotted.}
  \item{...}{ Other parameters to be passed through to generic \code{plot}.}
}
\details{
\code{plot.spls} provides the coefficient path plot of SPLS fits.
The plot shows how estimated coefficients change
as a function of the number of hidden components (\code{K}),
when \code{eta} is fixed at the value used by the original SPLS fit.
}
\value{NULL.}
\references{
Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection",
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25.}
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ \code{\link{print.spls}}, \code{\link{predict.spls}},
and \code{\link{coef.spls}}. }
\examples{
data(yeast)
# SPLS with eta=0.7 & 8 hidden components
f <- spls( yeast$x, yeast$y, K=8, eta=0.7 )
# Draw coefficient path plots for the first two responses
plot( f, yvar=c(1:2) )
}
\keyword{hplot}
\keyword{methods}
