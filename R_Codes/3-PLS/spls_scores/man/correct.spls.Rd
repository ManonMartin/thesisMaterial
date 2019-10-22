\name{correct.spls}
\alias{correct.spls}
\title{Correct the initial SPLS coefficient estimates
        based on bootstrapped confidence intervals}
\description{
 Correct initial SPLS coefficient estimates of the selected predictors
 based on bootstrapped confidence intervals and draw
 heatmap of original and corrected coefficient estimates.
}
\usage{
correct.spls( object, plot.it=TRUE )
}
\arguments{
  \item{object}{ An object obtained from the function \code{ci.spls}. }
  \item{plot.it}{ Draw the heatmap of
                    original coefficient estimates
                    and corrected coefficient estimates? }
}
\details{
The set of the selected variables is updated by setting
the coefficients with zero-containing confidence intervals to zero.
}
\value{ Invisibly returns a matrix of corrected coefficient estimates.
}
\references{
Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection",
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25.}
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ \code{\link{ci.spls}}. }
\examples{
data(mice)
# SPLS with eta=0.6 & 1 latent components
f <- spls( mice$x, mice$y, K=1, eta=0.6 )
# Calculate confidence intervals of coefficients
ci.f <- ci.spls(f)
# Corrected coefficient estimates
cf <- correct.spls( ci.f )
cf[20,1:5]
}
\keyword{regression}
\keyword{multivariate}
