\name{ci.spls}
\alias{ci.spls}
\title{Calculate bootstrapped confidence intervals of SPLS coefficients}
\description{
 Calculate bootstrapped confidence intervals of coefficients
 of the selected predictors and generate confidence interval plots.
}
\usage{
ci.spls( object, coverage=0.95, B=1000,
        plot.it=FALSE, plot.fix="y",
        plot.var=NA, K=object$K, fit=object$fit )
}
\arguments{
  \item{object}{ A fitted SPLS object. }
  \item{coverage}{ Coverage of confidence intervals.
                    \code{coverage} should have a number between 0 and 1.
                    Default is 0.95 (95\eqn{\%} confidence interval). }
  \item{B}{ Number of bootstrap iterations. Default is 1000. }
  \item{plot.it}{ Plot confidence intervals of coefficients? }
  \item{plot.fix}{ If \code{plot.fix="y"}, then plot confidence intervals
                of the predictors for a given response.
                If \code{plot.fix="x"}, then plot confidence intervals
                of a given predictor across all the responses.
                Relevant only when \code{plot.it=TRUE}. }
  \item{plot.var}{ Index vector of responses (if  \code{plot.fix="y"})
                or predictors (if  \code{plot.fix="x"}) to be fixed in \code{plot.fix}.
                The indices of predictors are defined
                among the set of the selected predictors.
                Relevant only when \code{plot.it=TRUE}. }
  \item{K}{ Number of hidden components.
            Default is to use the same \code{K} as in the original SPLS fit. }
  \item{fit}{ PLS algorithm for model fitting. Alternatives are
                    \code{"kernelpls"}, \code{"widekernelpls"},
                    \code{"simpls"}, or \code{"oscorespls"}.
                    Default is to use the same PLS algorithm
                    as in the original SPLS fit.}
}
\value{ Invisibly returns a list with components:
  \item{cibeta}{ A list with as many matrix elements as the number of responses.
                Each matrix element is p by 2, where i-th row of the matrix
                lists the upper and lower bounds of the bootstrapped confidence
                interval of the i-th predictor. }
  \item{betahat}{ Matrix of original coefficients of the SPLS fit. }
  \item{lbmat}{ Matrix of lower bounds of confidence intervals (for internal use). }
  \item{ubmat}{ Matrix of upper bounds of confidence intervals (for internal use). }
}
\references{
Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection",
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25.}
\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ \code{\link{correct.spls}} and \code{\link{spls}}. }
\examples{
data(mice)
# SPLS with eta=0.6 & 1 hidden components
f <- spls( mice$x, mice$y, K=1, eta=0.6 )
# Calculate confidence intervals of coefficients
ci.f <- ci.spls( f, plot.it=TRUE, plot.fix="x", plot.var=20 )
# Bootstrapped confidence intervals
cis <- ci.f$cibeta
cis[[20]]   # equivalent, 'cis$1422478_a_at'
}
\keyword{regression}
\keyword{multivariate}
