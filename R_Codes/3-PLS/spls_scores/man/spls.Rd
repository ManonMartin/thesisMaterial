\name{spls}
\alias{spls}
\title{Fit SPLS regression models}
\description{
 Fit a SPLS regression model.
}
\usage{
spls( x, y, K, eta, kappa=0.5, select="pls2", fit="simpls",
    scale.x=TRUE, scale.y=FALSE, eps=1e-4, maxstep=100, trace=FALSE)
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector or matrix of responses. }
  \item{K}{ Number of hidden components. }
  \item{eta}{ Thresholding parameter. \code{eta} should be between 0 and 1. }
  \item{kappa}{ Parameter to control the effect of
     the concavity of the objective function
     and the closeness of original and surrogate direction vectors.
     \code{kappa} is relevant only when responses are multivariate.
     \code{kappa} should be between 0 and 0.5. Default is 0.5. }
  \item{select}{ PLS algorithm for variable selection.
                    Alternatives are \code{"pls2"} or \code{"simpls"}.
                    Default is \code{"pls2"}.}
  \item{fit}{ PLS algorithm for model fitting. Alternatives are
                    \code{"kernelpls"}, \code{"widekernelpls"},
                    \code{"simpls"}, or \code{"oscorespls"}.
                    Default is \code{"simpls"}.}
  \item{scale.x}{ Scale predictors by dividing each predictor variable
                by its sample standard deviation? }
  \item{scale.y}{ Scale responses by dividing each response variable
                by its sample standard deviation? }
  \item{eps}{ An effective zero. Default is 1e-4. }
  \item{maxstep}{ Maximum number of iterations when fitting direction vectors.
                    Default is 100. }
  \item{trace}{ Print out the progress of variable selection? }
}
\details{
The SPLS method is described in detail in Chun and Keles (2010).
SPLS directly imposes sparsity on the dimension reduction step of PLS
in order to achieve accurate prediction and variable selection simultaneously.
The option \code{select} refers to the PLS algorithm for variable selection.
The option \code{fit} refers to the PLS algorithm for model fitting
and \code{spls} utilizes algorithms offered by the \pkg{pls} package for this purpose.
See help files of the function \code{plsr} in the \pkg{pls} package for more details.
The user should install the \pkg{pls} package before using \pkg{spls} functions.
The choices for \code{select} and \code{fit} are independent.
}
\value{A \pkg{spls} object is returned.
print, plot, predict, coef, ci.spls, coefplot.spls methods use this object.}
\references{
    Chun H and Keles S (2010), "Sparse partial least squares
    for simultaneous dimension reduction and variable selection",
    \emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25. 
}
\author{ 
    Dongjun Chung, Hyonho Chun, and Sunduz Keles. 
}
\seealso{ \code{\link{print.spls}}, \code{\link{plot.spls}}, \code{\link{predict.spls}},
\code{\link{coef.spls}}, \code{\link{ci.spls}}, and \code{\link{coefplot.spls}}. }
\examples{
    data(yeast)
    # SPLS with eta=0.7 & 8 hidden components
    (f <- spls(yeast$x, yeast$y, K=8, eta=0.7))

    # Print out coefficients
    coef.f <- coef(f)
    coef.f[,1]

    # Coefficient path plot
    plot(f, yvar=1)
    dev.new()

    # Coefficient plot of selected variables
    coefplot.spls(f, xvar=c(1:4))
}
\keyword{regression}
\keyword{multivariate}
