\name{cv.spls}
\alias{cv.spls}
\title{Compute and plot cross-validated mean squared prediction error for SPLS regression}
\description{
 Draw heatmap of
 v-fold cross-validated mean squared prediction error
 and return optimal eta (thresholding parameter)
 and K (number of hidden components).
}
\usage{
cv.spls( x, y, fold=10, K, eta, kappa=0.5,
        select="pls2", fit="simpls",
        scale.x=TRUE, scale.y=FALSE, plot.it=TRUE )
}
\arguments{
  \item{x}{ Matrix of predictors. }
  \item{y}{ Vector or matrix of responses. }
  \item{fold}{ Number of cross-validation folds. Default is 10-folds. }
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
  \item{plot.it}{ Draw heatmap of cross-validated mean squared prediction error? }
}
\value{
Invisibly returns a list with components:
  \item{mspemat}{ Matrix of cross-validated mean squared prediction error.
                Rows correspond to \code{eta} and
                columns correspond to the number of components (\code{K}). }
  \item{eta.opt}{ Optimal \code{eta}. }
  \item{K.opt}{ Optimal \code{K}. }
}
\references{
    Chun H and Keles S (2010), "Sparse partial least squares
    for simultaneous dimension reduction and variable selection",
    \emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25.
}

\author{ Dongjun Chung, Hyonho Chun, and Sunduz Keles. }
\seealso{ \code{\link{print.spls}}, \code{\link{plot.spls}}, \code{\link{predict.spls}},
and \code{\link{coef.spls}}. }
\examples{
    data(yeast)
    set.seed(1)

    # MSPE plot. eta is searched between 0.1 and 0.9 and
    # number of hidden components is searched between 1 and 10

    \dontrun{
        cv <- cv.spls(yeast$x, yeast$y, K = c(1:10), eta = seq(0.1,0.9,0.1))
        
        # Optimal eta and K
        cv$eta.opt
        cv$K.opt
        (spls(yeast$x, yeast$y, eta=cv$eta.opt, K=cv$K.opt))
    }
}
\keyword{regression}
\keyword{multivariate}
