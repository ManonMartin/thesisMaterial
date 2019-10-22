\name{prostate}
\docType{data}
\alias{prostate}
\title{Prostate Tumor Gene Expression Dataset}
\description{
  This is the Prostate Tumor Gene Expression dataset used in Chung and Keles (2010).
}
\usage{ data(prostate) }
\format{
  A list with two components:
  \describe{
  \item{x}{ Gene expression data. A matrix with 102 rows and 6033 columns.}
  \item{y}{ Class index. A vector with 102 elements.}
  }
}
\details{
  The prostate dataset consists of 52 prostate tumor and 50 normal samples.
  Normal and tumor classes are coded in 0 and 1, respectively, in \code{y} vector.
  Matrix \code{x} is gene expression data and
  arrays were normalized, log transformed, and standardized 
  to zero mean and unit variance across genes as described
  in Dettling (2004) and Dettling and Beuhlmann (2002).
  See Chung and Keles (2010) for more details. 
}
\source{
Singh D, Febbo P, Ross K, Jackson D, Manola J, Ladd C, Tamayo P, Renshaw A,
DAmico A, Richie J, Lander E, Loda M, Kantoff P, Golub T, and Sellers W (2002),
"Gene expression correlates of clinical prostate cancer behavior", \emph{Cancer Cell}, 
Vol. 1, pp. 203--209.
}
\references{
Chung D and Keles S (2010), 
"Sparse partial least squares classification for high dimensional data",
\emph{Statistical Applications in Genetics and Molecular Biology}, Vol. 9, Article 17.
  
Dettling M (2004), "BagBoosting for tumor classification with gene expression data",
\emph{Bioinformatics}, Vol. 20, pp. 3583--3593.

Dettling M and Beuhlmann P (2002), "Supervised clustering of genes", \emph{Genome Biology},
Vol. 3, pp. research0069.1--0069.15.
}
\examples{
data(prostate)
prostate$x[1:5,1:5]
prostate$y
}
\keyword{datasets}
