\name{lymphoma}
\docType{data}
\alias{lymphoma}
\title{Lymphoma Gene Expression Dataset}
\description{
  This is the Lymphoma Gene Expression dataset used in Chung and Keles (2010).
}
\usage{ data(lymphoma) }
\format{
  A list with two components:
  \describe{
  \item{x}{ Gene expression data. A matrix with 62 rows and 4026 columns.}
  \item{y}{ Class index. A vector with 62 elements.}
  }
}
\details{
  The lymphoma dataset consists of 42 samples of diffuse large B-cell lymphoma (DLBCL),
  9 samples of follicular lymphoma (FL),
  and 11 samples of chronic lymphocytic leukemia (CLL).  
  DBLCL, FL, and CLL classes are coded in 0, 1, and 2, respectively, in \code{y} vector.
  Matrix \code{x} is gene expression data and
  arrays were normalized, imputed, log transformed, and standardized 
  to zero mean and unit variance across genes as described
  in Dettling (2004) and Dettling and Beuhlmann (2002).
  See Chung and Keles (2010) for more details. 
}
\source{
Alizadeh A, Eisen MB, Davis RE, Ma C, Lossos IS, Rosenwald A, Boldrick JC,
Sabet H, Tran T, Yu X, Powell JI, Yang L, Marti GE, Moore T, Hudson J Jr, 
Lu L, Lewis DB, Tibshirani R, Sherlock G, Chan WC, Greiner TC, 
Weisenburger DD, Armitage JO, Warnke R, Levy R, Wilson W, Grever MR, Byrd JC,
Botstein D, Brown PO, and Staudt LM (2000), "Distinct types of diffuse large
B-cell lymphoma identified by gene expression profiling", \emph{Nature}, Vol. 403, pp. 503--511.
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
data(lymphoma)
lymphoma$x[1:5,1:5]
lymphoma$y
}
\keyword{datasets}
