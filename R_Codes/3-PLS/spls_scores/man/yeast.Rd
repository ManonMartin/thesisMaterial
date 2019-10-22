\name{yeast}
\docType{data}
\alias{yeast}
\title{Yeast Cell Cycle Dataset}
\description{
  This is the Yeast Cell Cycle dataset used in Chun and Keles (2010).
}
\usage{ data(yeast) }
\format{
  A list with two components:
  \describe{
  \item{x}{ ChIP-chip data. A matrix with 542 rows and 106 columns.}
  \item{y}{ Cell cycle gene expression data.
                A matrix with 542 rows and 18 columns.}
  }
}
\details{
  Matrix \code{y} is cell cycle gene expression data (Spellman et al., 1998)
  of 542 genes from an \eqn{\alpha} factor based experiment.
  Each column corresponds to mRNA levels
  measured at every 7 minutes during 119 minutes (a total of 18 measurements).
  Matrix \code{x} is the chromatin immunoprecipitation on chip (ChIP-chip) data of
  Lee et al. (2002) and it contains the binding information for 106
  transcription factors. See Chun and Keles (2010) for more details. 
}
\source{
Lee TI, Rinaldi NJ, Robert F, Odom DT, Bar-Joseph Z, Gerber GK, Hannett NM,
Harbison CT, Thomson CM, Simon I, Zeitlinger J, Jennings EG, Murray HL,
Gordon DB, Ren B, Wyrick JJ, Tagne JB, Volkert TL, Fraenkel E, Gifford DK,
and Young RA (2002), "Transcriptional regulatory networks in \emph{Saccharomyces cerevisiae}",
\emph{Science}, Vol. 298, pp. 799--804.

Spellman PT, Sherlock G, Zhang MQ, Iyer VR, Anders K, Eisen MB, Brown PO,
Botstein D, and Futcher B (1998), "Comprehensive identification of cell cycle-regulated genes of
the yeast \emph{Saccharomyces cerevisiae} by microarray hydrization",
\emph{Molecular Biology of the Cell}, Vol. 9, pp. 3273--3279.
}
\references{
Chun H and Keles S (2010), "Sparse partial least squares
for simultaneous dimension reduction and variable selection",
\emph{Journal of the Royal Statistical Society - Series B}, Vol. 72, pp. 3--25. }
\examples{
data(yeast)
yeast$x[1:5,1:5]
yeast$y[1:5,1:5]
}
\keyword{datasets}
