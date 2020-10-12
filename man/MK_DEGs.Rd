\name{MK_DEGs}
\alias{MK_DEGs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_DEGs: One-step identification of DEGs.
}
\description{
One-step identification of differentially expressed genes by limma R package.
}
\usage{
MK_DEGs <- function(x, y, filt = T, log2FC = 2, padj = 0.01, pval = 0.01, save = T, Order = T, name = "temp")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Experimental (positive) group matrix.
}
  \item{y}{
Control (negative) group matrix.
}
  \item{filt}{
Filter the results by custom thresholds (log2FC and padj). Default TRUE.
}
  \item{log2FC}{
Custom threshold of log2-foldchange. Default 2.
}
  \item{padj}{
Custom threshold of ajust p-value (BH method). Default 0.01.
}
  \item{Save}{
Automatically save the resulting files. Default TRUE.
}
  \item{Order}{
Sort the analysis results in descending order by log2FC. Default TRUE.
}
  \item{name}{
Custom naming of this process. Default "temp".
}
}
\details{
Please:
  citation("limma")
}
\value{
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Lianhao Song (CodeNight)
}
\note{
---- MikuGene Bioinformatics Ecological Community.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
MK_DEGs(Pos-matrix, Neg-matrix, name = "DEGs")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
