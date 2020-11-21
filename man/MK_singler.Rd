\name{MK_singler}
\alias{MK_singler}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_singler: Identification of the cell types in scRNA-seq data.
}
\description{
Supervised identification of the cell types in single-cell sequencing data by SingleR package.
}
\usage{
MK_singler(x, Ref = "HPCA", mode = "main", cluster = NULL, Cells = 10, name = NULL, Save = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A two-dimensional matrix, including sparse matrix. (Column: cells; Row: features)
}
  \item{Ref}{
Reference data set from celldex package, including: HPCA; BPED and DICE. Default 'HPCA'.
}
  \item{mode}{
Accuracy mode of cell type recognition, including: main and fine. Default 'main'.
}
  \item{cluster}{
Cluster identities of each cell in x.
}
  \item{Cells}{
The number of cells (k) used to run separately. Default 10 (10,000 cells).
}
  \item{name}{
Custom naming of this process. Default 'temp'.
}
  \item{Save}{
The process will automatically create a backup folder under the working path, and store the resulting object as a name.csv file. Default TRUE.
}
}
\details{
Please:
  citation("SingleR")
  citation("celldex")
}
\value{
A vector of cell-type identities.
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
CIden = MK_singler(x, Ref = "HPCA", mode = "main")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
