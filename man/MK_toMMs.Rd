\name{MK_toMMs}
\alias{MK_toMMs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_toMMs: Save and process scRNA-seq matrix.
}
\description{
Save and simply process the single-cell sequencing expression matrix to distributed sparse matrices.
}
\usage{
MK_toMMs(x, name = "temp", Cells = 10, verbose = T, HK_bm = F, Mito_rm = T, AC_rm = T, RP_rm = T, RPLS_rm = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A two-dimensional matrix, including sparse matrices.
}
  \item{name}{
Custom naming of this process, which creating a 'name' folder in working path. Default "temp".
}
  \item{Cells}{
The number of cells (k) used to store the matrix separately. Default 10 (10,000 cells).
}
  \item{verbose}{
Messages in the process will be displayed. Default TRUE.
}
  \item{HK_bm}{
Perform batch correction through housekeeping genes. Default FALSE.
}
  \item{Mito_rm}{
Delete mitochondrial genes. Default TRUE.
}
  \item{AC_rm}{
Delete lncRNA genes. Default TRUE.
}
  \item{RP_rm}{
Delete pseudogenes genes. Default TRUE.
}
  \item{RPLS_rm}{
Delete ribosome-related genes. Default TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
MK_toMMs(SC, name = "ScRNA")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
