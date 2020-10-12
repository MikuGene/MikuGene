\name{MK_read10X}
\alias{MK_read10X}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_read10X: Read 10x genomic scRNA-seq data.
}
\description{
Read 10x genomic single-cell sequencing data in batches by entering the key file name part.
}
\usage{
MK_read10X(MKdir = getwd(), IDin = NULL, Barfile = "barcode", Genefile = "gene", Exprfile = "matrix", View = T, Save = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MKdir}{
The file path where the 10X-files are located. Default getwd().
}
  \item{IDin}{
The key file name parts that can distinguish file batches uniquely. Default NULL (automatically identify all available data in the target path).
}
  \item{Barfile}{
The key file name part of the Barcode file. Default "barcode".
}
  \item{Genefile}{
The key file name part of the Feature file. Default "gene".
}
  \item{Exprfile}{
The key file name part of the Matrix file. Default "matrix".
}
  \item{View}{
Messages in the process will be displayed. Default TRUE.
}
  \item{Save}{
Save the meta-feature (such as antibody) data under the working path.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of sparse matrices (Column: cells; Row: features).
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
MK_read10X()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
