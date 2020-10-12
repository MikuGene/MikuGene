\name{MK_reads}
\alias{MK_reads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_reads: read the MK_toMMs() output.
}
\description{
Read the sparse matrix data processed by MK_toMMs().
}
\usage{
MK_reads(path, IDin = NULL, verbose = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
The path MK_toMMs() has output. Default getwd().
}
  \item{IDin}{
The key file name parts that can distinguish file batches uniquely. Default NULL (automatically identify all available data in the target path).
}
  \item{verbose}{
Messages in the process will be displayed. Default TRUE.
}
}
\details{
When reading, please ensure that the folder has not been modified.
}
\value{
A sparse matrix.
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
MK_toMMs()
}
\examples{
MK_reads()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
