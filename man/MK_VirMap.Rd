\name{MK_VirMap}
\alias{MK_VirMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_VirMap: Check virus in high-throughput sequencing.
}
\description{
Align the reads of high-throughput sequencing to the viral reference genome.
}
\usage{
MK_VirMap(path_r1, path_r2, name = "temp", maxMiss = 3, GTF = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path_r1}{
Fastq R1 (or .gz)
}
  \item{path_r2}{
Fastq R2 (or .gz)
}
  \item{name}{
Custom naming of this process. Default "temp".
}
  \item{maxMiss}{
Max miss-match in Map. Default 3.
}
  \item{GTF}{
Use GTF annotation files. Default TRUE.
}
}
\details{
Please:
  citation("Rsubread")
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
MK_BuildVirusRef()
MK_VirMap(Fq1, Fq2, name = "FQ")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
