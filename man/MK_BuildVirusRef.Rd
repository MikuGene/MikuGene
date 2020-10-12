\name{MK_BuildVirusRef}
\alias{MK_BuildVirusRef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_BuildVirusRef
}
\description{
Use the data of virusite.org to make reference files for virus.
}
\usage{
MK_BuildVirusRef(version = "2020.3")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{version}{
the version of the data.
}
  \item{OutVs}{
The SeqID of the virus genome not included in alignment. If none, input NULL . Default "default".
}
  \item{verbose}{
Messages in the process will be displayed. Default TRUE.
}
}
\details{
Please:
  citation("Biostrings")
  citation: http://www.virusite.org/
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
MK_BuildVirusRef()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
