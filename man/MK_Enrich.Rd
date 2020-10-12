\name{MK_Enrich}
\alias{MK_Enrich}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_Enrich: A function to simply enrich GO, KEGG, Reactome.
}
\description{
Simply enrich GO, KEGG and Reactome through the clusterProfiler R package.
}
\usage{
MK_Enrich(x, EnID = "temp", CutP = 0.01, Save = T, Wid = 8, Hig = 8.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Genes, such as c('CD19', 'MS4A1', 'CD79A', ...).
}
  \item{EnID}{
Custom naming of this process, such as file name for automatic storage, etc. Default "temp".
}
  \item{CutP}{
P-value cutoff. Default 0.01.
}
  \item{Save}{
Automatically create a EnID folder, and store the resulting files. Default TRUE.
}
  \item{Wid}{
The width of figures when automatically saved. Default 8.
}
  \item{Hig}{
The height of figures when automatically saved. Default 8.3.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The gene enrichment results (GO, KEGG and Reactome).
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
Gene = c("SPP1", "AKR1C1", "AKR1C2", "ADH7", "KLK6", "UGT1A7", "ALDH1A1")

Enrich = MK_Enrich(Gene, EnID = "test")
Enrich
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
