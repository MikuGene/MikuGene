\name{MK_WG_Tom}
\alias{MK_WG_Tom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_WG_Tom: One-step WGCNA analysis.
}
\description{
One-step WGCNA analysis based on WGCNA R package.
}
\usage{
MK_WG_Tom(x, name = "temp", nGene = 10000, Save = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A two-dimensional matrix. (Column: samples; Row: features)
}
  \item{name}{
Custom naming of this process. Default "temp".
}
  \item{nGene}{
The number of genes will be used for co-expression network analysis.
}
  \item{Save}{
Create a backup folder under the working path, and store the resulting object as a name_WGtom_backup.rds file. Meanwile the co-expression network diagram will be drawn and saved.
}
}
\details{
Please:
  citation("WGCNA")
}
\value{
A list (WG_Tom) of five objects: x, MEs (Column: each ME; Row: each sample), Colors (color corresponding to each gene), Power (the power of WGCNA soft-threshold), TOM (the calculated topological matrix).
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
