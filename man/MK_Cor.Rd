\name{MK_Cor}
\alias{MK_Cor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_Cor: Multi-item correlation analysis
}
\description{
Perform correlation analysis on multiple items cyclically.
}
\usage{
MK_Cor(x, y, method = "all", p_cut = 0.01, adj = T, name = NULL, Save = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Source group: item or multi-item.
}
  \item{y}{
Target group: item or multi-item.
}
  \item{method}{
Correlation analysis method to be used, including pearson and spearman (and all). Default 'all'.
}
  \item{p_cut}{
P-value cutoff considered significant. Default '0.01'
}
  \item{adj}{
The 'BH' method will be used to adjust the p-value. Default TRUE.
}
  \item{name}{
Custom naming of this process. Default 'temp'.
}
  \item{Save}{
Automatically create a backup folder under the working path, and store a name_cor.csv file. Default TRUE.
}
}
\details{
}
\value{
A data-frame contains correlation analysis results.
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
Please waiting ...
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line