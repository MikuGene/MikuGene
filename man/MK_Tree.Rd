\name{MK_Tree}
\alias{MK_Tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_Tree: Step by step estimation of cell composition.
}
\description{
Estimate the cell components in the bulk-seq expression profile according to the hierarchical relationship between cell types, and output the tree-results.
}
\usage{
MK_Tree(Bulk, Sigl, Cluster, SuClust, SuType = c("Cancer cell", "Immune cell"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bulk}{
Bulk-seq expression matrix.
}
  \item{Sigl}{
ScRNA-seq expression matrix.
}
  \item{Cluster}{
Names of all cells in first tree.
}
  \item{SuClust}{
Names of all cells in second tree.
}
  \item{SuType}{
Choose the second tree sub-types.
}
}
\details{
This function is still being optimized ...
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
Please waiting ...
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
