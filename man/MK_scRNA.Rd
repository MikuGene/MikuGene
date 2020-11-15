\name{MK_scRNA}
\alias{MK_scRNA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MK_scRNA: One-step scRNA-seq data analysis.
}
\description{
One-step single-cell RNA-seq data analysis based on Seurat v3 and harmony R package.
}
\usage{
MK_scRNA(x, name = NULL, Reso = 0.8, nGene = c(200, Inf), nVar = 2.5, Dim = 2, SCT = F, BatchRemove = F, Umap = F, Plot = T, Norm = T, save = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A two-dimensional matrix, including sparse matrix. (Column: cells; Row: features)
}
  \item{name}{
Custom naming of this process, such as file name for backup. Default "temp".
}
  \item{Reso}{
The resolution of cell clustering (Seurat::FindCluster() can be used to change the resolution of the output result). Default 0.8.
}
  \item{nGene}{
Screen cells by setting an interval for the number of genes contained in each cell in standard procedure. Default c(200, Inf).
}
  \item{nVar}{
Filter the number of highly variable genes (k). Default 2.5 (2500 variable genes).
}
  \item{Dim}{
The dimension of unsupervised clustering. Default 2.
}
  \item{SCT}{
The 'sctransform' process recommended by Seurat v3 will be used. Default FALSE.
}
  \item{BatchRemove}{
Harmony will be applied to remove batch effects according "orig.ident". Default FALSE.
}
  \item{Umap}{
"RunUMAP" function in Seurat will be applied. Default FALSE.
}
  \item{Plot}{
%%     ~~Describe \code{Plot} here~~
}
  \item{Norm}{
"LogNormalize" will be applied in standard procedure. Default TRUE.
}
  \item{save}{
Automatically create a backup folder under the working path, and store the resulting object as a name-time_backup.rds file. Default TRUE.
}
}
\details{
Please:
  citation("Seurat")
  citation("harmony")
}
\value{
a Seurat object.
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
SC = MK_reads("ScRNA")
SC = MK_scRNA(SC)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
