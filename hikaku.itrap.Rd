\name{hikaku.itrap}
\alias{hikaku.itrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process label-free comparative proteomics data
}
\description{
This function returns comparative log2 ratio analysis and mean peak area heatmap comparison from label-free proteomics or an exploratory graphical sample data analysis.
}
\usage{
hikaku.itrap("data.csv", samplen, nrep, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.csv}{A tab-delimited .csv file in Protein Discoverer (ThermoLife) output format.}
  \item{samplen}{A natural number >1 and <16 for samples/conditions}
  \item{nrep}{A natural number >0 for replicates of each sample/condition}
  \item{expl}{TRUE for exploratory analysis and FALSE for comparative label-free proteomics analysis. Default: expl=FALSE}
  \item{pep}{A natural number for a threshold of the minimal number of peptides detected. Default: pep=0}
  \item{psm}{A natural number for a threshold of minimal PSM value. Default:psm=0}
  \item{heat.qt}{A rational number (0<x<1) which means the fraction of upper quantile to not appear in the heatmap. Default:heat.qt=0.1}
  \item{pnames}{TRUE for protein id as names and FALSE for Accession number. Default:pnames=TRUE}
  \item{housek}{Accession number with "" for the housekeeping protein selected for internal sample normalization}
}
\details{

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
Tiago A. de Souza (tiagoantonio@gmail.com)
Fernando G. de Almeida (fernandoalmeida@usp.br)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
