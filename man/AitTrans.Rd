\name{AitTrans}
\alias{AitTrans}
\title{
Aitchison transformation
}
\description{
'AitTrans' calculates log-ratio transformed measurements to account for individual size-effects.
}
\usage{
AitTrans(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
}

\value{
the data frame with transformed variables. 
}
\references{
Aitchison J., 1986. The Statistical Analysis of Compositional Data. Chapman and Hall, London - New York, XII, 416 pp.

Peres-Neto P.R., Magnan P., 2004. The influence of swimming demand on phenotypic plasticity and morphological integration: a comparison of two polymorphic charr species. Oecologia 140, 36-45.

Shinn C. et al., 2015. Phenotypic variation as an indicator of pesticide stress in gudgeon: accounting for confounding factors in the wild. Science of the Total Environment 538, 733-742.
}
\author{
Blondeau Da Silva Stephane - Da Silva Anne.
}

\examples{
data(test)
AitTrans(test)
}
