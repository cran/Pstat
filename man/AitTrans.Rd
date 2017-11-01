\name{AitTrans}
\alias{AitTrans}
\title{
Aitchison transformation
}
\description{
'AitTrans' calculates log-ratio transformed measurements to account for individual size-effects.
}
\usage{
AitTrans(data)
}

\arguments{
  \item{data}{
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
# data(test)
# AitTrans(test)

## The function is currently defined as
function (data) 
{
    nonNa.clm <- function(data, clm) {
        n = dim(data)[1]
        c = 0
        for (i in 1:n) if (is.na(data[i, clm])) 
            c = c + 1
        return(n - c)
    }
    nb.ind = dim(data)[1]
    nb.va = dim(data)[2] - 1
    v = rep(0, nb.ind)
    for (i in 1:nb.ind) {
        v[i] = 1/(nonNa.clm(t(data[i, -1]), 1)) * sum(log(data[i, 
            -1], base = 10), na.rm = TRUE)
        for (j in 1:nb.va) data[i, (j + 1)] = log(data[i, (j + 
            1)], base = 10) - v[i]
    }
    return(data)
  }
}

