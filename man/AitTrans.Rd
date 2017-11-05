\name{AitTrans}
\alias{AitTrans}
\title{
Aitchison transformation
}
\description{
'AitTrans' calculates log-ratio transformed measurements to account for individual size-effects.
}
\usage{
AitTrans(data, Rp = 0, Ri = 0)
}

\arguments{
  \item{data}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{Rp}{
a vector containing the names of the populations to be deleted.
}
  \item{Ri}{
a vector containing each number of individual to be deleted. The vector Ri must contain existent individuals, each of them once.
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

## The function is currently defined as
function (data, Rp = 0, Ri = 0) 
{
    nb.var = dim(data)[2] - 1
    dat.rem.ind.pop <- function(data, ind = 0, pop = 0) {
        data = as.data.frame(data)
        dat.rem.ind <- function(dat, ind) {
            nb.rem.ind = length(ind)
            nb.ind = dim(dat)[1]
            for (i in 1:nb.rem.ind) dat = dat[row.names(dat)[1:(nb.ind - 
                i + 1)] != ind[i], ]
            return(dat)
        }
        dat.rem.pop <- function(dat, pop) {
            nb.rem.pop = length(pop)
            for (i in 1:nb.rem.pop) dat = dat[dat[, 1] != pop[i], 
                ]
            return(dat)
        }
        if (ind[1] != 0) 
            data = dat.rem.ind(data, ind)
        if (pop[1] != 0) 
            data = dat.rem.pop(data, pop)
        return(data)
    }
    nonNa.row <- function(data, row) {
        nb.na = 0
        for (i in 1:nb.var) if (is.na(data[row, i + 1])) 
            nb.na = nb.na + 1
        return(nb.var - nb.na)
    }
    data = dat.rem.ind.pop(data, ind = Ri, pop = Rp)
    nb.ind = dim(data)[1]
    cor.val = rep(0, nb.ind)
    for (i in 1:nb.ind) {
        cor.val[i] = 1/(nonNa.row(data, i)) * sum(log(data[i, 
            -1], base = 10), na.rm = TRUE)
        for (j in 1:nb.var) data[i, (j + 1)] = log(data[i, (j + 
            1)], base = 10) - cor.val[i]
    }
    return(data)
  }
}

