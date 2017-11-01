\name{ReistTrans}
\alias{ReistTrans}

\title{
Reist standardization
}
\description{
'ReistTrans' calculates residuals (size adjusted measurements) from Reist tranformations to eliminate any variation resulting from allometric growth. There is a single regressor (one of the quantitative traits).
}

\usage{
ReistTrans(data, reg)
}

\arguments{
  \item{data}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{reg}{
the name (or the rank) of the variable chosen as the explanatory variable.
}
}

\value{
the data frame of adjusted variables, the column containing the quantitative trait used as a regressor being deleted.
}

\references{
Reist J.D., 1985. An empirical evaluation of several univariate methods that adjust for size variation in morphometric data. Canadian Journal Zoology 63, 1429-1439.

Kaeuffer R. et al., 2012. Parallel and nonparallel aspects of ecological, phenotypic, and genetic divergence across replicate population pairs of lake and stream stickleback. Evolution 66(2), 402-418.

He Y. et al., 2013. Morphological Variation Among Wild Populations of Chinese Rare Minnow (Gobiocypris rarus): Deciphering the Role of Evolutionary Processes. Zoological Science 30, 475-483.
}

\author{
Blondeau Da Silva Stephane - Da Silva Anne. 
}

\note{
dispensable quantitative measures can easily be deleted in the main functions of R.
}


\examples{
data(test)
names(test)[9]
ReistTrans(test,reg=9)

## The function is currently defined as
function (data, reg) 
{
    Reitra.va <- function(dat, clm, re) {
        dat = dat[is.finite(dat[, re]), ]
        log.dat = dat
        m = mean(dat[is.finite(dat[, clm]), re])
        Q = dim(dat)[2]
        log.dat[, clm] = log(dat[, clm], base = 10)
        log.dat[, re] = log(dat[, re], base = 10)
        m.clm = mean(log.dat[is.finite(log.dat[, clm]), clm], 
            na.rm = TRUE)
        m.reg = mean(log.dat[is.finite(log.dat[, clm]), re], 
            na.rm = TRUE)
        a = sum((log.dat[is.finite(log.dat[, clm]), re] - m.reg) * 
            log.dat[is.finite(log.dat[, clm]), clm])/sum((log.dat[is.finite(log.dat[, 
            clm]), re] - m.reg) * (log.dat[is.finite(log.dat[, 
            clm]), re] - m.reg))
        dat[, clm] = log.dat[, clm] - a * (log.dat[, re] - log(m, 
            base = 10))
        return(dat)
    }
    for (i in 2:dim(data)[2]) {
        if (names(data)[i] == reg) 
            reg = i - 1
    }
    if (is.numeric(reg) == FALSE) 
        return("reg value does not exist!")
    Q = dim(data)[2]
    if (reg == 1) 
        for (i in 3:Q) data = Reitra.va(data, clm = i, re = 2)
    else {
        for (i in 2:reg) data = Reitra.va(data, clm = i, re = reg + 
            1)
        if (reg < (Q - 1)) 
            for (j in (reg + 2):Q) data = Reitra.va(data, clm = j, 
                re = reg + 1)
    }
    return(data[-(reg + 1)])
  }
}

