\name{Res}
\alias{Res}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residuals from a linear regression
}
\description{
'Res' calculates residuals from simple linear regressions (in particular to eliminate any variation resulting from allometric growth). These regression adjustments assume the existence of linear relationships between the dependent variables and the regressor (one of the column of the data frame).
}
\usage{
Res(data, reg)
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

\author{
Blondeau Da Silva Stephane - Da Silva Anne.
}
\note{
dispensable quantitative measures can easily be deleted in the main functions of R.
}

\examples{

# data(test)
# names(test)[9]
# Res(test,r=9)

## The function is currently defined as
function (data, reg) 
{
    Res.va <- function(dat, clm, re) {
        dat = dat[is.finite(dat[, re]), ]
        m.clm = mean(dat[is.finite(dat[, clm]), clm], na.rm = TRUE)
        m.reg = mean(dat[is.finite(dat[, clm]), re], na.rm = TRUE)
        a = sum((dat[is.finite(dat[, clm]), re] - m.reg) * dat[is.finite(dat[, 
            clm]), clm])/sum((dat[is.finite(dat[, clm]), re] - 
            m.reg) * (dat[is.finite(dat[, clm]), re] - m.reg))
        b = m.clm - a * m.reg
        dat[, clm] = dat[, clm] - b - a * dat[, re]
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
        for (i in 3:Q) data = Res.va(data, clm = i, re = 2)
    else {
        for (i in 2:reg) data = Res.va(data, clm = i, re = reg + 
            1)
        if (reg < (Q - 1)) 
            for (j in (reg + 2):Q) data = Res.va(data, clm = j, 
                re = reg + 1)
    }
    return(data[-(reg + 1)])
  }
}

