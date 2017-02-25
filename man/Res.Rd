\name{Res}
\alias{Res}
\title{
Residuals from a linear regression
}
\description{
'Res' calculates residuals from simple linear regressions (in particular to eliminate any variation resulting from allometric growth). These regression adjustments assume the existence of linear relationships between the dependent variables and the regressor (one of the column of the data frame).
}
\usage{
Res(M, r)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{r}{
the number of the column containing the quantitative trait used as a regressor. 
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

data(test)
names(test)[9]
Res(test,r=9)

}
