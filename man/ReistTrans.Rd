\name{ReistTrans}
\alias{ReistTrans}

\title{
Reist standardization
}
\description{
'ReistTrans' calculates residuals (size adjusted measurements) from Reist tranformations to eliminate any variation resulting from allometric growth. There is a single regressor (one of the quantitative traits).
}
\usage{
ReistTrans(M, r)
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

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

data(test)
names(test)[9]
ReistTrans(test,r=9)
}
