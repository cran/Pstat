\name{BootPstk}
\alias{BootPstk}

\title{
Bootstrap method
}
\description{
'BootPstk' performs a bootstrap resampling procedure with all individuals of the selected populations and calculates Pst values of quantitative measures considered. This function provides a confidence interval or the distribution of Pst.
}
\usage{
BootPstk(M, b = 1, csh = 1, boot = 1000, k, Ri = 0, Rp = 0, Pw = 0, pe = 0.95, 
bars = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{b}{
if b=1 all the boot values of Pst are returned, if b=2 these values (ordered) and the confidence interval are returned, if b=3 these values (ordered) and Pst distribution histogram are returned.
}
  \item{csh}{
the value of c/h^2, where c is the assumed additive genetic proportion of differences between populations and where h^2 is (narrow-sense heritability) the  assumed  additive  genetic  proportion of differences between individuals within populations.
}
  \item{boot}{
the number of data frames generated to determine the confidence interval or to construct the distribution (with the bootstrap method).
}
  \item{k}{
the number of the column containing the quantitative measure considered.
}
  \item{Ri}{
a vector containing each number of individual to be deleted.
The vector Ri must contain existent individuals, each of them once.
}
  \item{Rp}{
a vector containing the names of the populations to be deleted.
}
  \item{Pw}{
a vector containing the names of the two populations considered to obtain pairwise Pst.
}
  \item{pe}{
the confidence level of the calculated interval.
}
  \item{bars}{
the maximum number of bars the histogram may have. Note that unfilled bars are taken into account.
}
}

\value{
In any case, the sizes of each population considered.
If a=1 a vector containing the boot values of Pst.
If a=2 an ordered vector containing values of Pst and the confidence interval (also a vector).
If a=3 an ordered vector containing values of Pst and the Pst distribution histogram.
}

\author{
Blondeau Da Silva Stephane - Da Silva Anne.
}


\examples{

data(test)
# BootPstk(test,k=2)
# BootPstk(test,k=7,b=2,csh=0.8,boot=2000,Ri=18)
# BootPstk(test,k=12,b=2,Ri=c(22,27,195),Rp=c("A","B","E"),pe=0.9)
BootPstk(test,k=4,b=3,Ri=c(1:160),bars=50)
# BootPstk(test,k=4,b=3,Ri=c(3,7:17),Pw=c("C","D"),bars=50)

}
