\name{Pst}
\alias{Pst}

\title{
Pst values and Pst confidence intervals 
}
\description{
'Pst' calculates Pst values of the quantitative measures considered and also their confidence intervals, Pst value being an index that measures the level of phenotypic differentiation among populations (Leinonen et al., 2006).
}
\usage{
Pst(Mat, csh = 1, col = 1, a = 1, boot = 1000, Pw = 0, Rp = 0, Ri = 0, pe = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mat}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{csh}{
the value of c/h^2, where c is the assumed additive genetic proportion of differences between populations and where h^2 is (narrow-sense heritability) the  assumed  additive  genetic  proportion of differences between  individuals within populations.
}
  \item{col}{
a vector with the number of the selected columns (i.e. the quantitative measures considered).
}
  \item{a}{
if a=1 the values of Pst are returned, if a=2 their confidence intervals are added. 
}
  \item{boot}{
the number of data frames generated to determine the confidence interval with the bootstrap method.
}
  \item{Pw}{
a vector containing the names of the two populations considered to obtain pairwise Pst.
}
  \item{Rp}{
a vector containing the names of the populations to be deleted.
}
  \item{Ri}{
a vector containing each number of individual to be deleted.
The vector Ri must contain existent individuals, each of them once.
}
  \item{pe}{
the confidence level of the calculated interval.
}
}

\value{
The sizes of each population considered.
Pst values of the selected populations (for quantitative traits considered) and if a=2 their confidence intervals.
}
\references{
Spitze K., 1993. Population structure in Daphnia obtusa: Quantitative genetic and allozymic variation. Genetics 135: 367-374.

Leinonen T., Cano J.M., Makinen H., Merila J., 2006. Contrasting patterns of body shape and neutral genetic divergence in marine and lake populations of threespine sticklebacks. Journal of Evolutionary Biology  19: 1803-1812.

Brommer J.E., 2011. Whither Pst? The approximation of Qst by Pst in evolutionary and conservation biology. Journal Evolution Biology 24: 1160-1168.
}
\author{
Blondeau Da Silva Stephane - Da Silva Anne.
}





\examples{

data(test)
Pst(test)
# Pst(test,csh=0.2,a=2)
# Pst(test,col=c(2,12),a=2,Rp=c("C","D"))
# Pst(test,col=c(5,8:11),a=2,boot=2000,Ri=56,Rp="A",pe=0.9)
# Pst(test,a=2,Ri=c(7,55:59),Pw=c("A","D"))

}

