\name{TracePst}
\alias{TracePst}

\title{
Pst variations in function of c/h^2
}
\description{
'TracePst' plots the curves of the functions that map c/h^2 onto Pst (for chosen quantitative measures). Indeed, Pst depends on the value of c/h^2, where c is the assumed additive genetic proportion of differences between populations and where h^2 is (narrow-sense heritability) the  assumed additive genetic proportion of differences between individuals within populations.
}
\usage{
TracePst(Mat, col = 1, ci = 1, boot = 1000, pe = 0.95, Fst = -1, Pw = 0, Rp = 0, 
Ri = 0, xm = 2, pts = 30)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mat}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{col}{
a vector with the number of the selected columns (i.e. the quantitative measures considered).
}
  \item{ci}{
if ci=1 the confidence interval of Pst is plotted.
}
  \item{boot}{
the number of data frames generated to determine the confidence interval or to construct the dotted lines representing this confidence interval (using the bootstrap method).
}
  \item{pe}{
the confidence level of the calculated interval.
}
  \item{Fst}{
the value of Wright's Fst, if avalaible.
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
  \item{xm}{
the maximum on x-axis (values of c/h^2).
}
  \item{pts}{
number of points used to plot the curves.
}
}

\value{
In any case, the sizes of each population considered.
The expected curves.
}
\references{
Brommer J., 2011. Whither Pst? The approximation of Qst by Pst in evolutionary and conservation biology. Journal of Evolutionary Biology, 24:1160-1168.


Lima M.R. et al., 2012. Genetic and Morphometric Divergence of an Invasive Bird: The Introduced House Sparrow (Passer domesticus) in Brazil. PloS One 7 (12).

On Fst : 
Wright S., 1951. The genetical structure of populations. Annals of Eugenics 15, 323-354.
}
\author{
Blondeau Da Silva Stephane - Da Silva Anne.
}

\note{
The time required to construct the dotted lines associated with the confidence intervals might be fairly long depending on the user choices.
}


\examples{

data(test)
# TracePst(test)
# TracePst(test,boot=2000,col=7,Ri=18,pe=0.9,pts=40,xm=4)
TracePst(test,col=c(7,9:11),Fst=0.3,Ri=c(22,27,195),Rp=c("A","C","E"),ci=0)
# TracePst(test,Ri=c(3,7:17),Pw=c("C","D"))

}
