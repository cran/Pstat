\name{Pstat-package}
\alias{Pstat-package}
\alias{Pstat}
\docType{package}
\title{
Assessing Pst Statistics
}
\description{
This package aims at calculating Pst values to assess differentiation among populations from a set of quantitative traits (with the function 'Pst'). The bootstrap method provides confidence intervals and distribution histograms of Pst (with the function 'BootPst'). Variations of Pst in function of the parameter c/h^2 are studied (with the function 'TracePst') as well.
Moreover the package allows users to transform theirs datas in three different ways in particular to eliminate any variation resulting from allometric growth (calculation of residuals from linear regressions -with the function 'Res'-, Reist standardizations -with the function 'ReistTrans'- or Aitchison transformation -with the function 'AitTrans'-).
}

\details{

The DESCRIPTION file:
This package aims at calculating Pst values to assess differentiation among populations from a set of quantitative traits (with the function 'Pst'). The bootstrap method provides confidence intervals and distribution histograms of Pst (with the function 'BootPst'). Variations of Pst in function of the parameter c/h^2 are studied (with the function 'TracePst') as well.
Moreover the package allows users to transform theirs datas in three different ways in particular to eliminate any variation resulting from allometric growth (calculation of residuals from linear regressions -with the function 'Res'-, Reist standardizations -with the function 'ReistTrans'- or Aitchison transformation -with the function 'AitTrans'-).
}

\author{
Blondeau Da Silva Stephane [aut, cre] - Da Silva Anne [aut].

Maintainer: Blondeau Da Silva Stephane <blondeaudasilva@xlim.fr>
}

\references{
Spitze K., 1993. Population structure in Daphnia obtusa: Quantitative genetic and allozymic variation. Genetics 135: 367-374.

Merila J., Crnokrak P., 2001. Comparison of genetic differentiation at marker loci and quantitative traits. Journal Evolution Biology 14: 892-903.

Leinonen T., O'Hara R.B., Cano J.M., Merila J., 2008. Comparative studies of quantitative trait and neutral marker divergence: a meta-analysis. Journal Evolution Biology 21: 1-17.

Brommer J.E., 2011. Whither Pst? The approximation of Qst by Pst in evolutionary and conservation biology. Journal Evolution Biology 24: 1160-1168.
}

\examples{
# data(test)
# test1=AitTrans(test)
# Pst(test1,csh=0.2,ci=1)
# test2=ReistTrans(test,reg="Body_length")
# BootPst(test2,va="QM4",opt="hist",Ri=c(3,7:17),Pw=c("C","D"),bars=50)
# TracePst(test2,va=c(7,9:12),Fst=0.3,Ri=c(22,27,195),Rp=c("A","C","E"))
}
