\name{Pst}
\alias{Pst}

\title{
Pst values and Pst confidence intervals 
}
\description{
'Pst' calculates Pst values of the quantitative measures considered and also their confidence intervals.
}
\usage{
Pst(data, ci = 0, csh = 1, va = 0, boot = 1000, Pw = 0, Rp = 0, Ri = 0, pe = 0.95)
}

\arguments{
  \item{data}{
a dataframe with as many rows as individuals. The first column contains the name of the population to which the individual belongs, the others contain quantitative variables.
}
  \item{ci}{
if ci=1 the confidence intervals are added to Pst values.
}
  \item{csh}{
the value of c/h^2, where c is the assumed additive genetic proportion of differences between populations and where h^2 is (narrow-sense heritability) the  assumed  additive  genetic  proportion of differences between  individuals within populations.
}
  \item{va}{
a vector containing the selected variables names or numbers (i.e. those of the quantitative measures considered). If va=0 all the variables are selected. 
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
a vector containing each number of individual to be deleted. The vector Ri must contain existent individuals, each of them once.
}
  \item{pe}{
the confidence level of the calculated interval.
}
}

\value{
The sizes of each population considered.
Pst values of the selected populations (for quantitative traits considered) and if ci=1 their confidence intervals.
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
# Pst(test,csh=0.2,ci=1)
# Pst(test,va="QM2",ci=1,Rp=c("C","D"))
# Pst(test,va=c(5,8:11),ci=1,boot=2000,Ri=56,Rp="A",pe=0.9)
# Pst(test,ci=1,Ri=c(7,55:59),Pw=c("A","D"))

## The function is currently defined as
function (data, ci = 0, csh = 1, va = 0, boot = 1000, Pw = 0, 
    Rp = 0, Ri = 0, pe = 0.95) 
{
    nonNa.clm <- function(data, clm) {
        n = dim(data)[1]
        c = 0
        for (i in 1:n) if (is.na(data[i, clm])) 
            c = c + 1
        return(n - c)
    }
    Prep <- function(data) {
        c = dim(data)[2] - 1
        data = as.data.frame(data)
        data[, 1] = as.character(data[, 1])
        for (i in 1:c) {
            if (is.numeric(data[, i + 1]) == FALSE) 
                data[, i + 1] = as.numeric(as.character(data[, 
                  i + 1]))
        }
        dat.sta <- function(dat) {
            p = dim(dat)[2]
            sd = rep(0, p - 1)
            m = rep(0, p - 1)
            for (i in 1:(p - 1)) {
                N = nonNa.clm(dat, i + 1)
                sd[i] = sqrt((N - 1)/N) * sd(dat[, i + 1], na.rm = TRUE)
                m[i] = mean(dat[, i + 1], na.rm = TRUE)
            }
            for (j in 1:(p - 1)) dat[, j + 1] = (dat[, j + 1] - 
                m[j])/sd[j]
            return(dat)
        }
        data = dat.sta(data)
        return(data)
    }
    dat.rem.ind.pop <- function(data, ind = 0, pop = 0) {
        data = as.data.frame(data)
        dat.rem.ind <- function(dat, ind) {
            l = length(ind)
            n = dim(dat)[1]
            for (i in 1:l) dat = dat[row.names(dat)[1:(n - i + 
                1)] != ind[i], ]
            return(dat)
        }
        dat.rem.pop <- function(dat, pop) {
            l = length(pop)
            for (i in 1:l) dat = dat[dat[, 1] != pop[i], ]
            return(dat)
        }
        if (ind[1] != 0) 
            data = dat.rem.ind(data, ind)
        if (pop[1] != 0) 
            data = dat.rem.pop(data, pop)
        return(data)
    }
    dat.pw <- function(data, pw = 0) {
        if (pw[1] == 0) 
            return(data)
        else {
            data = data[data[, 1] == pw[1] | data[, 1] == pw[2], 
                ]
            return(data)
        }
    }
    nb.pop <- function(data) {
        data = data[order(data[, 1]), ]
        n = dim(data)[1]
        c = 1
        for (i in 1:(n - 1)) if (data[i, 1] != data[i + 1, 1]) 
            c = c + 1
        return(c)
    }
    Pop <- function(data) {
        nb.ind = dim(data)[1]
        c = 1
        dat.fra = as.data.frame(data)
        dat.fra = dat.fra[order(dat.fra[, 1]), ]
        for (i in 1:(nb.ind - 1)) if (dat.fra[i, 1] != dat.fra[i + 
            1, 1]) 
            c = c + 1
        vec = rep(1, c)
        name = rep(0, c)
        k = 1
        l = 2
        name[1] = as.character(dat.fra[1, 1])
        for (i in 2:nb.ind) if (dat.fra[i - 1, 1] == dat.fra[i, 
            1]) 
            vec[k] = vec[k] + 1
        else {
            name[l] = as.character(dat.fra[i, 1])
            l = l + 1
            k = k + 1
        }
        names(vec) = name
        return(vec)
    }
    Pst.val <- function(data, csh = 1) {
        data = Prep(data)
        nbpop = nb.pop(data)
        va = dim(data)[2]
        data = data[order(data[, 1]), ]
        if (nbpop == 1) 
            return(rep(0, va - 1))
        else {
            v = Pop(data)
            Pst.clm <- function(dat, clm) {
                m = mean(dat[, clm], na.rm = TRUE)
                nna.clm = nonNa.clm(dat, clm)
                Sst = (nna.clm - 1) * var(dat[, clm], na.rm = TRUE)
                mp = rep(0, nbpop)
                vef = rep(0, nbpop)
                vef[1] = nonNa.clm(dat[1:(v[1]), ], clm)
                q = 0
                if (vef[1] == 0) 
                  q = 1
                else mp[1] = mean(dat[1:(v[1]), clm], na.rm = TRUE)
                for (i in 2:nbpop) {
                  vef[i] = nonNa.clm(dat[(sum(v[1:(i - 1)]) + 
                    1):(sum(v[1:i])), ], clm)
                  if (vef[i] != 0) 
                    mp[i] = mean(dat[(sum(v[1:(i - 1)]) + 1):(sum(v[1:i])), 
                      clm], na.rm = TRUE)
                  else q = q + 1
                }
                Ssb = sum(vef * (mp - m)^2)
                Ssw = Sst - Ssb
                if ((nna.clm - nbpop + q) * (nbpop - q - 1) != 
                  0) {
                  Msw = Ssw/(nna.clm - nbpop + q)
                  Msb = Ssb/(nbpop - q - 1)
                  return(csh * Msb/(csh * Msb + 2 * Msw))
                }
                else {
                  if ((nna.clm - nbpop + q) == 0) 
                    return(1)
                  else return(0)
                }
            }
            pst = rep(0, va - 1)
            for (j in 1:(va - 1)) pst[j] = Pst.clm(data, j + 
                1)
            return(pst)
        }
    }
    boot.pst.va <- function(data, csh, boot, clm) {
        data = data[, c(1, clm)]
        data = Prep(data)
        v = rep(0, boot)
        n = dim(data)[1]
        nbpop = nb.pop(data)
        Psts <- function(dat, csh) {
            dat = dat[order(dat[, 1]), ]
            Po = nb.pop(dat)
            if (Po == 1) 
                return(0)
            else {
                vec = c(Pop(dat), rep(0, nbpop - Po))
                m = mean(dat[, 2], na.rm = TRUE)
                nna.clm = nonNa.clm(dat, 2)
                Sst = (nna.clm - 1) * var(dat[, 2], na.rm = TRUE)
                mp = rep(0, nbpop)
                vef = rep(0, nbpop)
                c = 0
                vef[1] = nonNa.clm(dat[1:(vec[1]), ], 2)
                if (vef[1] == 0) 
                  c = 1
                else mp[1] = mean(dat[1:(vec[1]), 2], na.rm = TRUE)
                for (i in 2:Po) {
                  vef[i] = nonNa.clm(dat[(sum(vec[1:(i - 1)]) + 
                    1):(sum(vec[1:i])), ], 2)
                  if (vef[i] != 0) 
                    mp[i] = mean(dat[(sum(vec[1:(i - 1)]) + 1):(sum(vec[1:i])), 
                      2], na.rm = TRUE)
                  else c = c + 1
                }
                Ssb = sum(vef * (mp - m)^2)
                Ssw = Sst - Ssb
                if ((Po - c - 1) * (nna.clm + c - Po) != 0) {
                  Msb = Ssb/(Po - c - 1)
                  Msw = Ssw/(nna.clm + c - Po)
                  return(csh * Msb/(csh * Msb + 2 * Msw))
                }
                else {
                  if ((nna.clm + c - Po) == 0) 
                    return(1)
                  else return(0)
                }
            }
        }
        for (i in 1:boot) {
            da = data[sample(1:n, n, T), ]
            v[i] = Psts(da, csh)
        }
        return(v)
    }
    ConInt.pst.va <- function(data, csh, boot, clm, per) {
        v = boot.pst.va(data = data, csh = csh, boot = boot, 
            clm = clm)
        v = sort(v)
        return(c(v[floor(boot * (1 - per)/2 + 1)], v[ceiling(boot * 
            (per + 1)/2)]))
    }
    l = length(va)
    for (i in 1:l) {
        for (j in 2:dim(data)[2]) {
            if (names(data)[j] == va[i]) 
                va[i] = j - 1
        }
    }
    va = as.numeric(va)
    if (is.na(sum(va)) == TRUE) 
        return("va is not valid!")
    data = Prep(data)
    data = dat.rem.ind.pop(data, ind = Ri, pop = Rp)
    data = dat.pw(data, Pw)
    print("Populations sizes are:")
    print(Pop(data))
    if (ci != 1) {
        if (va[1] == 0) {
            Q = dim(data)[2]
            D = data.frame(Quant_Varia = names(data)[2:Q], Pst_Values = Pst.val(data, 
                csh = csh), row.names = NULL)
        }
        else {
            D = data.frame(Quant_Varia = names(data)[va + 1], 
                Pst_Values = Pst.val(data, csh = csh)[va], row.names = NULL)
        }
        return(D)
    }
    if (ci == 1) {
        if (va[1] == 0) {
            Q = dim(data)[2]
            v = rep(0, Q - 1)
            w = rep(0, Q - 1)
            for (i in 1:(Q - 1)) {
                v[i] = ConInt.pst.va(data, csh = csh, boot = boot, 
                  clm = i + 1, per = pe)[1]
                w[i] = ConInt.pst.va(data, csh = csh, boot = boot, 
                  clm = i + 1, per = pe)[2]
            }
            dat = data.frame(Quant_Varia = names(data)[2:Q], 
                Pst_Values = Pst.val(data, csh = csh), LowBoundCI = v, 
                UpBoundCI = w, row.names = NULL)
            names(dat)[3] = paste(100 * pe, "\%_LowBoundCI")
            names(dat)[4] = paste(100 * pe, "\%_UpBoundCI")
        }
        else {
            n = length(va)
            v = rep(0, n)
            w = rep(0, n)
            for (i in 1:n) {
                v[i] = ConInt.pst.va(data, csh = csh, boot = boot, 
                  clm = va[i] + 1, per = pe)[1]
                w[i] = ConInt.pst.va(data, csh = csh, boot = boot, 
                  clm = va[i] + 1, per = pe)[2]
            }
            dat = data.frame(Quant_Varia = names(data)[va + 1], 
                Pst_Values = Pst.val(data, csh = csh)[va], LowBoundCI = v, 
                UpBoundCI = w, row.names = NULL)
            names(dat)[3] = paste(100 * pe, "\%_LowBoundCI")
            names(dat)[4] = paste(100 * pe, "\%_UpBoundCI")
        }
        return(dat)
    }
  }
}

