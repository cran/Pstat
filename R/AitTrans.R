AitTrans <-
function(M){
Nk<-function(M,j){n=dim(M)[1];c=0;
for (i in 1:n) {if (is.na(M[i,j])) c=c+1;}
return(n-c);}
N=dim(M)[1];Q=dim(M)[2]-1;
v=rep(0,N);
for (i in 1:N) {v[i]=1/(Nk(t(M[i,-1]),1))*sum(log(M[i,-1],base=10),na.rm=TRUE);
for (j in 1:Q) M[i,(j+1)]=log(M[i,(j+1)],base=10)-v[i];}
return(M);}
