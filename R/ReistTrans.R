ReistTrans <-
function(M,r){
ReistTransk<-function(M,k,r){M=M[is.finite(M[,r]),];L=M;m=mean(M[is.finite(M[,k]),r]);Q=dim(M)[2];
L[,k]=log(M[,k],base=10);L[,r]=log(M[,r],base=10);
mk=mean(L[is.finite(L[,k]),k],na.rm=TRUE);
mr=mean(L[is.finite(L[,k]),r],na.rm=TRUE);
a=sum((L[is.finite(L[,k]),r]-mr)*L[is.finite(L[,k]),k])/sum((L[is.finite(L[,k]),r]-mr)*(L[is.finite(L[,k]),r]-mr));
M[,k]=L[,k]-a*(L[,r]-log(m,base=10));
return(M);}
Q=dim(M)[2]; if (r==2) for (i in 3:Q) M=ReistTransk(M,i,2)
else {for (i in 2:(r-1)) M=ReistTransk(M,i,r); if (r<Q) for (j in (r+1):Q) M=ReistTransk(M,j,r);};
return(M[-r]);}
