Res <-
function(M,r){
Resk<-function(M,k,r){M=M[is.finite(M[,r]),];
mk=mean(M[is.finite(M[,k]),k],na.rm=TRUE);
mr=mean(M[is.finite(M[,k]),r],na.rm=TRUE);
a=sum((M[is.finite(M[,k]),r]-mr)*M[is.finite(M[,k]),k])/sum((M[is.finite(M[,k]),r]-mr)*(M[is.finite(M[,k]),r]-mr));
b=mk-a*mr;
M[,k]=M[,k]-b-a*M[,r];
return(M);}
Q=dim(M)[2]; if (r==2) for (i in 3:Q) M=Resk(M,i,2)
else {for (i in 2:(r-1)) M=Resk(M,i,r); if (r<Q) for (j in (r+1):Q) M=Resk(M,j,r);}
return(M[-r]);}
