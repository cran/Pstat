AitTrans <-
function(data){

# function returning the number of values different from Na in a selected column 

nonNa.clm<-function(data,clm){
n=dim(data)[1];c=0;
for (i in 1:n) if (is.na(data[i,clm])) c=c+1;
return(n-c);}

nb.ind=dim(data)[1];nb.va=dim(data)[2]-1;v=rep(0,nb.ind);
for (i in 1:nb.ind) {v[i]=1/(nonNa.clm(t(data[i,-1]),1))*sum(log(data[i,-1], base=10),na.rm=TRUE);
for (j in 1:nb.va) data[i,(j+1)]=log(data[i,(j+1)],base=10)-v[i];}
return(data);}
