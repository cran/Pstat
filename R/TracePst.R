TracePst <-
function(data,va=0,ci=1,boot=1000,pe=0.95,Fst=-1,Pw=0,Rp=0,Ri=0,xm=2,pts=30){

# function returning the number of values differents from Na in a selected column 

nonNa.clm<-function(data,clm){
n=dim(data)[1];c=0;
for (i in 1:n) if (is.na(data[i,clm])) c=c+1;
return(n-c);}

# function preparing the data frame: the first column being a column of characters, the others being numerical columns 

Prep<-function(data){
c=dim(data)[2]-1;data=as.data.frame(data);data[,1]=as.character(data[,1]);
for (i in 1:c) {if(is.numeric(data[,i+1])==FALSE) data[,i+1]=as.numeric(as.character(data[,i+1]))};
# function standardizing columns (except the first one): Pst values are not impacted by such a transformation
dat.sta<-function(dat){
p=dim(dat)[2];sd=rep(0,p-1);m=rep(0,p-1);
for (i in 1:(p-1)) {N=nonNa.clm(dat,i+1);sd[i]=sqrt((N-1)/N)*sd(dat[,i+1], na.rm=TRUE); m[i]=mean(dat[,i+1],na.rm=TRUE);}
for (j in 1:(p-1)) dat[,j+1]=(dat[,j+1]-m[j])/sd[j];
return(dat);}
data=dat.sta(data); return(data);}

# function removing unwanted individuals or populations 

dat.rem.ind.pop<-function(data,ind=0,pop=0){
data=as.data.frame(data);
dat.rem.ind<-function(dat,ind){
l=length(ind);n=dim(dat)[1];
for(i in 1:l) dat=dat[row.names(dat)[1:(n-i+1)]!=ind[i],];
return(dat)};
dat.rem.pop<-function(dat,pop){
l=length(pop);
for(i in 1:l) dat=dat[dat[,1]!=pop[i],]; 
return(dat);}
if (ind[1]!=0) data=dat.rem.ind(data,ind);
if (pop[1]!=0) data=dat.rem.pop(data,pop);
return(data);}

# function selecting two populations to make pairwise comparisons

dat.pw<-function(data,pw=0){
if (pw[1]==0) return(data)
else {data=data[data[,1]==pw[1] | data[,1]==pw[2],]; return(data)}}

# function returning the number of populations

nb.pop<-function(data){
data=data[order(data[,1]),];n=dim(data)[1]; c=1; 
for(i in 1:(n-1)) if(data[i,1]!=data[i+1,1]) c=c+1; 
return(c);}

# function returning a vector containing the number of individuals of each population 

Pop<-function(data){
nb.ind=dim(data)[1];c=1;dat.fra=as.data.frame(data);dat.fra=dat.fra[order(dat.fra[,1]),];
for(i in 1:(nb.ind-1)) if(dat.fra[i,1]!=dat.fra[i+1,1]) c=c+1;
vec=rep(1,c);name=rep(0,c);k=1;l=2;name[1]=as.character(dat.fra[1,1]);
for(i in 2:nb.ind) if(dat.fra[i-1,1]==dat.fra[i,1]) vec[k]=vec[k]+1 
else {name[l]=as.character(dat.fra[i,1]);l=l+1;k=k+1;}
names(vec)=name;
return(vec);}

# function returning a vector containing the Pst values of the data frame variables  

Pst.val<-function(data,csh=1){
data=Prep(data);nbpop=nb.pop(data);va=dim(data)[2];data=data[order(data[,1]),];
if(nbpop==1) return(rep(0,va-1))
else {v=Pop(data);
# function returning Pst value of the selected quantitative variable  of a data frame
Pst.clm<-function(dat,clm){
m=mean(dat[,clm],na.rm=TRUE); nna.clm=nonNa.clm(dat,clm);
Sst=(nna.clm-1)*var(dat[,clm],na.rm=TRUE);mp=rep(0,nbpop); vef=rep(0,nbpop); vef[1]=nonNa.clm(dat[1:(v[1]),],clm); q=0;
if (vef[1]==0) q=1 
else mp[1]=mean(dat[1:(v[1]),clm],na.rm=TRUE); 
for (i in 2:nbpop) {vef[i]=nonNa.clm(dat[(sum(v[1:(i-1)])+1):(sum(v[1:i])),],clm);
if(vef[i]!=0) mp[i]=mean(dat[(sum(v[1:(i-1)])+1):(sum(v[1:i])),clm],na.rm=TRUE) else q=q+1;}
Ssb=sum(vef*(mp-m)^2);Ssw=Sst-Ssb;
if ((nna.clm-nbpop+q)*(nbpop-q-1)!=0) {Msw=Ssw/(nna.clm-nbpop+q);Msb=Ssb/(nbpop-q-1);
return(csh*Msb/(csh*Msb+2*Msw));}
else {if ((nna.clm-nbpop+q)==0) return(1)
else return(0);};} 
pst=rep(0,va-1);for(j in 1:(va-1)) pst[j]=Pst.clm(data,j+1); 
return(pst);}}

# function returning a vector containing the Pst values of the bootstrap data frames built from the selected quantitative variable

boot.pst.va<-function(data,csh,boot,clm){
data=data[,c(1,clm)];data=Prep(data);v=rep(0,boot);n=dim(data)[1];nbpop=nb.pop(data);
# function returning the Pst value of a bootstrap data frame containing only the selected variable. Attention should be paid to data frames with a single population and also to populations not containing enough non-Na values
Psts<-function(dat,csh){
dat=dat[order(dat[,1]),];Po=nb.pop(dat);
if (Po==1) return(0)
else {vec=c(Pop(dat),rep(0,nbpop-Po)); m=mean(dat[,2],na.rm=TRUE);nna.clm=nonNa.clm(dat,2);
Sst=(nna.clm-1)*var(dat[,2],na.rm=TRUE);mp=rep(0,nbpop); vef=rep(0,nbpop); c=0;vef[1]=nonNa.clm(dat[1:(vec[1]),],2);
if (vef[1]==0) c=1 
else mp[1]=mean(dat[1:(vec[1]),2],na.rm=TRUE); 
for (i in 2:Po) {vef[i]=nonNa.clm(dat[(sum(vec[1:(i-1)])+1):(sum(vec[1:i])),],2);
if(vef[i]!=0) mp[i]=mean(dat[(sum(vec[1:(i-1)])+1):(sum(vec[1:i])),2],na.rm=TRUE) else c=c+1;}
Ssb=sum(vef*(mp-m)^2);Ssw=Sst-Ssb;
if ((Po-c-1)*(nna.clm+c-Po)!=0) {Msb=Ssb/(Po-c-1); Msw=Ssw/(nna.clm+c-Po); 
return(csh*Msb/(csh*Msb+2*Msw));}
else {if ((nna.clm+c-Po)==0) return(1)
else return(0);};}}
for (i in 1:boot) {da=data[sample(1:n,n,T),];v[i]=Psts(da,csh);}
return(v);}

# function providing the confidence interval associated with bootstrap values of Pst

ConInt.pst.va<-function(data,csh,boot,clm,per){
v=boot.pst.va(data=data,csh=csh,boot=boot,clm=clm);v=sort(v); 
return(c(v[floor(boot*(1-per)/2+1)],v[ceiling(boot*(per+1)/2)]));}

# function drawing the chosen plots

Trace<-function(data,pts,boot,Fst,xm=2,ci){
# function drawing the curve associated with Pst values of a selected variable (and also the Fst line, if requested)
tra.pst.va<-function(data,pts=30,clm,Fst,xm=2){
data=data[,c(1,clm)];
points<-function(n){
P=Pst.val(data,0); 
for (i in 1:n) P=c(P,Pst.val(data,xm*i/n));return(P);}
p=points(pts);d=xm*c(0:pts)/pts;
plot(p~d,type="l",xlab="c/h^2",ylab="Pst",main=c("Pst variations:",names(data)[2]),ylim=c(0,1),col="firebrick1");
if (Fst!=-1) {
abline(h=Fst,col="green",lty=4);text(0.05*a-0.06,Fst+0.04*a-0.01,"Fst",col="green");};}
# function drawing the two dotted curves associated with confidence intervals of a selected variable
tra.confint.va<-function(clm){
point<-function(n){
R=ConInt.pst.va(data,csh=0,boot=boot,clm=clm,per=pe);P=R[2];Q=R[1];
for (i in 1:n) R=c(R,ConInt.pst.va(data,csh=xm*i/n,boot=boot,clm=clm,per=pe));
for (i in 1:n) P=c(R[2+2*i],P);
for (i in 1:n) Q=c(Q,R[1+2*i]);return(c(P,Q));} 
p=point(pts);d=xm*c(0:pts)/pts;e=rev(d);
plot(p~c(e,d),type="l",xlab="c/h^2",ylab="Pst",main=c("Pst variations:",names(data)[clm]), ylim=c(0,1),col="chocolate4",lty=2)}
Q=dim(data)[2];a=ceiling(sqrt(Q-1));par(mfrow=c(a,a));
for (i in 2:Q) {tra.pst.va(data,pts=pts,Fst=Fst,clm=i,xm=xm); if (ci==1){par(new=TRUE);tra.confint.va(clm=i);};}}

# to obtain the output:

l=length(va);for (i in 1:l) {for (j in 2:dim(data)[2]) {if (names(data)[j]==va[i]) va[i]=j-1}};
va=as.numeric(va);if (is.na(sum(va))==TRUE) return("va is not valid!");
data=Prep(data);data=dat.rem.ind.pop(data,ind=Ri,pop=Rp);data=dat.pw(data,Pw);
print("Populations sizes are:"); print(Pop(data)); 
if (va[1]==0) {dev.new(); Trace(data,pts=pts,boot=boot,Fst=Fst,xm=xm,ci=ci)}
else {data=data[,c(1,va+1)]; dev.new(); Trace(data,pts=pts,boot=boot,Fst=Fst,xm=xm,ci=ci);};}
