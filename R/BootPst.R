BootPst <-
function(data,va,opt=0,csh=1,boot=1000,Ri=0,Rp=0,Pw=0,pe=0.95,bars=20){

# function returning the number of values different from Na in a selected column 

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
print(c(v[floor(boot*(1-per)/2+1)],v[ceiling(boot*(per+1)/2)]));return(v);}

# function drawing the histogram of Pst ditribution (from the bootstrap values)

dis.pst.va<-function(data,csh,boot,clm,bars){
psts.va=boot.pst.va(data=data,csh=csh,boot=boot,clm=clm);
hist(psts.va,breaks=c(0:bars)/bars,xlab="Pst",ylab="Frequency",main=c("Pst distribution:",names(data)[clm]),col="gray88");
return(sort(psts.va));}

# to obtain the output:

for (i in 2:dim(data)[2]) {if (names(data)[i]==va) va=i-1};
if (is.numeric(va)==FALSE) return("va value does not exist!");
data=dat.rem.ind.pop(data,ind=Ri,pop=Rp);data=dat.pw(data,pw=Pw);
print("Populations sizes are:"); print(Pop(data)); 
if (opt!="ci" & opt!="hist") {print(paste(boot,"bootstrap values:")); return(boot.pst.va(data,csh=csh,boot=boot,clm=va+1));}
if (opt=="ci") {print(paste(100*pe,"% confidence interval determined by",boot,"bootstrap values:")); return(ConInt.pst.va(data,csh=csh,boot=boot,clm=va+1,per=pe));};
if (opt=="hist") {print(paste(boot,"bootstrap values and","Pst distribution:"));dev.new(); dis.pst.va(data,csh,boot,va+1,bars);}}
