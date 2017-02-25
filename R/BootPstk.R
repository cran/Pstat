BootPstk <-
function(M,b=1,csh=1,boot=1000,k,Ri=0,Rp=0,Pw=0,pe=0.95,bars=20){
Prep<-function(M){
Nk<-function(M,j){n=dim(M)[1];c=0;
for (i in 1:n) {if (is.na(M[i,j])) c=c+1;}
return(n-c);}
c=dim(M)[2]-1;M=as.data.frame(M);
M[,1]=as.character(M[,1]);
for (i in 1:c) {if(is.numeric(M[,i+1])==FALSE) M[,i+1]=as.numeric(as.character(M[,i+1]))};
Mst<-function(Ma){p=dim(Ma)[2];sd=rep(0,p-1);m=rep(0,p-1);
for (i in 1:(p-1)) {N=Nk(Ma,i+1); sd[i]=sqrt((N-1)/N)*sd(M[,i+1],na.rm=TRUE);m[i]=mean(Ma[,i+1],na.rm=TRUE);}
for (j in 1:(p-1)) Ma[,j+1]=(Ma[,j+1]-m[j])/sd[j];
return(Ma);}
M=Mst(M); return(M);}
MaReIP<-function(M,I=0,P=0){Ma=as.data.frame(M);
MaReI<-function(Mat,I){l=length(I);n=dim(Mat)[1];
for(i in 1:l) Mat=Mat[row.names(Mat)[1:(n-i+1)]!=I[i],];return(Mat)};
MaReP<-function(Mat,P){l=length(P);
for(i in 1:l) Mat=Mat[Mat[,1]!=P[i],]; 
return(Mat);}
if (I[1]!=0) Ma=MaReI(Ma,I);
if (P[1]!=0) Ma=MaReP(Ma,P);
return(Ma);}
MaPw<-function(M,Pw=0){if (Pw[1]==0) return(M)
else {M=M[M[,1]==Pw[1] | M[,1]==Pw[2],]; return(M)}}
Pop<-function(M){n=dim(M)[1];c=1; Ma=as.data.frame(M);Ma=Ma[order(Ma[,1]),];
for(i in 1:(n-1)) if(Ma[i,1]!=Ma[i+1,1]) c=c+1;
vec=rep(1,c);name=rep(0,c);l=1;c=2;name[1]=as.character(Ma[1,1]);
for(i in 2:n) if(Ma[i-1,1]==Ma[i,1]) vec[l]=vec[l]+1 else {name[c]=as.character(Ma[i,1]);c=c+1;l=l+1;}
names(vec)=name;
return(vec);}
Nk<-function(M,j){n=dim(M)[1];c=0;
for (i in 1:n) {if (is.na(M[i,j])) c=c+1;}
return(n-c);}
NbPop<-function(M){M=M[order(M[,1]),];
n=dim(M)[1]; c=1; 
for(i in 1:(n-1)) if(M[i,1]!=M[i+1,1]) c=c+1; 
return(c);}
Bootk<-function(M,csh,boot,k){
M=M[,c(1,k)];M=Prep(M);
v=rep(0,boot);n=dim(M)[1];P=NbPop(M);
Psts<-function(Ma,csh){Ma=Ma[order(Ma[,1]),];Po=NbPop(Ma);
if(Po==1) return(0)
else {vec=c(Pop(Ma),rep(0,P-Po));
Pstks<-function(Man){m=mean(Man[,2],na.rm=TRUE);nk=Nk(Man,2);
Sst=(nk-1)*var(Man[,2],na.rm=TRUE);
mp=rep(0,P); vef=rep(0,P); c=0; vef[1]=Nk(Man[1:(vec[1]),],2);
if (vef[1]==0) c=1 else mp[1]=mean(Man[1:(vec[1]),2],na.rm=TRUE); 
for (i in 2:Po) {vef[i]=Nk(Man[(sum(vec[1:(i-1)])+1):(sum(vec[1:i])),],2);
if(vef[i]!=0) mp[i]=mean(Man[(sum(vec[1:(i-1)])+1):(sum(vec[1:i])),2],na.rm=TRUE) else c=c+1;}
Ssb=sum(vef*(mp-m)^2);Ssw=Sst-Ssb;
if ((Po-c-1)*(nk+c-Po)!=0) {Msb=Ssb/(Po-c-1); Msw=Ssw/(nk+c-Po); 
return(csh*Msb/(csh*Msb+2*Msw));}
else {if ((nk+c-Po)==0) return(1)
else return(0);};}
return(Pstks(Ma));}}
for (i in 1:boot) {B=M[sample(1:n,n,T),];v[i]=Psts(B,csh);}
return(v);}
M=MaReIP(M,Ri,Rp);M=MaPw(M,Pw);
if (b==1) {print("Populations sizes are:"); print(Pop(M)); print(Bootk(M,csh,boot,k));}
if (b==2) {print("Populations sizes are:"); print(Pop(M));
ConfinterPstk<-function(M,csh,boot,k,per){
v=Bootk(M,csh,boot,k);v=sort(v); print(v);
return(c(v[floor(boot*(1-per)/2+1)],v[ceiling(boot*(per+1)/2)]));}
print(ConfinterPstk(M,csh,boot,k,per=pe))};
if (b==3) {print("Populations sizes are:"); print(Pop(M));
DistPstk<-function(M,csh,boot,k,bars){psts=Bootk(M,csh,boot,k);
hist(psts,breaks=c(0:bars)/bars,xlab="Pst",ylab="Frequency",main=c("Pst distribution:",names(M)[k]),col="gray88"); return(sort(psts));}
dev.new(); DistPstk(M,csh,boot,k,bars);}}
