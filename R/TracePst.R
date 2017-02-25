TracePst <-
function(Mat,col=1,ci=1,boot=1000,pe=0.95,Fst=-1,Pw=0,Rp=0,Ri=0,xm=2,pts=30){
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
Nk<-function(M,j){n=dim(M)[1];c=0;
for (i in 1:n) {if (is.na(M[i,j])) c=c+1;}
return(n-c);}
NbPop<-function(M){M=M[order(M[,1]),];
n=dim(M)[1]; c=1; 
for(i in 1:(n-1)) if(M[i,1]!=M[i+1,1]) c=c+1; 
return(c);}
Pop<-function(M){n=dim(M)[1];c=1; Ma=as.data.frame(M);Ma=Ma[order(Ma[,1]),];
for(i in 1:(n-1)) if(Ma[i,1]!=Ma[i+1,1]) c=c+1;
vec=rep(1,c);name=rep(0,c);l=1;c=2;name[1]=as.character(Ma[1,1]);
for(i in 2:n) if(Ma[i-1,1]==Ma[i,1]) vec[l]=vec[l]+1 else {name[c]=as.character(Ma[i,1]);c=c+1;l=l+1;}
names(vec)=name;
return(vec);}
Pstc<-function(M,csh=1){
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
Nk<-function(M,j){n=dim(M)[1];c=0;
for (i in 1:n) {if (is.na(M[i,j])) c=c+1;}
return(n-c);}
Pop<-function(M){n=dim(M)[1];c=1; Ma=as.data.frame(M);Ma=Ma[order(Ma[,1]),];
for(i in 1:(n-1)) if(Ma[i,1]!=Ma[i+1,1]) c=c+1;
vec=rep(1,c);name=rep(0,c);l=1;c=2;name[1]=as.character(Ma[1,1]);
for(i in 2:n) if(Ma[i-1,1]==Ma[i,1]) vec[l]=vec[l]+1 else {name[c]=as.character(Ma[i,1]);c=c+1;l=l+1;}
names(vec)=name;
return(vec);}
NbPop<-function(M){M=M[order(M[,1]),];
n=dim(M)[1]; c=1; 
for(i in 1:(n-1)) if(M[i,1]!=M[i+1,1]) c=c+1; 
return(c);}
M=Prep(M);
P=NbPop(M);p=dim(M)[2];M=M[order(M[,1]),];
if(P==1) return(rep(0,p-1))
else {v=Pop(M);
Pstk<-function(M,k){m=mean(M[,k],na.rm=TRUE);nk=Nk(M,k);
Sst=(nk-1)*var(M[,k],na.rm=TRUE);
mp=rep(0,P); vef=rep(0,P); vef[1]=Nk(M[1:(v[1]),],k); q=0;
if (vef[1]==0) q=1 else mp[1]=mean(M[1:(v[1]),k],na.rm=TRUE); 
for (i in 2:P) {vef[i]=Nk(M[(sum(v[1:(i-1)])+1):(sum(v[1:i])),],k);
if(vef[i]!=0) mp[i]=mean(M[(sum(v[1:(i-1)])+1):(sum(v[1:i])),k],na.rm=TRUE) else q=q+1;}
Ssb=sum(vef*(mp-m)^2);Ssw=Sst-Ssb;
if ((nk-P+q)*(P-q-1)!=0) {Msw=Ssw/(nk-P+q);Msb=Ssb/(P-q-1);
return(csh*Msb/(csh*Msb+2*Msw));}
else {if ((nk-P+q)==0) return(1)
else return(0);};} 
pst=rep(0,p-1);for(j in 1:(p-1)) pst[j]=Pstk(M,j+1); 
return(pst);}}
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

ConfinterPstk<-function(M,csh,boot,k,per){
v=Bootk(M,csh=csh,boot=boot,k=k);v=sort(v);
return(c(v[floor(boot*(1-per)/2+1)],v[ceiling(boot*(1+per)/2)]));};
TraceCIk<-function(k){point<-function(n){R=ConfinterPstk(A,csh=0,boot=boot,k=k,per=pe);P=R[2];Q=R[1];
for (i in 1:n) R=c(R,ConfinterPstk(A,csh=xm*i/n,boot=boot,k=k,per=pe));
for (i in 1:n) P=c(R[2+2*i],P);
for (i in 1:n) Q=c(Q,R[1+2*i]);return(c(P,Q));} 
p=point(pts);d=xm*c(0:pts)/pts;e=rev(d);
plot(p~c(e,d),type="l",xlab="c/h^2",ylab="Pst",main=c("Pst variations:",names(A)[k]),ylim=c(0,1),col="chocolate4",lty=2)}
Trace<-function(A,pts=30,boot=1000,xm=2,ci){
TracePstk<-function(A,pts=30,k,xm=2){A=A[,c(1,k)];points<-function(n){P=Pstc(A,0); 
for (i in 1:n) P=c(P,Pstc(A,xm*i/n));return(P);} 
p=points(pts);d=xm*c(0:pts)/pts;
plot(p~d,type="l",xlab="c/h^2",ylab="Pst",main=c("Pst variations:",names(A)[2]),ylim=c(0,1),col="firebrick1")}
TracePstCIk<-function(A,pts=30,boot=1000,k,xm=2){
TracePstk(A,pts=pts,k=k,xm=xm); if (ci==1) {par(new=TRUE);TraceCIk(k);}}
Q=dim(A)[2];a=ceiling(sqrt(Q-1));par(mfrow=c(a,a));
for (i in 2:Q) TracePstCIk(A,pts=pts,boot=boot,k=i,xm=xm);}
TraceF<-function(A,pts=30,boot,Fst,xm=2,ci){
TraceFQk<-function(A,pts=30,k,Fst,xm=2){points<-function(n){P=Pstc(A,0)[k-1]; 
for (i in 1:n) P=c(P,Pstc(A,xm*i/n)[k-1]);return(P);} 
p=points(pts);d=xm*c(0:pts)/pts;a=ceiling(sqrt(dim(A)[2]-1));
plot(p~d,type="l",xlab="c/h^2",ylab="Pst",main=c("Pst variations:",names(A)[k]),ylim=c(0,1),col="firebrick1");
abline(h=Fst,col="green",lty=4);text(0.05*a-0.06,Fst+0.04*a-0.01,"Fst",col="green");}
TraceFQCIk<-function(A,pts=30,boot=1000,k,Fst,xm=2){
TraceFQk(A,pts=pts,k=k,Fst=Fst,xm=xm); if (ci==1) {par(new=TRUE);TraceCIk(k);}}
Q=dim(A)[2];a=ceiling(sqrt(Q-1));par(mfrow=c(a,a));
for (i in 2:Q) TraceFQCIk(A,pts=pts,boot=boot,k=i,Fst=Fst,xm=xm);}
A=Prep(Mat);
A=MaReIP(A,Ri,Rp);A=MaPw(A,Pw);
if (Fst==-1) {print("Populations sizes are:"); print(Pop(A)); if (col[1]==1) {dev.new(); Trace(A,pts=pts,boot=boot,xm=xm,ci=ci)}
else {A=A[,c(1,col)]; dev.new(); Trace(A,pts=pts,boot=boot,xm=xm,ci=ci);};};
if (Fst!=-1) {print("Populations sizes are:"); print(Pop(A)); dev.new(); if (col[1]==1) TraceF(A,pts=pts,boot=boot,Fst=Fst,xm=xm,ci=ci)
else {A=A[,c(1,col)]; dev.new(); TraceF(A,pts=pts,boot=boot,Fst=Fst,xm=xm,ci=ci);};};}
