ReistTrans <-
function(data,reg){

# function returning one of the quantitative variables transformed following Reist's transformation,  the regressor being previously chosen by the user

Reitra.va<-function(dat,clm,re){
dat=dat[is.finite(dat[,re]),];log.dat=dat;m=mean(dat[is.finite(dat[,clm]),re]);Q=dim(dat)[2];
log.dat[,clm]=log(dat[,clm],base=10);log.dat[,re]=log(dat[,re],base=10);
m.clm=mean(log.dat[is.finite(log.dat[,clm]),clm],na.rm=TRUE);
m.reg=mean(log.dat[is.finite(log.dat[,clm]),re],na.rm=TRUE);
a=sum((log.dat[is.finite(log.dat[,clm]),re]-m.reg)*log.dat[is.finite(log.dat[,clm]),clm])/sum((log.dat[is.finite(log.dat[,clm]),re]-m.reg)*(log.dat[is.finite(log.dat[,clm]),re]-m.reg));
dat[,clm]=log.dat[,clm]-a*(log.dat[,re]-log(m,base=10));
return(dat);}

for (i in 2:dim(data)[2]) {if (names(data)[i]==reg) reg=i-1};
if (is.numeric(reg)==FALSE) return("reg value does not exist!");
Q=dim(data)[2]; 
if (reg==1) for (i in 3:Q) data=Reitra.va(data,clm=i,re=2)
else {for (i in 2:reg) data=Reitra.va(data,clm=i,re=reg+1); 
if (reg<(Q-1)) for (j in (reg+2):Q) data=Reitra.va(data,clm=j,re=reg+1);};
return(data[-(reg+1)]);}
