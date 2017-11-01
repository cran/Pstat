Res <-
function(data,reg){

# function returning the residuals of a simple linear regression; a quantitative variable is transformed, the regressor being previously chosen by the user

Res.va<-function(dat,clm,re){
dat=dat[is.finite(dat[,re]),];m.clm=mean(dat[is.finite(dat[,clm]),clm],na.rm=TRUE);
m.reg=mean(dat[is.finite(dat[,clm]),re],na.rm=TRUE);
a=sum((dat[is.finite(dat[,clm]),re]-m.reg)*dat[is.finite(dat[,clm]),clm])/sum((dat[is.finite(dat[,clm]),re]-m.reg)*(dat[is.finite(dat[,clm]),re]-m.reg));
b=m.clm-a*m.reg;
dat[,clm]=dat[,clm]-b-a*dat[,re];
return(dat);}

for (i in 2:dim(data)[2]) {if (names(data)[i]==reg) reg=i-1};
if (is.numeric(reg)==FALSE) return("reg value does not exist!");
Q=dim(data)[2]; 
if (reg==1) for (i in 3:Q) data=Res.va(data,clm=i,re=2)
else {for (i in 2:reg) data=Res.va(data,clm=i,re=reg+1); 
if (reg<(Q-1)) for (j in (reg+2):Q) data=Res.va(data,clm=j,re=reg+1);}
return(data[-(reg+1)]);}
