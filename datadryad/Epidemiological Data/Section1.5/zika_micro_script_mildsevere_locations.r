
#zika_analysis by location for nominally 53 places with zika>100

rm(list=ls())
setwd("~/desktop/")

zm=read.table("new_Zika_ms.txt",header=F)
micros=read.table("new_micros_ms.txt",header=F)

pdf("mildsevere_zika_micro_location.pdf")
par(mfrow=c(1,2))
plot(seq(0,50,length=10),seq(0,20,length=10),type="n",bty="l",xlab="Time (Weeks)",ylab="sqrt(Cases)",main="Zika cases")
for(i in 1:53){
    lines(seq(1:49), sqrt(as.numeric(zm[i,2:50])),lwd=2,col=rpois(1,200))}

zikas=function(kk){

time=seq(1:49)
zika_cases=c(as.numeric(zm[kk,2:50]))
micro_cases=(as.numeric(micros[kk,5:53]))

zm1=cbind(time,micro_cases,zika_cases)

ltt=array(0,dim=c(19))
ltq=array(0,dim=c(19))
lte=array(0,dim=c(19))
ltr=array(0,dim=c(19))

zika_nums=array(0,dim=c(30))

wk_range=seq(1,19)

for(j in 1:19){

    res=zm1[20:49,1]-wk_range[j]
    
    for(i in 1:30){
        zika_nums[i]=zm1[zm1[,1]==res[i],3]}
    
    wk_dat=cbind(zika_nums,zm1[20:49,2])
    
    
gaus_fn=function(p){
    nlen=length(wk_dat[,1])
    mod=array(0,dim=c(nlen))
    tmp=0
    
    for(i in 1:nlen){
        mod[i]=p[1]+p[2]*wk_dat[i,1]
        likl=(wk_dat[i,2]-mod[i])^2
        tmp=tmp+likl	}
    (nlen)*log(2*pi*p[3]^2)+tmp/(2*p[3]^2)
}


#define the initial parameter set
p=c(1.0,1.0,1.0)

#use appropriate optimization routine (in R - default is Nelder-Mead)
out_gau=optim(p,gaus_fn)

#displays parameter estimates
ltt[j]=out_gau$value
ltq[j]=out_gau$par[3]
lte[j]=out_gau$par[1]
ltr[j]=out_gau$par[2]
}

rtt=cbind(wk_range,ltt,ltq,lte,ltr)
    rtt[rtt[,2]==min(rtt[,2]),]}

rest=array(0,dim=c(53,5))
for(i in 1:53){
    a=zikas(i)
    for(kj in 1:5){
        rest[i,kj]=as.numeric(a[kj])}
}

res=rest[micros[,55]>0,]

vd=res[1:16,]

mean(20-vd[,1])
sqrt(var(20-vd[,1])/16)*1.96

hist(20-vd[,1],col="grey",xlab="Weeks post conception",ylab="Frequency",main="Zika week predicting MC")

dev.off()