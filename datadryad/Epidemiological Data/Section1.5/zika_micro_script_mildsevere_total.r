#new total zika-microcephaly analysis across all locations

rm(list=ls())
setwd("~/desktop/")

zm=read.table("milsev_mc_zika.txt",header=T)

pdf("mildsevere_zika_micro_epid.pdf")
plot(zm[,1],sqrt(zm[,2]),type="l",xlab="Time (Weeks)",ylab="Sqrt(Cases)",bty="l",lwd=2,col="purple",main="Zika - Microcephaly cases")
lines(zm[,1],sqrt(zm[,3]),lwd=2)
dev.off()

ltt=array(0,dim=c(15))
ltq=array(0,dim=c(15))
lte=array(0,dim=c(15))
ltr=array(0,dim=c(15))

zika_nums=array(0,dim=c(20))

wk_range=seq(1,15)

for(j in 1:15){

res=zm[23:42,1]-wk_range[j]

for(i in 1:20){
    zika_nums[i]=zm[zm[,1]==res[i],2]}

    wk_dat=cbind(zika_nums,zm[23:42,3])

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
rtt[rtt[,2]==min(rtt[,2]),]