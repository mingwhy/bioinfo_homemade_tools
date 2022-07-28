library(ggplot2);library(gridExtra)
#BiocManager::install('ReIns')
library(ReIns)
# Plot of the PDF
x <- seq(0, 10, 0.01)
plot(x, dtexp(x, rate = 1, endpoint=10), xlab="x", ylab="PDF", type="l",col='darkred',lwd=2)
# Plot of the CDF
x <- seq(0, 10, 0.01)
plot(x, ptexp(x, rate = 2, endpoint=5), xlab="x", ylab="CDF", type="l")

###########################################################
## closed-form of mean and variance compared with simulation
theo_mean<-function(u,t){
  1/u-t/(exp(u*t)-1)
}
theo_var<-function(u,t){
  -(t^2*u^2+2*t*u)/(u^2*(exp(t*u)-1)) - (1/u-t/(exp(t*u)-1))^2 + 2/u^2
}

u=1/c(0.1,1,10,100,1000); #rate, xx per day
right.b=c(1,10,100,1000,10000); #days
df=expand.grid(u,right.b)
colnames(df)=c('rate','right.truncted')
out=t(sapply(1:nrow(df),function(i){
  c(theo_mean(df[i,1],df[i,2]),
    theo_var(df[i,1],df[i,2]))
}))
df$analytical_mean=out[,1];
df$analytical_var=out[,2];

n=1e5; #sample
out2=t(sapply(1:nrow(df),function(i){
  tmp=rtexp(n, rate = df[i,1], endpoint=df[i,2])
  c(mean(tmp),var(tmp))
}))
df$simu_mean=out2[,1];
df$simu_var=out2[,2];

df$rate=factor(df$rate)
p1<-ggplot(df,aes(x=right.truncted,y=analytical_mean,group=rate,col=rate))+
  geom_line()+xlab('right trunctated time point (day)')+ylab('analytical or simulated mean')+
  geom_point(aes(x=right.truncted,y=simu_mean))+theme_bw(base_size = 20)
p1+scale_y_log10()+scale_x_log10()

p2<-ggplot(df,aes(x=right.truncted,y=analytical_var,group=rate,col=rate))+
  geom_line()+xlab('right trunctated time point (day)')+ylab('analytical or simulated var')+
  geom_point(aes(x=right.truncted,y=simu_var))+theme_bw(base_size = 20)
p2+scale_y_log10()+scale_x_log10()

grid.arrange(p1+scale_y_log10()+scale_x_log10(),p2+scale_y_log10()+scale_x_log10(),
             ncol=2)
###################################################
age_slope<-function(k,u,t){
  tmp=exp(t*u)
  dvar=t*tmp * (t*u+tmp*(t*u-2)+2) / (tmp-1)^3 #https://www.wolframalpha.com/input?i=D%5BDivide%5B2%2CPower%5Bu%2C2%5D%5D-Divide%5B%5C%2840%29Power%5Bt%2C2%5D*Power%5Bu%2C2%5D%2B2*t*u%5C%2841%29%2CPower%5Bu%2C2%5D%5C%2840%29Power%5Be%2Cu*t%5D-1%5C%2841%29%5D+-+Power%5B%5C%2840%29Divide%5B1%2Cu%5D-Divide%5Bt%2CPower%5Be%2Cu*t%5D-1%5D%5C%2841%29%2C2%5D%2Ct%5D&assumption=%22UnitClash%22+-%3E+%7B%22u%22%2C+%7B%22AtomicMassUnit%22%7D%7D&assumption=%7B%22C%22%2C+%22u%22%7D+-%3E+%7B%22Variable%22%7D&assumption=%7B%22C%22%2C+%22u%22%7D+-%3E+%7B%22Variable%22%7D&assumption=%7B%22MC%22%2C+%22u*t%22%7D+-%3E+%7B%22Variable%22%7D
  dmean= (tmp*(t*u-1) +1) / (tmp-1)^2 #https://www.wolframalpha.com/input?i2d=true&i=D%5B%5C%2840%29Divide%5B1%2Cu%5D-Divide%5Bt%2C%5C%2840%29Power%5Be%2Cu*t%5D-1%5C%2841%29%5D%5C%2841%29%2Ct%5D
  k*dmean + k^2*dvar
}


u=1/c(0.1,1,10,100,1000); #rate, xx per day
right.b=c(1,10,100,1000,10000); #days
#k=0.01;
plots=list();
for(k in c(0.001,0.01)){
  df=expand.grid(u,right.b)
  colnames(df)=c('rate','right.truncted')
  out=lapply(1:nrow(df),function(i){
    age_slope(k,df[i,1],df[i,2])
  })
  df$slope=unlist(out)
  
  df$rate=factor(df$rate) 
  #when exp(tu) becomes too big, r shows Inf
  p1<-ggplot(df,aes(x=right.truncted,y=slope,group=rate,col=rate))+
    geom_line()+xlab('right trunctated time point (day)')+ylab('d(Var(Xt))/dt')+
    theme_bw(base_size = 20)+ggtitle(paste0('k=',k))
  plots[[as.character(k)]]<-p1+scale_x_log10()+
    scale_y_continuous(trans = scales::pseudo_log_trans(1e-20, 10),
                       breaks=c(0, 1e-20,1e-10,0.01, 0.1, 1))
}
grid.arrange(grobs=plots,ncol=2)

###################################################
# time: days
x <- seq(0, 30*30, 1)
lambda=1/c(0.1,0.3,1,3,10,30,90,900,3000,10000,1e6)
library(RColorBrewer)
mycols=RColorBrewer::brewer.pal(length(lambda),'RdYlBu')
b=c(3,6,9,12,18,24)*30 #month->days
par(mfrow=c(2,3))
for(bt in b){
  out=lapply(lambda,function(i){
    ptexp(x, rate = i, endpoint=bt)
  })
  sapply(out,summary)
  plot(0,0, xlab="x", ylab="PDF",main=bt,
       xlim=c(0,30*30),ylim=c(0,1),type='n')
  for(i in 1:length(out)){
    points(x,out[[i]], type="l",col=mycols[i])
  }
}

## generate random number and calcualte mean, var
library(ggplot2);library(RColorBrewer)
mycols=RColorBrewer::brewer.pal(length(lambda),'RdYlBu')

lambda=1/c(0.1,0.3,1,3,10,30,90,300,900,3000,10000)
b=c(3,6,9,12,18,24)*30 #month->days
n=10000;
out=lapply(b,function(bt){
  out=lapply(lambda,function(i){
    tmp=rtexp(n, rate = i, endpoint=bt)
    #hist(tmp)
    tmp
  })
  meanx=sapply(out,mean)
  varx=sapply(out,var)
  df=data.frame(trunc_age=bt,rate=lambda,mean=meanx,var=varx)
  return(df)
})

df.out=as.data.frame(Reduce(`rbind`,out))
head(df.out)

df.out$lifespan=1/df.out$rate
df.out$lifespan=factor(df.out$lifespan)
p1<-ggplot(df.out,aes(x=trunc_age,y=mean,col=lifespan,group=lifespan))+
  geom_point()+geom_line()+theme_bw()
p1

p2<-ggplot(df.out,aes(x=trunc_age,y=var,col=lifespan,group=lifespan))+
  geom_point()+geom_line()+theme_bw()
p2
