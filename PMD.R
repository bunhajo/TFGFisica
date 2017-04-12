#### libraries ####
require(reshape)
require(ggplot2)
require(gtable)
require(lattice)
require(gridExtra)
require(grid)
require(dplyr)
require(scales)
require(grDevices)
require(treemap)
require(data.tree)
require(stringr)
library(RColorBrewer)  
library(reshape2) 
library(lomb)
#### Style ####
style = theme(#panel.grid.major = element_blank(),
                    #panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.title.x = element_text(size = 13,
                                                hjust=1),
                    axis.title.y = element_text(size = 13,
                                                vjust=1.5),
                    axis.text.x = element_text(size= 8,
                                               angle = 0,
                                               face = 'bold',
                                               vjust = 1,
                                               hjust=0.5),
                    axis.text.y = element_text(size= 8,
                                               face = 'bold',
                                               vjust = 0.5,
                                               hjust=1),
                    plot.title = element_text(hjust = 0),
                    legend.text = element_text(size=15)) +
  theme(plot.title = element_text(size=50)
  )  
#### Lomb Scargle ####

#lomb=lsp(df[,1:2], from = 0, to = 0.6,  ofac = 100)
#saveRDS(lomb, "Lomb-Scargle.rds")


lomb=readRDS("Lomb-Scargle.RDS")

l_plot=data.frame(Freq=lomb$scanned, Power=lomb$power)

ggplot(l_plot, aes(x=Freq, y=Power))+geom_line()+
  xlab(expression(paste("Frequency [", day^{-1}, "]")))+
  ylab(expression(paste("-log"[10], "(Probability)")))+
  scale_x_continuous(breaks=c(seq(0,0.6,0.05)),labels=abbreviate)+
  style


#### Computation ####




#Read data
df=read.delim("C:/Users/Tornero/Documents/Uni/TFGFisica/HESS_LS5039.dat", header = FALSE, sep= "")
names(df)=c("time", "data", "error")
#today=57817 05/03/2017

#Initial Plot
ggplot(df, aes(x=time, y=10^12*data, ymin=10^12*(data-error), ymax=10^12*(data+error)))+
  geom_point()+
  geom_smooth()+
  geom_errorbar()+
  annotate("text", x=c(53200, 53600), y=c(3.5,3.5), label=c("2004", "2005"), size=10)+
  xlab("MJD")+ylab(expression(paste(Delta,"F > 1 TeV [", 10^{-12},"ph. ", cm^{-2}, s^{-1},"]")))+
  theme(axis.title.x=element_text(hjust = 1))

#time conversion
t1=min(df$time)
df$time=df$time-t1
df$time=df$time*86400


#Statistic Definition
# samples_stat=function(df, P=3.9, t1=53148.08){
# 
# 
# Per=86400*P
# phi=(df$time)/Per-trunc((df$time)/Per)
# aux=data.frame(time=phi, data=df$data)
# s=double()
# n=double()
# j=0
# for(i in seq(0,0.9,0.1)){
#   j=j+1
#   s[j]=var(aux$data[which(aux$time >= i & aux$time < i+0.1)])
#   n[j]=length(which(aux$time >= i & aux$time < i+0.1))
# }
# s_sample=sum((n-1)*s)/(sum(n)-10)
# 
# return(s_sample^2/var(df$data)^2)
# 
# }


#Statistic Definition
samples_stat=function(df, P=3.9, N_b=5, N_c=5){
  
  Per=86400*P
  phi=(df$time)/Per-trunc((df$time)/Per)
  aux=data.frame(time=phi, data=df$data)
  s=double()
  n=double()
  j=0
  for (k in 0:(N_c-1)){
    adj= k/(N_c*N_b)
    for(i in seq(-adj,1,1/N_b)){
      points=which(aux$time >= i & aux$time < i+1/N_b)
      if(length(points)>0){
        j=j+1
        s[j]=sum((aux$data[points]-mean(aux$data[points]))^2)/(length(points)-1)
        n[j]=length(points)
        n=n[which(!is.na(s))]
        s=s[which(!is.na(s))]
      }
    }
  }  
  s_sample=sum((n-1)*s, na.rm = T)/(sum(n, na.rm = T)-N_b*N_c)
  var_data=sum((aux$data-mean(aux$data))^2)/(nrow(aux)-1)
  
  return(s_sample/var_data)
}
#Statistic Calculation
result=double()
periods_prova=seq(3,5,0.001)
j=0
# for(N_c in 1:3)
for(i in periods_prova){
  j=j+1
  result[j]= samples_stat(df, P=i,N_b=25, N_c = 2)
}

#Statistic Results
stat_results=data.frame(period=periods_prova, theta=result)
minim_theta=min(na.omit(stat_results$theta))
minim_period=stat_results$period[which(stat_results$theta==minim_theta)]



### Results
copy_10.1 
copy_10.10  
copy_5.5 
copy_25.2 = stat_results # Exact Period
copy_50.1 # Exact Period  

full_results = full_join(copy_10.1, full_join(copy_10.10, full_join(copy_5.5, full_join(copy_25.2, copy_50.1))))
full_results$n_result=rownames(full_results)
full_results$n_result=as.double(trunc(as.double(full_results$n_result)/(length(periods_prova)+1)))
minim_periods=results %>% group_by(n_result) %>% summarise(minim_period= period[which(theta==min(theta))])
minim_theta=full_results %>% group_by(n_result) %>% summarise(minim_theta= min(theta))
saveRDS(full_results, "Results.rds")


#Period Determination


results=readRDS("Results.rds")

#b_period = mean(minim_periods$minim_period)
b_period=3.9082
sqrt(var(minim_periods$minim_period))

Per=86400*b_period
#Per=86400*3.9063
phi=(df$time)/Per-trunc((df$time)/Per)
b_per=data.frame(time=phi, data=df$data, error=df$error)
saveRDS(b_per, "best_period.rds")
#### Plots ####

# Plot Stat
p=ggplot(results[which(results$n_result==0),], aes(x=period, y=theta, colour=factor(n_result)))+
  geom_line()+
  ggtitle(expression(paste("N",scriptscriptstyle(b),"=10","   N",scriptscriptstyle(c),"=1")))+
#  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+scale_color_discrete(h=c(90,270), guide=FALSE)+
  geom_segment(aes(x=3.8, xend=4, y=0.4, yend=0.4),colour="red")+
  geom_segment(aes(x=3.8, xend=4, y=0.85, yend=0.85),colour="red")+
  geom_segment(aes(x=3.8, xend=3.8, y=0.4, yend=0.85),colour="red")+
  geom_segment(aes(x=4, xend=4, y=0.4, yend=0.85),colour="red")
  #xlim(c(3.9,3.92))+
  #scale_color_discrete(high = "green", low="blue")

q=ggplot(results[which(results$n_result==1),], aes(x=period, y=theta, colour=factor(n_result)))+
  geom_line()+
  ggtitle(expression(paste("N",scriptscriptstyle(b),"=10","   N",scriptscriptstyle(c),"=10")))+
  #  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+scale_color_discrete(h=c(180,270), guide=FALSE)+
  geom_segment(aes(x=3.8, xend=4, y=0.4, yend=0.4),colour="red")+
  geom_segment(aes(x=3.8, xend=4, y=0.85, yend=0.85),colour="red")+
  geom_segment(aes(x=3.8, xend=3.8, y=0.4, yend=0.85),colour="red")+
  geom_segment(aes(x=4, xend=4, y=0.4, yend=0.85),colour="red")

r=ggplot(results[which(results$n_result==2),], aes(x=period, y=theta, colour=factor(n_result)))+
  geom_line()+
  #  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+scale_color_discrete(h=c(90,180), guide=FALSE)+
  geom_segment(aes(x=3.8, xend=4, y=0.4, yend=0.4),colour="red")+
  geom_segment(aes(x=3.8, xend=4, y=0.85, yend=0.85),colour="red")+
  geom_segment(aes(x=3.8, xend=3.8, y=0.4, yend=0.85),colour="red")+
  geom_segment(aes(x=4, xend=4, y=0.4, yend=0.85),colour="red")

s=ggplot(results[which(results$n_result==3),], aes(x=period, y=theta, colour=factor(n_result)))+
  geom_line()+
  ggtitle(expression(paste("N",scriptscriptstyle(b),"=25","   N",scriptscriptstyle(c),"=2")))+
  #  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+scale_color_discrete(h=c(90,180), guide=FALSE)+
  geom_segment(aes(x=3.8, xend=4, y=0.3, yend=0.3),colour="red")+
  geom_segment(aes(x=3.8, xend=4, y=0.85, yend=0.85),colour="red")+
  geom_segment(aes(x=3.8, xend=3.8, y=0.3, yend=0.85),colour="red")+
  geom_segment(aes(x=4, xend=4, y=0.3, yend=0.85),colour="red")

t=ggplot(results[which(results$n_result==4),], aes(x=period, y=theta, colour=factor(n_result)))+
  geom_line()+
  ggtitle(expression(paste("N",scriptscriptstyle(b),"=50","   N",scriptscriptstyle(c),"=1")))+
  #  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+scale_color_discrete(h=c(180,270), guide=FALSE)+
  geom_segment(aes(x=3.8, xend=4, y=0.3, yend=0.3),colour="red")+
  geom_segment(aes(x=3.8, xend=4, y=0.85, yend=0.85),colour="red")+
  geom_segment(aes(x=3.8, xend=3.8, y=0.3, yend=0.85),colour="red")+
  geom_segment(aes(x=4, xend=4, y=0.3, yend=0.85),colour="red")

multiplot(plotlist = list(p,q,s,t), cols=2)

#Plot Best Period
ggplot(b_per, aes(x=100*time, y=data/mean(data), ymin=(data-error)/mean(data), ymax=(data+error)/mean(data)))+
  geom_point()+
  geom_smooth(se=var(b_per$data))+
  geom_errorbar()+
  xlab("Phase")+
  ylab(expression(paste(Delta,"F > 1 TeV [", 10^{-12},"ph. ", cm^{-2}, s^{-1},"]")))


