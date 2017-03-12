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
#### Style ####
style = theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_blank(),
                    axis.title.x = element_text(size = 26,
                                                vjust=0),
                    axis.title.y = element_text(size = 26,
                                                vjust=1.5),
                    axis.text.x = element_text(size= 16,
                                               angle = 45,
                                               face = 'bold',
                                               vjust = 1,
                                               hjust=1),
                    axis.text.y = element_text(size= 16,
                                               face = 'bold',
                                               vjust = 0.5,
                                               hjust=1),
                    plot.title = element_text(hjust = 0),
                    legend.text = element_text(size=15)) +
  theme(plot.title = element_text(size=50)
  )  
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
samples_stat=function(df, P=3.9, N_b=10, N_c=1){
  
  Per=86400*P
  phi=(df$time)/Per-trunc((df$time)/Per)
  aux=data.frame(time=phi, data=df$data)
  s=double()
  n=double()
  j=0
  for (i in 0:(N_c-1)){
    adj=i/(N_c*N_b)
    for(i in seq(adj,0.9+adj,1/N_b)){
      j=j+1
      s[j]=var(aux$data[which(aux$time >= i & aux$time < i+1/N_b)])
      n[j]=length(which(aux$time >= i & aux$time < i+1/N_b))
    }
  }  
  s_sample=sum((n-1)*s, na.rm = T)/(sum(n, na.rm = T)-N_b*N_c)
  
  return(s_sample^2/var(df$data)^2)
  
}
#Statistic Calculation
result=double()
periods_prova=seq(3,5,0.001)
j=0
for(N_c in 1:3)
for(i in periods_prova){
  j=j+1
  result[j]= samples_stat(df, P=i,N_c = N_c)
}


#Statistic Results
stat_results=data.frame(period=rep(periods_prova,3), theta=result, N_c=c(rep(1,length(periods_prova)),rep(2,length(periods_prova)),rep(3,length(periods_prova))))
minim_period[1]=stat_results$period[which(stat_results$theta==min(na.omit(stat_results$theta[which(stat_results$N_c==1)])))]
minim_period[2]=stat_results$period[which(stat_results$theta==min(na.omit(stat_results$theta[which(stat_results$N_c==2)])))]
minim_period[3]=stat_results$period[which(stat_results$theta==min(na.omit(stat_results$theta[which(stat_results$N_c==3)])))]

minim_stat=min(na.omit(stat_results$theta))

#Best period 
Per=86400*minim_period[1]
phi=(df$time)/Per-trunc((df$time)/Per)
b_per=data.frame(time=phi, data=df$data, error=df$error)
#### Plots ####

# Plot Stat
ggplot(stat_results, aes(x=period, y=theta, colour=N_c))+
  geom_line()+
#  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+
  scale_color_gradient(high = "red", low="blue", guide = FALSE)
  
#Plot Best Period
ggplot(b_per, aes(x=100*time, y=data/mean(data), ymin=(data-error)/mean(data), ymax=(data+error)/mean(data)))+
  geom_point()+
  geom_smooth(se=var(b_per$data))+
  geom_errorbar()+
  xlab("Phase")+
  ylab(expression(paste(Delta,"F > 1 TeV [", 10^{-12},"ph. ", cm^{-2}, s^{-1},"]")))


