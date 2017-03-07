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
metri_style = theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_blank(),
                    axis.title.x = element_text(size = 26,
                                                vjust=0.1),
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
ggplot(df, aes(x=time, y=data))+
  geom_point()+
  geom_smooth()+
  geom_errorbar(data=df, aes(x=time, ymin=data-error, ymax=data+error))+
  metri_style

#time conversion
t1=min(df$time)
df$time=df$time-t1
df$time=df$time*86400

#Statistic Definition
samples_stat=function(df, P=3.9){

Per=86400*P
phi=(df$time)/Per-trunc((df$time)/Per)
aux=data.frame(time=phi, data=df$data)
s=double()
n=double()
j=0
for(i in seq(0,0.9,0.1)){
  j=j+1
  s[j]=var(aux$data[which(aux$time >= i & aux$time < i+0.1)])
  n[j]=length(which(aux$time >= i & aux$time < i+0.1))
}
s_sample=sum((n-1)*s)/(sum(n)-10)

return(s_sample^2/var(df$data)^2)

}

#Statistic Calculation
result=double()
periods_prova=seq(3,8,0.01)
j=0
for(i in periods_prova){
  j=j+1
  result[j]= samples_stat(df, P=i)
}


#Statistic Results
stat_results=data.frame(period=periods_prova, theta=result)
minim_period=stat_results$period[which(stat_results$theta==min(na.omit(stat_results$theta)))]
minim_stat=min(na.omit(stat_results$theta))

#Best period 
Per=86400*minim_period
phi=(df$time)/Per-trunc((df$time)/Per)
b_per=data.frame(time=phi, data=df$data, error=df$error)
#### Plots ####

# Plot Stat
ggplot(stat_results, aes(x=period, y=theta))+
  geom_line()+
#  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+metri_style
  
#Plot Best Period
ggplot(b_per, aes(x=100*time, y=data/mean(data), ymin=(data-error)/mean(data), ymax=(data+error)/mean(data)))+
  geom_point()+
  geom_smooth(se=var(b_per$data))+
  geom_errorbar()+
  xlab("Phase")+ylab("Energy")+metri_style


