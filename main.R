source("PMD.R")
library(dplyr)

parameter = data.frame(N_b=c(5,5,5,10,10,25,50), N_c=c(2,4,5,2,5,2,1))
periods_prova=seq(1,6,0.001)
k=0
# for(N_c in 1:3)

for(i in 1:nrow(parameter)){
  j=0  
  for(j in periods_prova){
    result[k]= samples_stat(df, P=j,N_b=parameter$N_b[i], N_c = parameter$N_c[i])
    k=k+1
  }
}

aux=as.double(unlist(result))
result=data.frame(period= rep(periods_prova,7), theta=aux, n_result=c(1:length(aux)))
result$n_result=as.double(trunc(as.double(result$n_result)/(length(periods_prova)+1)))
minim_periods=result %>% group_by(n_result) %>% summarise(minim_period= period[which(theta==min(theta))])
minim_theta=result %>% group_by(n_result) %>% summarise(minim_theta= min(theta))
saveRDS(result, "Results.rds")

p=list()
j=0
for(i in c(0,2,4,5)){
j=j+1
  p[[j]]=ggplot(result[which(result$n_result==i),], aes(x=period, y=theta, colour=factor(n_result)))+
  geom_line()+
  #ggtitle(expression(paste("N",scriptscriptstyle(b),"=10","   N",scriptscriptstyle(c),"=10")))+
    
  #  annotate("text", x=c(minim_period, 1.5*minim_period), y=c(0,0), label=c("Period", "Period and half"))+
  ylab(expression(theta))+scale_color_discrete(h=c(180,270), guide=FALSE)+
  geom_segment(aes(x=3.8, xend=4, y=0.4, yend=0.4),colour="red")+
  geom_segment(aes(x=3.8, xend=4, y=0.85, yend=0.85),colour="red")+
  geom_segment(aes(x=3.8, xend=3.8, y=0.4, yend=0.85),colour="red")+
  geom_segment(aes(x=4, xend=4, y=0.4, yend=0.85),colour="red")
}

i=4
p[[i]]=p[[i]]+ggtitle(expression(paste("N",scriptscriptstyle(b),"=25","   N",scriptscriptstyle(c),"=2")))
saveRDS(p, "plots_theta.rds")

multiplot(plotlist = p, cols = 2)


