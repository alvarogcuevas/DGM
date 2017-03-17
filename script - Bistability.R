library(rootSolve)
library(ggplot2)
library(reshape2)
library(pushoverr)

source('~/Concordia/Research/Git/Djehiche_Giesecke/functions - Copy.R')

set.seed(5678)
x0<-0.75
varepsilon_vec<-2e-03#c(0.004)
Omega_a<-0.5
Omega_b<-1
mark_set<-c(-1,1)
T_vec<-2500
low_t<-0
N_steps<-62500
h<-(T_vec-low_t)/N_steps
#sample_size<-207
sys_cap<-100
gamma_tilde<-1

#####################################################
# Program #
###########

initial<-data.frame(Z=x0,s=low_t,A=0,E_n=rexp(n = 1,rate = 1),i = 0,L = 1)
time<-low_t
big_chain_Z<-big_chain_s<-indicators<-c()

while(time<T_vec){
  small_chain<-chain(Z = initial$Z,s = initial$s,A = initial$A,E_n = initial$E_n,i = initial$i,L = initial$L,T_ = T_vec,t0 = low_t,epsilon = varepsilon_vec,Nsteps = N_steps,c_ = 0,new = F)
  big_chain_Z<-c(big_chain_Z,small_chain$Z[-length(small_chain$Z)])
  big_chain_s<-c(big_chain_s,small_chain$s[-length(small_chain$s)])
  indicators<-c(indicators,small_chain$s[length(small_chain$s)]-small_chain$s[1])
  initial<-data.frame(Z =  small_chain$Z[length(small_chain$Z)],s =  small_chain$s[length(small_chain$s)],A =  small_chain$A[length(small_chain$A)],E_n =  small_chain$E_n[length(small_chain$E_n)],i =  small_chain$i[length( small_chain$i)],L =  small_chain$L[length( small_chain$L)])
  Omega_a<-0.5*(small_chain$Z[length(small_chain$Z)]>0.5)
  Omega_b<-0.5+0.5*(small_chain$Z[length(small_chain$Z)]>0.5)
  time<-small_chain$s[length(small_chain$s)]
  print(c(time,"//",T_vec))
}

old_indicators<-indicators
old_big_chain<-big_chain

big_chain<-data.frame(s=big_chain_s,Z=big_chain_Z)
p<-ggplot(big_chain,aes(x = s,y = Z))
p+geom_line()+ ylab("Z(s)")
ggsave("original.png")


min(which(old_big_chain$Z < 0.5))
index<-match(old_indicators[1],big_chain$s)
first_chain<-data.frame(s=big_chain$s[1:index],Z=big_chain$Z[1:index])
q<-ggplot(first_chain,aes(x = s,y = Z))
q+geom_line()+ ylab("Z(s)")
ggsave("original_detail.png")

r<-ggplot(big_chain, aes(Z))
r+geom_density()+xlab("Z(s)")
ggsave("original_hist.png")

indicators.df<-data.frame(indic=indicators)
s<-ggplot(indicators.df,aes(indic))
s+geom_histogram()+xlab("Histogram of frequencies of time between states")

pushover(message = "Hieeee")



###############################################################################################################  
# Intermediate #
################
# 
# c_vec<-rep(0,length(parameters[,2]))
# for(k in 1:length(parameters[,2])){
#   c_vec[k]<-optimize(min_max,interval =c(c_H(),100000),maximum =TRUE ,cap_t=parameters[k,2], low_t=low_t,a=Omega_a,b_=Omega_b)$maximum
# }
#
###############################################################################################################  
# After that #
##############

initial<-data.frame(Z=x0,s=low_t,A=0,E_n=rexp(n = 1,rate = 1),i = 0,L = 1)
time<-low_t
big_chain_Z<-big_chain_s<-indicators<-c()

while(time<T_vec){
  c_<-optimize(min_max,interval=c(c_H(),1000),maximum =TRUE ,cap_t=T_vec,low_t=initial$s,a=Omega_a,b_=Omega_b,x_0=initial$Z)$maximum
  small_chain<-chain(Z = initial$Z,s = initial$s,A = initial$A,E_n = initial$E_n,i = initial$i,L = initial$L,T_ = T_vec,t0 = low_t,epsilon = varepsilon_vec,Nsteps = N_steps,c_ = c_,new = T)
  big_chain_Z<-c(big_chain_Z,small_chain$Z[-length(small_chain$Z)])
  big_chain_s<-c(big_chain_s,small_chain$s[-length(small_chain$s)])
  indicators<-c(indicators,small_chain$s[length(small_chain$s)]-small_chain$s[1])
  initial<-data.frame(Z =  small_chain$Z[length(small_chain$Z)],s =  small_chain$s[length(small_chain$s)],A =  small_chain$A[length(small_chain$A)],E_n =  small_chain$E_n[length(small_chain$E_n)],i =  small_chain$i[length( small_chain$i)],L =  small_chain$L[length( small_chain$L)])
  Omega_a<-0.5*(small_chain$Z[length(small_chain$Z)]>0.5)
  Omega_b<-0.5+0.5*(small_chain$Z[length(small_chain$Z)]>0.5)
  time<-small_chain$s[length(small_chain$s)]
  print(c(time,"//",T_vec))
}

big_chain<-data.frame(s=big_chain_s,Z=big_chain_Z)
p2<-ggplot(big_chain,aes(x = s,y = Z))
p2+geom_line()

big_chain2<-data.frame(s=big_chain_s,Z=big_chain_Z)
p2<-ggplot(big_chain2,aes(x = s,y = Z))
p2+geom_line()+ ylab("Z(s)")
ggsave("new.png")

index2.<-min(which(big_chain$Z < 0.5))
index2<-match(indicators[1],big_chain2$s)
first_chain2<-data.frame(s=big_chain2$s[1:index2.],Z=big_chain2$Z[1:index2.])
q2<-ggplot(first_chain2,aes(x = s,y = Z))
q2+geom_line()+ ylab("Z(s)")
ggsave("new_detail.png")

r2<-ggplot(big_chain2, aes(Z))
r2+geom_density()+xlab("Z(s)")
ggsave("new_hist.png")

indicators2.df<-data.frame(indic=indicators)
s2<-ggplot(indicators2.df,aes(indic))
s2+geom_histogram()+xlab("Histogram of frequencies of time between states")

pushover(message = "Hieeee")