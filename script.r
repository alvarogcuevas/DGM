library(rootSolve)
source('C:/Users/User/Google Drive/Research/Git/Djehiche_Giesecke/functions.R')

#library(pushoverr)

set.seed(1234)
x0<-100
lambda_0<-5
lambda_1<-50
varepsilon_vec<-c(0.05,0.06,0.07,0.08)
Omega_a<-90
Omega_b<-110
mark_set<-c(0.01, 0.015, 0.02, 0.025, 0.03)
T_vec<-0.04
parameters<-expand.grid(varepsilon_vec,T_vec)
low_t<-0
N_steps<-100
sample_size<-10^3

#####################################################
# Program #
###########

for(k in 1:length(parameters[,1])){
  time_trials<-indicators<-prom<-c()
  for(j in 1:sample_size){
    ZandL<-chain(x0 = x0,t0 = low_t,T_ = parameters[k,2],epsilon = parameters[k,1],Nsteps = N_steps,c_=0,new = FALSE)
    indicators<-c(indicators,1*(ZandL$s[length(ZandL$s)]<parameters[k,2]))
    print(c(k,j))
  }
  estim<-c(mean(indicators),(1/sample_size)*sqrt((1/mean(indicators))-1))
  write.table(x=c(date(),"## BEFORE ##",sample_size,N_steps,parameters[k,],estim),file=paste("C:/Users/User/Google Drive/Research/Code/[Djehiche+Giesecke] Epsilon Jump Diffusion v",vers,".txt"), append=TRUE,row.names = F, col.names = F, quote=F, sep="//")
}

# time_used.old<-proc.time()-start_time.old

###############################################################################################################  
# Intermediate #
################

c_vec<-rep(0,length(parameters[,2]))
for(k in 1:length(parameters[,2])){
  c_vec[k]<-optimize(min_max,interval =c(c_H(),100000),maximum =TRUE ,cap_t=parameters[k,2], low_t=low_t,a=Omega_a,b=Omega_b)$maximum
}

###############################################################################################################  
# After that #
##############

for(k in 1:length(parameters[,1])){
  time_trials<-indicators<-prom<-c()
  #c_<-c_vec[k]
  for(j in 1:sample_size){
    ZandL<-chain(x0 = x0,t0 = low_t,T_ = parameters[k,2],epsilon = parameters[k,1],Nsteps = N_steps,c_=c_vec[k],new = TRUE)    #prom<-c(prom,exp(-(1/parameters[k,1])*g_function(ZandL$Z[length(ZandL$Z)])))
    indicators<-c(indicators,1*(ZandL$s[length(ZandL$s)]<parameters[k,2])*(ZandL$L[length(ZandL$L)])^-1)
    print(c(k,j))
  }
  estim<-c(mean(indicators),(1/sample_size)*sqrt((1/mean(indicators))-1))
  write.table(x=c(date(),"## AFTER ##",sample_size,N_steps,parameters[k,],estim),file=paste("C:/Users/User/Google Drive/Research/Code/[Djehiche+Giesecke] Epsilon Jump Diffusion v",vers,".txt"), append=TRUE,row.names = F, col.names = F, quote=F, sep="//")
}
#pushover(message="Your shit is fresh")
