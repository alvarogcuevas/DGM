library(rootSolve)
library(ggplot2)
library(reshape2)
library(pushoverr)

source('~/Concordia/Research/Git/DGM/functions - Bistability - April.R')

set.seed(5678)

parameters<-list(x0=0.25,low_t=0,cap_t=2500, omega_a=0, omega_b=0.5,mark_set=c(-1,1),varepsilon=0.02, N=1500, c_=0, h=0.1, new=F,stop.at.tau=T)

b<-function(x){
  -x + (1- x)-(16/3)*(1-x)*x+(32/3)*x^2*(1-x)
}

sigma_fn<-function(x){
  sqrt(x*(1-x))
}

lambda<-function(x){
  0.00001#(1/2)*varepsilon*x*(1-x)*N^2
}

gamma_fn<-function(y){
  0*y
}

random_mark<-function(){
  sample(x =mark_set,replace = FALSE,size = 1)
}

giant_vector<-c()
prueba<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn)
prueba$s[length(prueba$s)]
c_<-optimize(min_max,interval=c(c_H(),1000),maximum =TRUE, cap_t=parameters$cap_t,low_t=parameters$low_t,omega_a=parameters$omega_a,omega_b=parameters$omega_b,x0=parameters$x0)$maximum
parameters$c_<-c_
prueba.new<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn)
prueba.new$s[length(prueba.new$s)]

results_matrix<-sapply(X=1:3,FUN = function(iterNumber){
  one_path<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn)$Z
  giant_vector<<-c(giant_vector,one_path)
  one_path
})

try1<-melt(results_matrix[,1])
try2<-melt(results_matrix[,2])
try3<-melt(results_matrix[,3])
smmry<-melt(results_matrix)

r<-ggplot(try1, aes(value))
r+geom_density()+xlab("Z(s)")

r<-ggplot(try2, aes(value))
r+geom_density()+xlab("Z(s)")

r<-ggplot(try3, aes(value))
r+geom_density()+xlab("Z(s)")

r<-ggplot(smmry, aes(value))
r+geom_density()+xlab("Z(s)")

#With h=0.01 it takes 6 hours to do one chain

p<-ggplot(prueba,aes(x = s,y = Z))
p+geom_line()+ ylab("Z(s)")+geom_hline(yintercept = 0.25,colour = "red")+geom_hline(yintercept = 0.75,colour = "red")
#ggsave("benchmark2_ts.png")

r<-ggplot(prueba, aes(Z))
r+geom_density()+xlab("Z(s)")
#ggsave("benchmark2_hist.png")

#write.csv(x = prueba, file = "benchmark2.csv")

pushover(message = "HIEEEE, benchmark2 is done")


##########################################################
#Trying the new independent chain function which is able to stop the simulation 
#once it exits the domain
##########################################################

prueba<-chain.independent(Z = x0,s = low_t,T_ = 100,epsilon = varepsilon_vec,stop.at.tau = FALSE)
min(prueba$Z)

prueba.stop<-chain.independent(Z = x0,s = low_t,T_ = 100,epsilon = varepsilon_vec,stop.at.tau = TRUE)
min(prueba.stop$Z)
prueba.stop$s[length(prueba.stop$Z)]

N_sample<-1000

tau<-c()
for(i in 1:N_sample){
  stopped.chain<-chain.independent(Z=x0,s = low_t,T_ = T_vec,epsilon = varepsilon_vec,stop.at.tau=TRUE)
  tau<-c(tau,stopped.chain$s[length(stopped.chain$s)])
}


benchmarks<-as.table(summary(tau.F1.25))+1
#benchmarks<-c(1.5,2,2.5,3)
p_tau<-rep(0,length(benchmarks))
for (i in 1:length(benchmarks)){
  p_tau[i]<-mean(1*(tau.F1.25<benchmarks[i]))
}

relative.error<-sqrt(((1/p_tau)-1)/N_sample)

###############################################################################################################  
# Intermediate #
################

c_<-optimize(min_max,interval=c(c_H(),1000),maximum =TRUE ,cap_t=T_vec,low_t=0,a=Omega_a,b_=Omega_b,x_0=x0)$maximum


###############################################################################################################  
# After that #
##############

tau.new<-L.final<-c()
for(i in 1:N_sample){
  stopped.chain.new<-chain.independent(Z=x0,s = low_t,L=1,T_ = T_vec,epsilon = varepsilon_vec,stop.at.tau=TRUE,new = TRUE)
  tau.new<-c(tau.new,stopped.chain.new$s[length(stopped.chain.new$s)])
  L.final<-c(L.final,stopped.chain.new$L[length(stopped.chain.new$s)])
}

# benchmarks<-as.table(summary(tau.new))+1
p_tau.new<-rep(0,length(benchmarks))
for (i in 1:length(benchmarks)){
  p_tau.new[i]<-mean(1*(tau.F1.25_new[1,]<benchmarks[i])*tau.F1.25_new[2,])
}

relative.error.new<-sqrt(((1/p_tau.new)-1)/N_sample)

curve<-seq(from=1,to=300,by = 0.1)
p.curve<-c()
for (i in 1:length(curve)){
  p.curve[i]<-mean(1*(tau.F1.25_new[1,]<curve[i])*tau.F1.25_new[2,])
}

fittingcurve<-data.frame(x=curve,y=log(1-p.curve))
plot(x=curve,y=log(1-p.curve))
smmry<-lm(y~x,data = fittingcurve)
summary(smmry)
par(mfrow = c(2, 2))
plot(smmry)



smmry.res = resid(smmry) 
qqnorm(smmry.res)
qqline(smmry.res,col=2)
shapiro.test(x = smmry.res)

results<-rbind(p_tau.new,relative.error,relative.error.new)
colnames(x = results)<-benchmarks
save(results, file="results.RData")
