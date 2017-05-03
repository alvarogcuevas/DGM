library(rootSolve)
library(ggplot2)
library(reshape2)
library(pushoverr)

source('~/Concordia/Research/Git/DGM/functions - Bistability - April.R')

set.seed(5678)

parameters<-list(x0=0.25,low_t=0,cap_t=2500, omega_a=0.5, omega_b=1,mark_set=c(-1,1),varepsilon=0.03, N=1500, c_=0, h=0.1, new=F,stop.at.tau=F)

b<-function(x){
  -x + (1- x)-(16/3)*(1-x)*x+(32/3)*x^2*(1-x)
}

sigma_fn<-function(x){
  sqrt(x*(1-x))
}

lambda<-function(x){
  cap_lambda<-(1/2)*varepsilon*x*(1-x)*N^2
  r_plus<-1-x+(32/3)*x^2*(1-x)
  r_minus<-x+(16/3)*x*(1-x)
  cap_lambda+r_plus+r_minus
}

gamma_fn<-function(y){
  1*y
}

random_mark<-function(x){
  cap_lambda<-(1/2)*varepsilon*x*(1-x)*N^2
  r_plus<-1-x+(32/3)*x^2*(1-x)
  r_minus<-x+(16/3)*x*(1-x)
  sample(x =mark_set,replace = FALSE,size = 1,prob = c(r_minus+0.5*cap_lambda,r_plus+0.5*cap_lambda))
}

giant_vector<-c()

Sys.time()
beg<-proc.time()
sample_path<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)
end<-proc.time()-beg

p<-ggplot(sample_path,aes(x = s,y = Z))
p+geom_line()+ ylab("Z(s)")+geom_hline(yintercept = 0.25,colour = "red")+geom_hline(yintercept = 0.75,colour = "red")

write.csv(x=sample_path,file="sample_path.csv")

results_matrix<-sapply(X=1:100,FUN = function(iterNumber){
  one_path<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)$Z
  giant_vector<<-c(giant_vector,one_path)
  one_path
})

smmry<-melt(results_matrix)
r<-ggplot(smmry, aes(value))
r+geom_density()+xlab("Z(s)")

pushover(message = "Main occupation density and sample path is done")

write.csv(x = sample_path, file = "sample_path.csv")
write.csv(x=smmry, file="results_matrix.csv")

############################################################

parameters<-list(x0=0.25,low_t=0,cap_t=2500, omega_a=0, omega_b=0.5,mark_set=c(-1,1),varepsilon=0.03, N=1500, c_=0, h=0.1, new=F,stop.at.tau=T)
N_sample<-1000
tau<-c()

tau.F1.25<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)
  stopped.chain$s[length(stopped.chain$s)]
})

write.csv(x=tau.F1.25, file="tau.F1.25.csv")
pushover(message = "tau.F1.25.csv is done")

parameters$c_<-optimize(min_max,interval=c(c_H(),1000),maximum =TRUE, cap_t=parameters$cap_t,low_t=parameters$low_t,omega_a=parameters$omega_a,omega_b=parameters$omega_b,x0=parameters$x0)$maximum
parameters$new<-T

tau.F1.25_new<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)
  c(stopped.chain$s[length(stopped.chain$s)],stopped.chain$L[length(stopped.chain$s)])
})

write.csv(x=tau.F1.25_new, file="tau.F1.25_new.csv")
pushover(message = "tau.F1.25_new.csv is done")

benchmarks.25<-seq(from=10,to=100, by=10)
p_tau.25<-rep(0,length(benchmarks.25))
for (i in 1:length(benchmarks.25)){
  p_tau.25[i]<-mean(1*(tau.F1.25<benchmarks.25[i]))
}

relative.error.25<-sqrt(((1/p_tau.25)-1)/N_sample)

p_tau.new.25<-rep(0,length(benchmarks.25))
for (i in 1:length(benchmarks.25)){
  p_tau.new.25[i]<-mean(1*(tau.F1.25_new[1,]<benchmarks.25[i])*tau.F1.25_new[2,])
}

relative.error.new.25<-sqrt(((1/p_tau.new.25)-1)/N_sample)
results25<-rbind(p_tau.new.25,relative.error.25,relative.error.new.25)
colnames(x = results25)<-benchmarks.25

pushover(message = "The exit time when x_0=0.25 is done")

############################################
#Fitting the first curve

curve<-seq(from=1,to=250,by = 0.1)
p.curve25<-c()
for (i in 1:length(curve)){
  p.curve25[i]<-mean(1*(tau.F1.25_new[1,]<curve[i])*tau.F1.25_new[2,])
}

fittingcurve25<-data.frame(x=curve,y=log(1-p.curve25))
model25<-lm(y~x-1,data = fittingcurve25)
model25$coefficients

#########################################
#

parameters<-list(x0=0.75,low_t=0,cap_t=2500, omega_a=0.5, omega_b=1,mark_set=c(-1,1),varepsilon=0.03, N=1500, c_=0, h=0.1, new=F,stop.at.tau=T)
N_sample<-1000

tau.F1.75<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)
  stopped.chain$s[length(stopped.chain$s)]
})

write.csv(x=tau.F1.75, file="tau.F1.75.csv")
pushover(message = "tau.F1.75.csv is done")

parameters$c_<-optimize(min_max,interval=c(c_H(),1000),maximum =TRUE, cap_t=parameters$cap_t,low_t=parameters$low_t,omega_a=parameters$omega_a,omega_b=parameters$omega_b,x0=parameters$x0)$maximum
parameters$new<-T

tau.F1.75_new<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)
  c(stopped.chain$s[length(stopped.chain$s)],stopped.chain$L[length(stopped.chain$s)])
})

write.csv(x=tau.F1.75_new, file="tau.F1.75_new.csv")
pushover(message = "tau.F1.75_new.csv is done")

benchmarks.75<-seq(from=10,to=100, by=10)
p_tau.75<-rep(0,length(benchmarks.75))
for (i in 1:length(benchmarks.75)){
  p_tau.75[i]<-mean(1*(tau.F1.75<benchmarks.75[i]))
}

relative.error.75<-sqrt(((1/p_tau.75)-1)/N_sample)

p_tau.new.75<-rep(0,length(benchmarks.75))
for (i in 1:length(benchmarks.75)){
  p_tau.new.75[i]<-mean(1*(tau.F1.75_new[1,]<benchmarks.75[i])*tau.F1.75_new[2,])
}

relative.error.new.75<-sqrt(((1/p_tau.new.75)-1)/N_sample)
results75<-rbind(p_tau.new.75,relative.error.75,relative.error.new.75)
colnames(x = results75)<-benchmarks.75

pushover(message = "The exit time when x_0=0.75 is done")

############################################
#Fitting the second curve

curve<-seq(from=1,to=250,by = 0.1)
p.curve75<-c()
for (i in 1:length(curve)){
  p.curve75[i]<-mean(1*(tau.F1.75_new[1,]<curve[i])*tau.F1.75_new[2,])
}

fittingcurve75<-data.frame(x=curve,y=log(1-p.curve75))

model2<-lm(y~x-1,data = fittingcurve75)
model2$coefficients

############################################
