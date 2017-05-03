library(rootSolve)
library(ggplot2)
library(reshape2)
library(pushoverr)

#source('~/Concordia/Research/Git/DGM/functions - Bistability - April.R')

b<-function(x){
  0
}

sigma_fn<-function(x){
  0
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

new_step_original<-function(Z,s,A,E_n,i,L){
  A_temp<-A+lambda(Z)*(low_t+(i+1)*h-s)
  if(A_temp >= E_n){
    tau_n<-s+(E_n-A)/(lambda(Z))
    Z_tau_minus<-Z+b(Z)*(tau_n-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(tau_n-s)*rnorm(n = 1)
    Z<-max(0,min(Z_tau_minus+1*gamma_fn(y = random_mark(Z_tau_minus)),1))
    s<-tau_n
    A<-E_n
    E_n<-E_n+rexp(n = 1,rate = 1)
  }
  else{
    # Z<-Z+b(Z)*(low_t+(i+1)*h-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(low_t+(i+1)*h-s)*rnorm(n = 1)
    # Z_tau_minus<-max(0,min(Z,1))
    # s<-low_t+(i+1)*h
    # A<-A_temp
    i<-i+1
  }
  data.frame(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
}


chain<-function(parameters,b,sigma_fn,lambda,gamma_fn,random_mark){
  x0<<-Z<-parameters$x0
  low_t<<-s<-parameters$low_t
  cap_t<<-parameters$cap_t
  omega_a<<-parameters$omega_a
  omega_b<<-parameters$omega_b
  varepsilon<<-parameters$varepsilon
  mark_set<<-parameters$mark_set
  N<<-parameters$N
  c_<<-parameters$c_
  h<<-parameters$h
  
  A<-0
  E_n<-rexp(n = 1,rate = 1)
  i<-0
  L<-1
  Z_t<-data.frame(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
  
  stopping.condition<-ifelse(parameters$stop.at.tau==TRUE,"(s<cap_t)&&(omega_a<Z)&&(Z<omega_b)","s<cap_t")
  
  #if(parameters$new==FALSE){
    while(eval(parse(text=stopping.condition))){
      current_step<-new_step_original(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
      Z<-max(0,min(current_step$Z,1))
      s<-current_step$s
      #print(c(Z,s))
      A<-current_step$A
      E_n<-current_step$E_n
      i<-current_step$i
      Z_t<-rbind(Z_t,current_step)
    }
  #}
  # else{
  #   while(eval(parse(text=stopping.condition))){
  #     current_step<-new_step_after(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
  #     Z<-max(0,min(current_step$Z,1))
  #     s<-current_step$s
  #     #print(c(Z,s))
  #     A<-current_step$A
  #     E_n<-current_step$E_n
  #     i <-current_step$i
  #     L<-current_step$L  
  #     Z_t<-rbind(Z_t,current_step)
  #   }
  #}
  return(Z_t)
}

set.seed(5678)

parameters<-list(x0=0.25,low_t=0,cap_t=100, omega_a=0, omega_b=0.5,mark_set=c(-1/500,1/500),varepsilon=2e-4, N=500, c_=0, h=0.1, new=F,stop.at.tau=F)

Sys.time()
beg<-proc.time()
prueba<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn,random_mark = random_mark)
end<-proc.time()-beg

p<-ggplot(prueba,aes(x = s,y = Z))
p+geom_line()+ ylab("Z(s)")+geom_hline(yintercept = 0.25,colour = "red")+geom_hline(yintercept = 0.75,colour = "red")

parameters$h<-0.01
parameters$cap_t<-2500

beg2<-proc.time()
prueba2<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn,random_mark = random_mark)
end2<-proc.time()-beg2

p<-ggplot(prueba2,aes(x = s,y = Z))
p+geom_line()+ ylab("Z(s)")+geom_hline(yintercept = 0.25,colour = "red")+geom_hline(yintercept = 0.75,colour = "red")









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

r<-ggplot(prueba, aes(Z))
r+geom_density()+xlab("Z(s)")

pushover(message = "HIEEEE, benchmark2 is done")

N_sample<-1000

tau<-c()
for(i in 1:N_sample){
  stopped.chain<-chain.independent(Z=x0,s = low_t,T_ = T_vec,epsilon = varepsilon_vec,stop.at.tau=TRUE)
  tau<-c(tau,stopped.chain$s[length(stopped.chain$s)])
}

benchmarks<-as.table(summary(tau.F1.25))+1
benchmarks<-c(1.5,2,2.5,3)
p_tau<-rep(0,length(benchmarks))
for (i in 1:length(benchmarks)){
  p_tau[i]<-mean(1*(tau.F1.25<benchmarks[i]))
}

relative.error<-sqrt(((1/p_tau)-1)/N_sample)


##############################Alternative algorithm########################


parameters<-list(x0=0.25,low_t=0,cap_t=2500, omega_a=0, omega_b=0.5,mark_set=c(-1/500,1/500),varepsilon=2e-4, N=500, c_=0, h=0.1, new=F,stop.at.tau=F)

regular_chain<-function(parameters){
  x<<-parameters$x0
  low_t<<-s<-parameters$low_t
  cap_t<<-parameters$cap_t
  omega_a<<-parameters$omega_a
  omega_b<<-parameters$omega_b
  varepsilon<<-parameters$varepsilon
  mark_set<<-parameters$mark_set
  N<<-parameters$N

  stopping.condition<-ifelse(parameters$stop.at.tau==TRUE,"(s<cap_t)&&(omega_a<x)&&(x<omega_b)","s<cap_t")
  
  X_s<-data.frame(x,s)
  while(eval(parse(text=stopping.condition))){
    cap_lambda<-(1/2)*varepsilon*x*(1-x)*N^2
    r_plus<-1-x+(32/3)*x^2*(1-x)
    r_minus<-x+(16/3)*x*(1-x)
    ds <- rexp(1,rate=cap_lambda+r_plus+r_minus)
    event <- sample(x=mark_set,size = 1,replace = FALSE,prob=c(r_minus+0.5*cap_lambda,r_plus+0.5*cap_lambda))
    s <- s + ds
    x <- min(max(0,x+event),1)
    current_step<-data.frame(x,s)
    X_s<-rbind(X_s,current_step)
  }
  return(X_s)
} 

prueba2<-regular_chain(parameters)

p<-ggplot(prueba2,aes(x = s,y = x))
p+geom_line()+ ylab("Z(s)")+geom_hline(yintercept = 0.25,colour = "red")+geom_hline(yintercept = 0.75,colour = "red")

#################Exit probability when x0=0.25################

parameters$stop.at.tau<-T
N_sample<-1000
Sys.time()
begtime<-proc.time()
tau<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<<-chain(parameters=parameters,b=b,sigma_fn=sigma_fn,lambda=lambda,gamma_fn=gamma_fn, random_mark=random_mark)
  stopped.chain$s[length(stopped.chain$s)]
})
endtime<-proc.time()-begtime
#Si lo hago asi van a ser 43.88889 horas MINIMO



Sys.time()
begtime<-proc.time()
tau.alternative<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<<-regular_chain(parameters=parameters)
  stopped.chain$s[length(stopped.chain$s)]
})
endtime<-proc.time()-begtime
#Asi son 10 horas MINIMO

write.csv(x=tau.alternative,file="tau.alternative.csv")



benchmarks.fig2<-seq(from=200,to=1000,by=100)

p_tau.fig2<-rep(0,length(benchmarks.fig2))
for (i in 1:length(benchmarks.fig2)){
  p_tau.fig2[i]<-mean(1*(tau.alternative<benchmarks.fig2[i]))
}

curve<-seq(from=200, to=2000, by=1)

p.curvefig225<-rep(0,length(curve))
for (i in 1:length(curve)){
  p.curvefig225[i]<-mean(1*(tau.alternative<curve[i]))
}

fittingcurvefig225<-data.frame(x=curve,y=-log(1-p.curvefig225))
modelfig225<-lm(y~x-1,data = fittingcurvefig225)
modelfig225$coefficients


#################Exit probability when x0=0.75################

parameters$x0<-0.75
parameters$stop.at.tau<-T
N_sample<-1000

Sys.time()
begtime<-proc.time()
tau.alternative75<-sapply(X=1:N_sample,FUN=function(iterNumber){
  stopped.chain<<-regular_chain(parameters=parameters)
  stopped.chain$s[length(stopped.chain$s)]
})
endtime<-proc.time()-begtime
Sys.time()
#Asi son 10 horas MINIMO

write.csv(x=tau.alternative75,file="tau.alternative75.csv")



benchmarks.fig2<-seq(from=200,to=1000,by=100)

p_tau.fig275<-rep(0,length(benchmarks.fig275))
for (i in 1:length(benchmarks.fig275)){
  p_tau.fig275[i]<-mean(1*(tau.alternative75<benchmarks.fig275[i]))
}

curve<-seq(from=200, to=2000, by=1)

p.curvefig275<-rep(0,length(curve))
for (i in 1:length(curve)){
  p.curvefig275[i]<-mean(1*(tau.alternative75<curve[i]))
}

fittingcurvefig275<-data.frame(x=curve,y=-log(1-p.curvefig275))
modelfig275<-lm(y~x-1,data = fittingcurvefig275)
modelfig275$coefficients

