#############
# Functions #       
#############

b<-function(x){
  0
}

sigma_fn<-function(x){
  x
}

lambda<-function(x){
  lambda_0+lambda_1*x
}

gamma_fn<-function(y){
  1*y
}

random_mark<-function(){
  sample(x =mark_set,replace = FALSE,size = 1)
}

H.px<-function(p,x){
  exponent<-diag(x=p,nrow=length(p),ncol=length(p))%*%t(replicate(length(p),mark_set))
  b(x)*p+(0.5*sigma_fn(x)^2)*p^2+rowSums(exp(exponent)-1-exponent)
}

c_H<-function(){
  inside<-function(x){
    optimize(f=H.px,interval = c(-10000,10000),maximum = FALSE,x=x)$objective
  }
  optimize(f=inside,interval = c(-10000,10000),maximum = TRUE)$objective
}


p.fn<-function(z,c_,max){
  H.xpc<-function(p,xx,c_){
    exponent<-diag(x=p,nrow=length(p),ncol=length(p))%*%t(replicate(length(p),mark_set))
    b(xx)*p+(0.5*sigma_fn(xx)^2)*p^2+rowSums(exp(exponent)-1-exponent)-c_
  }
  pp<-sapply(X = z,FUN=function(x) uniroot.all(f = H.xpc,interval = c(-10,10),xx=x,c_=c_))
  if (max==TRUE){fff<-"max"}
  else {fff<-"min"}
  if (is.null(dim(pp))){
    return(pp)
  }
  else{
    return(apply(X = pp,MARGIN = 2,FUN = fff))
  }
}

mane<-function(c_,x,y){
  if(x<y){
    integrate(p.fn,lower = x,upper=y,c_=c_,max=TRUE)
  }
  else{
    integrate(p.fn,lower = x,upper=y,c_=c_,max=FALSE)
  }
}

min_max<-function(c_,cap_t,low_t,a,b_){
  min(mane(c_,x0,a)$value,mane(c_,x0,b_)$value)-c_*(cap_t-low_t)
}

Theta<-function(x,c_){
  if(x>x0){
    sigma_fn(x)*p.fn(x,c_,max=TRUE)
  }
  else{
    sigma_fn(x)*p.fn(x,c_,max=FALSE)
  }
}

b_tilde<-function(x,varepsilon,c_){
  b(x)+sqrt(varepsilon)*sigma_fn(x)*Theta(x,c_)
}

new_step_original<-function(Z,s,A,E_n,i,L,h,varepsilon,TT){
  A_temp<-A+lambda(Z)*((i+1)*h-s)
  N<-rnorm(n = 1,mean = 0,sd = 1)
  if(A_temp >= E_n){
    tau_n<-s+(E_n-A)/(lambda(Z)) #time the jump actually happened
    Z_tau_minus<-Z+b(Z)*(tau_n-s)+varepsilon*sigma_fn(Z)*sqrt(tau_n-s)*N
    Z<-Z_tau_minus+varepsilon*gamma_fn(y = random_mark())
    s<-tau_n
    A<-E_n
    E_n<-E_n+rexp(n = 1,rate = 1)
  }
  else{
    Z<-Z+b(Z)*((i+1)*h-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt((i+1)*h-s)*N
    Z_tau_minus<-Z
    s<-(i+1)*h
    A<-A_temp
    i<-i+1
  }
  data.frame(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
}

new_step_after<-function(Z,s,A,E_n,i,L,h,varepsilon,TT,c_){
  A_temp<-A+lambda(Z)*((i+1)*h-s)
  N<-rnorm(n = 1,mean = 0,sd = 1)
  if(A_temp >= E_n){
    tau_n<-s+(E_n-A)/(lambda(Z)) #time the jump actually happened
    Z_tau_minus<-Z+b_tilde(Z,varepsilon,c_)*(tau_n-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(tau_n-s)*N
    L<-L+(Theta(Z,c_)*L*sqrt(tau_n-s)*N)/sqrt(varepsilon)
    Z<-Z_tau_minus+varepsilon*gamma_fn(y = random_mark())
    s<-tau_n
    A<-E_n
    E_n<-E_n+rexp(n = 1,rate = 1)
  }
  else{
    Z_tau_minus<-Z+b_tilde(Z,varepsilon,c_)*((i+1)*h-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt((i+1)*h-s)*N
    L<-L+(Theta(Z,c_)*L*sqrt((i+1)*h-s)*N)/sqrt(varepsilon)
    Z<-Z_tau_minus
    s<-(i+1)*h
    A<-A_temp
    i<-i+1
  }
  data.frame(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
}


chain<-function(x0,t0,T_,epsilon,Nsteps,c_,new){
  h<-(T_-t0)/Nsteps
  Z<-x0
  s<-t0
  A<-0
  E_n<-rexp(n = 1,rate = 1)
  i<-0
  L<-1
  Z_t<-data.frame(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
  if(new==FALSE){
    while((s<T_)&&(Omega_a<Z)&&(Z<Omega_b)){
      #while(s<T_){
      current_step<-new_step_original(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L,h=h,varepsilon=epsilon,TT=T_)
      Z<-current_step$Z
      s<-current_step$s
      A<-current_step$A
      E_n<-current_step$E_n
      i<-current_step$i
      Z_t<-rbind(Z_t,current_step)
    }
  }
  else{
    while((s<T_)&&(Omega_a<Z)&&(Z<Omega_b)){
      #while(s<T_){
      current_step<-new_step_after(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L,h=h,varepsilon=epsilon,TT=T_,c_=c_)
      Z<-current_step$Z
      s<-current_step$s
      A<-current_step$A
      E_n<-current_step$E_n
      i <-current_step$i
      L<-current_step$L  
      Z_t<-rbind(Z_t,current_step)
    }
  }
  return(Z_t)
}

