#############
# Functions #       
#############

H.px<-function(p,x){
  exponent<-diag(x=p,nrow=length(p),ncol=length(p))%*%t(replicate(length(p),mark_set))
  b(x)*p+(0.5*sigma_fn(x)^2)*p^2+rowSums(exp(exponent)-1-exponent)
}

c_H<-function(){
  inside<-function(x){
    optimize(f=H.px,interval = c(0,1000),maximum = FALSE,x=x)$objective
  }
  optimize(f=inside,interval = c(0,1.05696),maximum = TRUE)$objective
}

p.fn<-function(z,c_,max){
  H.xpc<-function(p,xx,c_){
    exponent<-diag(x=p,nrow=length(p),ncol=length(p))%*%t(replicate(length(p),mark_set))
    b(xx)*p+(0.5*sigma_fn(xx)^2)*p^2+rowSums(exp(exponent)-1-exponent)-c_
  }
  pp<-sapply(X = z,FUN=function(x) uniroot.all(f = H.xpc,interval = c(-10,10),xx=x,c_=c_))
  fff<-ifelse(max==TRUE,"max","min")
  if (is.null(dim(pp))){
    return(pp)
  }
  else{
    return(apply(X = pp,MARGIN = 2,FUN = fff))
  }
}

mane<-function(c_,x,y){
  ifelse(x<y,integrate(p.fn,lower = x,upper=y,c_=c_,max=TRUE)$value,integrate(p.fn,lower = x,upper=y,c_=c_,max=FALSE)$value)
  # if(x<y){
  #   integrate(p.fn,lower = x,upper=y,c_=c_,max=TRUE)
  # }
  # else{
  #   integrate(p.fn,lower = x,upper=y,c_=c_,max=FALSE)
  # }
}

min_max<-function(c_,cap_t,low_t,omega_a,omega_b,x0){
  min(mane(c_,x0,omega_a),mane(c_,x0,omega_b))-c_*(cap_t-low_t)
}

Theta<-function(x,c_){
  ifelse(x>x0,sigma_fn(x)*p.fn(x,c_,max=TRUE),sigma_fn(x)*p.fn(x,c_,max=FALSE))
  # if(x>x0){
  #   sigma_fn(x)*p.fn(x,c_,max=TRUE)
  # }
  # else{
  #   sigma_fn(x)*p.fn(x,c_,max=FALSE)
  # }
}

b_tilde<-function(x,c_){
  b(x)+sqrt(varepsilon)*sigma_fn(x)*Theta(x,c_)
}

new_step_original<-function(Z,s,A,E_n,i,L){
  A_temp<-A+lambda(Z)*(low_t+(i+1)*h-s)
  if(A_temp >= E_n){
    tau_n<-s+(E_n-A)/(lambda(Z))
    Z_tau_minus<-Z+b(Z)*(tau_n-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(tau_n-s)*rnorm(n = 1)
    Z<-max(0,min(Z_tau_minus+varepsilon*gamma_fn(y = random_mark(Z_tau_minus)),1))
    s<-tau_n
    A<-E_n
    E_n<-E_n+rexp(n = 1,rate = 1)
  }
  else{
    Z<-Z+b(Z)*(low_t+(i+1)*h-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(low_t+(i+1)*h-s)*rnorm(n = 1)
    Z_tau_minus<-max(0,min(Z,1))
    s<-low_t+(i+1)*h
    A<-A_temp
    i<-i+1
  }
  data.frame(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
}

new_step_after<-function(Z,s,A,E_n,i,L){
  A_temp<-A+lambda(Z)*(low_t+(i+1)*h-s)
  if(A_temp >= E_n){
    tau_n<-s+(E_n-A)/(lambda(Z)) #time the jump actually happened
    Z_tau_minus<-Z+b_tilde(Z,c_)*(tau_n-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(tau_n-s)*rnorm(n = 1)
    L<-L+(Theta(Z,c_)*L*sqrt(tau_n-s)*rnorm(n = 1))/sqrt(varepsilon)
    Z<-max(0,min(Z_tau_minus+varepsilon*gamma_fn(y = random_mark(Z_tau_minus)),1))
    s<-tau_n
    A<-E_n
    E_n<-E_n+rexp(n = 1,rate = 1)
  }
  else{
    Z_tau_minus<-Z+b_tilde(Z,c_)*(low_t+(i+1)*h-s)+sqrt(varepsilon)*sigma_fn(Z)*sqrt(low_t+(i+1)*h-s)*rnorm(n = 1)
    L<-L+(Theta(Z,c_)*L*sqrt((i+1)*h-s)*rnorm(n = 1))/sqrt(varepsilon)
    Z<-max(0,min(Z_tau_minus,1))
    s<-low_t+(i+1)*h
    A<-A_temp
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

  if(parameters$new==FALSE){
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
  }
  else{
      while(eval(parse(text=stopping.condition))){
      current_step<-new_step_after(Z=Z,s=s,A=A,E_n=E_n,i=i,L=L)
      Z<-max(0,min(current_step$Z,1))
      s<-current_step$s
      #print(c(Z,s))
      A<-current_step$A
      E_n<-current_step$E_n
      i <-current_step$i
      L<-current_step$L  
      Z_t<-rbind(Z_t,current_step)
    }
  }
  return(Z_t)
}


####################################### TESTING #####
# 
# 
# chain.independent<-function(Z,s,L,T_,epsilon,new,stop.at.tau){
#   x_0<-x0
#   if(stop.at.tau==TRUE){stopping.condition<-"(s<T_)&&(omega_a<Z)&&(Z<omega_b)"}
#   else {stopping.condition<-"s<T_"}
#   Z_t<-data.frame(Z=Z,s=s,L=L)
#   if(new==FALSE){
#     while(eval(parse(text=stopping.condition))){
#       number.jumps.interval<-rpois(n = 1,lambda = lambda(x=Z)*h)
#       Z<-Z+b(Z)*h+sqrt(epsilon)*sigma_fn(Z)*sqrt(h)*rnorm(n = 1)+epsilon*sum(sample(x =mark_set,replace =TRUE,size = number.jumps.interval))
#       Z<-max(0,min(1,Z))
#       s<-round(s+h,decimalplaces(h))
#       Z_t<-rbind(Z_t,data.frame(Z=Z,s=s))
#     }  
#   }
#   else{
#     while(eval(parse(text=stopping.condition))){
#       number.jumps.interval<-rpois(n = 1,lambda = lambda(x=Z)*h)
#       Z<-Z+b_tilde(Z,varepsilon,c_,x_0)*h+sqrt(epsilon)*sigma_fn(Z)*sqrt(h)*rnorm(n = 1)+epsilon*sum(sample(x =mark_set,replace =TRUE,size = number.jumps.interval))
#       L<-L+(Theta(Z,c_,x_0)*L*sqrt(h)*rnorm(n = 1))/sqrt(varepsilon)
#       Z<-max(0,min(1,Z))
#       s<-round(s+h,decimalplaces(h))
#       Z_t<-rbind(Z_t,data.frame(Z=Z,s=s,L=L))
#     }
#   }
#   return(Z_t)
# }

decimalplaces <- function(x) {
  if ((x %% 1) != 0) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
