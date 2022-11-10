### Function for minimization test

rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <-  function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}


machineprec<-2.220e-16


hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}







newt<-function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  # take the derivative to get the gradient
  ftheta<-grad(theta)
  thetanew<-theta 
  iter=0
  
  if(is.null(hess)){ 
    while (ftheta[1]>0.001|ftheta[1]< (-0.001)|ftheta[2]>0.001|ftheta[2]<(-0.001)){
      
      ftheta<-gb(thetanew) 
      prec<-sqrt(2.220e-16)
      
      #Need to build a hessian matrix
      # Finite derivative
      ffd1d1<-(gb(thetanew+c(prec,0))[1]-gb(thetanew)[1])/prec
      ffd2d2<-(gb(thetanew+c(0,prec))[2]-gb(thetanew)[2])/prec
      ffd1d2<-(gb(thetanew+c(0,prec))[1]-gb(thetanew)[1])/prec
      ffd2d1<-(gb(thetanew+c(prec,0))[2]-gb(thetanew)[2])/prec
      h <- matrix(0,2,2)
      h[1,1]<-ffd1d1
      h[2,2]<-ffd2d2
      h[1,2]<-ffd1d2
      h[2,1]<-ffd2d1
      
      fftheta<-c(ffd1d1,ffd2d2)
      
      #Using xi+1=xi-f'x/f''x
      thetanew<-thetanew-ftheta/fftheta
      value<-func(thetanew)
      iter<-iter+1
      
      
      
    }
    
    
    list(f=value,theta=thetanew,iter=iter,g=ftheta,Hi=chol2inv(chol(h)))
    
  } else {
  while (ftheta[1]>0.001|ftheta[1]< (-0.001)|ftheta[2]>0.001|ftheta[2]<(-0.001)){
  
  
  ftheta<-grad(thetanew)
  fftheta<-c(hess(thetanew)[1,1],hess(thetanew)[2,2])
  
  thetanew<-thetanew-ftheta/fftheta
  value<-func(thetanew)
  iter<-iter+1
  }
    list(f=value,theta=thetanew,iter=iter,g=ftheta,Hi=chol2inv(chol(hess(thetanew))))
  }
  
  
  
  
 
  
}


newt(c(1.05,1.05),rb,gb)
newt(c(1.05,1.05),rb,gb,hb)









