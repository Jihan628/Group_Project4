### Function for minimization test
install.packages('corpcor')
library(corpcor)
rb <- function(th,k=2) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}

gb <-  function(th,k=2) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}




hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


# Example 2 (Lecture notes, still doesn't work)


rll <- function(theta,t,y) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i)) ## theta = (alpha,beta)
  mu <- theta[1] * exp(theta[2] * t) ## mu = E(y)
  -sum(dpois(y,mu,log=TRUE)) ## the negative log likelihood 
} 

gll <- function(theta,t,y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model 
  alpha <- theta[1]
  beta <- theta[2] ## enhances readability 
  ebt <- exp(beta*t) ## avoid computing twice 
  -c(sum(y)/alpha - sum(ebt), ## -dl/dalpha
     sum(y*t) - alpha*sum(t*ebt)) ## -dl/dbeta 
} 

hll <- function(theta,t,y) {
  ## Hessian of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1]
  beta <- theta[2] ## enhances readability 
  ebt <- exp(beta*t) ## avoid computing twice
  H <- matrix(0,2,2) ## matrix for Hessian of -ve ll
  H[1,1] <- sum(y)/alpha^2
  H[2,2] <- alpha*sum(t^2*ebt)
  H[1,2] <- H[2,1] <- sum(t*ebt)
  H
}


#############################################################################
### Code


newt<-function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  
  ftheta<-grad(theta,...) # gradient
  intvalue<-func(theta,...) # value of the function
  value<-intvalue
  thetanew<-theta 
  iter=0
  dim<-length(theta) #dimension of theta
  
  if(all(is.finite(value) & is.finite(ftheta))){ #Both objective and derivative is finite
  
  if(is.null(hess)){ #If without hess function
    # Compute the f'x and f''x before entering the loop  
    
    
    ftheta<- grad(thetanew,...) 
    H<- matrix(0,dim,dim)
    
    for (i in 1:dim) { ## loop over parameters
      th1 <- thetanew;
      th1[i] <- th1[i] + eps ## increase th0[i] by eps 
      grad1 <- grad(th1,...) ## compute resulting nll
      H[i,] <- (grad1 - ftheta)/eps ## approximate second derivs
    }
    fftheta<-diag(H)
    while (any(abs(ftheta)>=tol*abs(value)+fscale)){ # Criteria listed by the pratical 
      
      
      #Using xi+1=xi-f'x/f''x
      thetanew<-thetanew-ftheta/fftheta
      value<-func(thetanew,...)
      
      #Need to build a hessian matrix
      
      ftheta<- grad(thetanew,...) 
      H<- matrix(0,dim,dim)
      
      for (i in 1:dim) { ## loop over parameters
        th1 <- thetanew;
        th1[i] <- th1[i] + eps ## increase th0[i] by eps 
        grad1 <- grad(th1,...) ## compute resulting nll
        H[i,] <- (grad1 - ftheta)/eps ## approximate second derivs
      }
      fftheta<-diag(H)
      iter<-iter+1
      
      if (iter=max.half & value>=intvalue){
        warning("")
      }
      
      
      
      if (iter > maxit){
        stop("maximum iterations is reached without convergence")
      }
      
      
    }
    
    
    list(f=value,theta=thetanew,iter=iter,g=ftheta,Hi=chol2inv(chol(H)))
    
  } else { # With hess matrix
    fftheta<-diag(hess(thetanew,...)) # Get the diagonal entries for hess
    while (any(abs(ftheta)>=tol*abs(value)+fscale)){
      
      
      thetanew<-thetanew-ftheta/fftheta
      value<-func(thetanew,...)
      ftheta<-grad(thetanew,...)
      fftheta<-diag(hess(thetanew,...))
      iter<-iter+1
      if (iter=max.half & value>=intvalue){
        warning("")
      }
      
      if (iter > maxit){
        stop("maximum iterations is reached without convergence")
      }
    }
    
    
    list(f=value,theta=thetanew,iter=iter,g=ftheta,Hi=chol2inv(chol(hess(thetanew))))
  }

  
  }  
  else {
    stop('the objective or derivatives are not finite at the initial theta')
  }
}

# Test this, works well, results are very close!
newt(c(1.2,1.2),rb,gb)
newt(c(1.2,1.2),rb,gb,hb)

newt(c(1.8,1.3),rb,gb,hb)#test
newt(c(Inf,1),rb,gb)#test

newt(c(-0.5,-0.3),rb,gb)#test
#Error in chol.default(H) :
#the leading minor of order 2 is not positive definite
newt(c(-0.5,-0.3),rb,gb,hb)#test


### But it still shoots off for this one... Not sure why??
t80 <- 1:13 ## years since 1980
y <- c(12,14,33,50,67,74,123,141,165,204,253,246,240)
newt(theta=c(10,0.1),func=rll,grad=gll,hess=hll,t=t80,y=y)




#### What is done
### Should work for any dimensions of theta
### Returns all outputs required
### Hessian matrix computation done


### Still needs to be done
### add maxit and max.half into the function
### Check whether hessian is positive definite at convergence
### check the objective and derivatives are finite at the initial theta
### Other weird scenarios!!!

#-------------------------------------------------------------------------
## Jihan:
## Have add maxit but not sure how to add max.half correctly
## I try to use is.positive.definite but it need package 'corpcor'
## Have add it by is.finite 
## while (any(abs(ftheta) >= tol * abs(value) + fscale)) I don't know how to 
## solve it, maybe NaNs produced in dpois
