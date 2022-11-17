## This file contains the solution to Practical 4 by Group34 ::
## Evangelos Antypas (s2449453), Daniel Kwok (s2308472), Jihan Li(s2322347)
## The github repository containing the solution can be found at ::
## "https://github.com/Jihan628/Group_Project4"

##-----------------------------------------------------------------------------

## Contributions :: To be completed

##-----------------------------------------------------------------------------

## Our goal is to create a function 'newt' implementing Newton optimization,
## which is based on the idea of fitting a quadratic approximation to our 
## objective function f at the trial value xk, having the same slope and 
## curvature as f at that point, and then proceeding to the minimum of the 
## quadratic approximation, in this way we create a sequence of values that, we 
## hope, is converging to the minimum of f. 

##-----------------------------------------------------------------------------

## More specifically, the procedure described above is theoretically justified 
## as follows ::
## Given a twice differentiable function f (note that f could be a 
## multidimensional function), we can obtain a quadratic approximation for f 
## by computing the Taylor expansion at a point xk::
## f(xk+t) = f(xk) + f'(xk)t + (1/2)f''(xk)t^2 + error, so now we can use this
## approximation to obtain the next value in our sequence, by finding the min
## for the quadratic::
## d/dt(f(xk) + f'(xk)t + (1/2)f''(xk)t^2) = 0 -> t_min = -(f''(xk))^(-1)f'(xk)
## so, the recursive formula for our method is given by::
## xk+1 = xk - (f''(xk))^(-1)f'(xk), in higher dimensions f'' is the Hessian
## matrix, and f' is the gradient of f.

##-----------------------------------------------------------------------------

## This function is used to check the positive definiteness of the Hessian 
## in the newt function.

## Defining the is.posdef function.

is.posdef<-function(A){
  
  ## Inputs :: A square matrix A 
  ##---------------------------------------------------------------------------
  ## Outputs :: The function checks if A is positive definite by trying the 
  ## Cholesky decomposition. If the decomposition exists then A is positive   
  ## definite and the function returns A, otherwise we perturb A until it 
  ## becomes positive definite and return the perturbed version.
  ##---------------------------------------------------------------------------
  
  # Checking if Α is positive definite 
  
  result<-tryCatch(chol(A),error=function(e) e) 
  bool<-inherits(result,"matrix") ## Boolean
  dim<-dim(A)[1]
  
  if (is.null(dim)==FALSE){
  temp<- 1e-6*norm(A,type = 'F')*diag(1,dim) ## Used for perturbing A 
  }
  ## If A is not positive definite we perturb it by adding a 
  ## multiple of the identity, until it becomes positive definite
  
  ## We store the initial A because in case we want to perturb A multiple 
  ## times, we don't want to perturb the already perturbed version but A
  ## itself.
  
  A_init<-A
  
  
  while(bool==FALSE){
    
    A<-A_init+temp
    
    result<-tryCatch(chol(A),error=function(e) e) 
    bool<-inherits(result,"matrix") ## Boolean
    temp<-10*temp ## If our perturbation didn't work we multiply the 
    ## magnitude by 10
    
  }
  
  return(A)
  
}

##-----------------------------------------------------------------------------

## Defining the newt function

newt<-function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
               maxit=100,max.half=20,eps=1e-6){
  
  
  ## Inputs :: theta: our initial point, func: the objective function to 
  ## be minimized, grad: the gradient of ff, hess: the Hessian of obj, which by
  ## default is not given and will be computed inside the function by finite
  ## differences, tol: is our tolerance for the magnitude of the components 
  ## of the gradient(we need that to check for convergence), fscale: a rough 
  ## estimate of the magnitude of f near the optimum(used in convergence testing)
  ## maxit: the maximum number of iterations, max.half:  the maximum number of 
  ## times a step should be halved before concluding that the step has failed to
  ## improve obj, eps: the length of the finite difference intervals for the 
  ## computation of the Hessian, any other necessary parameters are passed using
  ## the '...' notation.
  
  ##-----------------------------------------------------------------------------
  
  ## Outputs :: a list containing: f: the value of the objective at the minimum
  ## theta: the value of the parameters at the minimum, iter: the number of 
  ## iterations taken to reach the minimum, g: the gradient vector at the minimum
  ## (this is used to judge closeness to the 0 of the machine), Hi: the inverse 
  ## of the Hessian matrix at the minimum. 
  
  ##-----------------------------------------------------------------------------
  
  ## Approach :: In the default case where the Hessian is not provided by the 
  ## user, we compute it by using the method of finite differences, otherwise we
  ## use the input, in any case the next step is to check for positive
  ## definiteness using the Cholesky decomposition, if our Hessian is not 
  ## positive definite we perturb it by adding l*norm(H)*I, where l is a small 
  ## real value and I is the identity matrix, in order to make it positive 
  ## definite. As soon as we have confirmed positive definiteness we know that 
  ## the method converges.
  
  ##-----------------------------------------------------------------------------
  
  ## Initializations
  
  ftheta<-grad(theta,...)   ## Initial value of the gradient
  intvalue<-func(theta,...) ## Initial value of the objective
  
  value<-intvalue
  thetanew<-theta 
  
  iter=0             ## Number of iterations
  dim<-length(theta) ## Dimension of theta
  
  ##-----------------------------------------------------------------------------
  
  if(all(is.finite(value) & is.finite(ftheta))){ ## Checking for finiteness of 
    ## objective and gradient
    
    if(is.null(hess)){ ## Default case 
      
      
      ## Computing the Hessian at the initial value
      
      
      H<- matrix(0,dim,dim)
      
      for (i in 1:dim) { ## Looping over parameters
        th1 <- theta;
        th1[i] <- th1[i] + eps             ## Increase th0[i] by eps 
        grad1 <- grad(th1,...)             
        H[i,] <- (grad1 - ftheta)/eps      ## Second derivatives approximation
      }
      
      ## Checking if the Hessian is positive definite at the initial value and
      ## if is not we perturb it until it becomes positive definite.
      
      Hnew<-is.posdef(H)
      
      while (any(abs(ftheta)>=(tol*(abs(value)+fscale)))){ ## Gradient components  
        ## sufficiently small
        H<-Hnew
        ## Updating with x_i+1 = x_i + inv_hess*grad
        
        thetanew<-theta-backsolve(chol(H),forwardsolve(t(chol(H)),ftheta)) ## Obtaining the next value
        value<-func(thetanew,...)
        
        ## Check if we are actually improving, otherwise half the step 
        
        half_it <- 0 ## Counter for halvings
        
        while(intvalue<value|is.finite(value)==FALSE){
          half_it<-half_it+1
          
          if (half_it > max.half){
            
            if (is.finite(value)==FALSE){
              stop("maximum halving iterations reached with an
                   infinite value")
            } else {
            
            stop("maximum halving iterations reached without convergence")
            }
          }
          
        
          thetanew<-theta-(1/2)^(half_it+1)*backsolve(chol(H),forwardsolve(t(chol(H)),ftheta)) ## Half the step
          value<-func(thetanew,...)
          
          
          
        } 
        
        ## We have obtained the next value
        
        ## We compute the Hessian for the next value
        
        ftheta<- grad(thetanew,...) 
        H<- matrix(0,dim,dim)
        
        for (i in 1:dim) { ## Looping over parameters
          th1 <- thetanew;
          th1[i] <- th1[i] + eps        ## Increase th0[i] by eps 
          grad1 <- grad(th1,...) 
          H[i,] <- (grad1 - ftheta)/eps ## Approximate second derivatives
        }
        
        ## We are checking if the Hessian is positive definite, if not we perturb
        
        Hnew<-is.posdef(H)
        
        iter<-iter+1
        intvalue<-value
        theta<-thetanew
        
        
        if (iter > maxit){
          stop("maximum iterations is reached without convergence")
        }
        
        
      }
      
      ## Checking whether Hessian is positive definite at convergence
      
      Outcome<-tryCatch(chol(H),error=function(e) e) 
      Bool_H<-inherits(Outcome,"matrix") ## Boolean
      
      if (Bool_H == FALSE){
        warning('Hessian is not
        positive definite at convergence.')
        list(f=value,theta=thetanew,iter=iter,g=ftheta)
        
      } else {
        
        
        list(f=value,theta=thetanew,iter=iter,g=ftheta,Hi=chol2inv(chol(H)))
      }
      
    } else { ## If Hessian is given, proceed as before without finite differences
      
      H<-hess(theta,...)
      
      Hnew<-is.posdef(H)
      
      while (any(abs(ftheta)>=(tol*(abs(value)+fscale)))){
        
        H<-Hnew
        thetanew<-theta-backsolve(chol(H),forwardsolve(t(chol(H)),ftheta))
        value<-func(thetanew,...)
        half_it <- 0
        while(intvalue<value|is.finite(value)==FALSE){
          
          half_it<-half_it+1
          
          if (half_it > max.half){
            
            if (is.finite(value)==FALSE){
              stop("maximum halving iterations reached with an
                   infinite value")
            } else {
              
              stop("maximum halving iterations reached without convergence")
            }
          }
          
          thetanew<-theta-(1/2)^(half_it+1)*backsolve(chol(H),forwardsolve(t(chol(H)),ftheta)) ## Half the step
          value<-func(thetanew,...)
          
          
        } 
        ## Computing gradient & Hessian for the next value
        
        ftheta<-grad(thetanew,...)
        H<-hess(thetanew,...)
        
        ## Check if the Hessian is positive definite
        
        Hnew<-is.posdef(H)
        
        iter<-iter+1
        intvalue<-value
        theta<-thetanew
        
        if (iter > maxit){
          stop("maximum iterations reached without convergence")
        }
      }
      
      ## Checking if the Hessian is positive definite at convergence
      
      Outcome<-tryCatch(chol(H),error=function(e) e) 
      Bool_H<-inherits(Outcome,"matrix") ## Boolean
      
      if (Bool_H == FALSE){
        warning('Hessian is not
        positive definite at convergence')
        
        list(f=value,theta=thetanew,iter=iter,g=ftheta)
        
      } else {
        
        
        list(f=value,theta=thetanew,iter=iter,g=ftheta,Hi=chol2inv(chol(H)))
      }
    }
    
    
  }  
  else {
    stop('the objective or derivatives are not finite at the initial theta')
  }
}