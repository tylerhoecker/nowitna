# Load Ryan Kelly functions for ZIL procedure ----------------------------------
# Wrapper function returning kernel weights (K), standardized weights (w), and the effective sample size (n). 
ry.kernel = function(x, x0, h, kernel.fn) {
  dists = abs(x-x0)
  
  K  = kernel.fn(dists,h)
  #		K[which(K<1e-10)] = 0; # Minimum wernel weight to be considered (think machine precision)
  K0 = kernel.fn(0,h)
  
  if(sum(K)==0) { 
    # The total weight of the data is 0 (e.g. x0 is too far away). Assign manually because can't divide by sum(K).
    w = numeric(length(K));
  } else {
    # Otherwise, just divide by sum(K) so weights add to 1
    w = K/sum(K);
  }
  
  n = sum(K/K0)
  # Derived independently but equivalently to:
  #		Chaudhuri, P., and J. S. Marron. 1999. SiZer for exploration of structures in curves. Journal of the American Statistical Association:807–823. PAGE 812
  
  return(list(w=w,n=n,K=K))
}

#------------------ Individual kernels
# Uniform kernel with half-window width = h
kernel.unif = function(dists, h) {
  inWin = dists <= h
  
  K = numeric(length(dists))
  K[inWin] = 1
  
  return(K)
}

# Gaussian kernel with SD=h
kernel.gauss = function(dists, h) {
  K = dnorm(dists,0,h)
  return(K)
}

#------------------ Parameter order for all funtions below
# Zero-inflated lognormal
# parms[1] = p 			= probability of zero
# parms[2] = mu 		= mean of lognormal dist
# parms[3] = sigma  = sd of lognormal dist

#------------------ ZIlnorm.randomdraw
# Produce a random draw from the zero-inflated lognormal distribution
ZIlnorm.randomdraw = function(parms,n) {
  (runif(n,0,1)>parms[1]) * rlnorm(n, parms[2], parms[3])
}
#------------------ ZIlnorm.mle

# Analytical solution for MLE parameters of zero-inflated log-normal distribution
ZIlnorm.mle = function(sampledata) {
  n=length(sampledata)
  mle.lnorm.p    = sum(sampledata==0)/n
  mle.lnorm.mean = mean(log(sampledata[sampledata>0])) 
  mle.lnorm.sd   = 
    sqrt( sum( (log(sampledata[sampledata>0]) - mle.lnorm.mean)^2 )/sum(sampledata>0) ) 
  mleparms 			 = c(mle.lnorm.p, mle.lnorm.mean, mle.lnorm.sd)
  return(mleparms)
}
#------------------ ZIlnorm.llk

# Log-likelihood for zero-inflated lognormal
ZIlnorm.llk = function(parms,x,w) {
  # Weight vector w can be omitted, in which case values are given equal weight
  if( missing(w) ) { w = rep(1,length(x)) }
  
  # This code works in several separated steps:
  #	ind0 = which(x==0)
  #	ind1 = which(x>0)
  #
  #	llike0 = sum( w[ind0] * log(parms[1]) )
  #	llike1 = sum( w[ind1] * (log(1-parms[1]) + dlnorm(x[ind1],parms[2],parms[3], log=TRUE)) )
  #
  #	return( -1* (llike0 + llike1) )
  
  # But this seems a bit faster:
  return( -1* (
    sum( w[x==0] * log(parms[1]) ) +
      sum( w[x>0] * (log(1-parms[1]) + dlnorm(x[x>0],parms[2],parms[3], log=TRUE)) )
  ))
}
ZIlnorm.llk2 = function(parms,x,w) {
  # Weight vector w can be omitted, in which case values are given equal weight
  if( missing(w) ) { w = rep(1,length(x)) }
  
  return( -1*sum( w[x>0] * dlnorm(x[x>0],parms[1],parms[2], log=TRUE) ) )
}
#------------------ ZIlnorm.mean

# Calculate the mean of a specified ZIlnorm distribution. Works on a 3-parameter vector or an array with 3 parameters per row. 
ZIlnorm.mean = function(parms) {
  if( length(parms) == 3) { # it's just a vector of three parms
    return(  (1-parms[1]) * (exp(parms[2] + (parms[3]^2)/2))  )
  } else { # it's an array with 3 parms per row
    n.rows = dim(parms)[1]
    output = numeric(n.rows)
    for( i in 1:n.rows) {
      output[i] = (1-parms[i,1]) * (exp(parms[i,2] + (parms[i,3]^2)/2))
    }
    return(output)
  }
}
ZIlnorm.median = function(parms) {
  if(parms[1]>=0.5) { # at least 50% of points == 0
    return(0)
  } else {		
    return(  qlnorm( ((0.5-parms[1])/(1-parms[1])) ,parms[2],parms[3])  )
  }
}()

# REPEAT THESE FOR EACH RUN -TH 10/18
# Parameters --------------------------------------------------------------
# Temporal params
ageLim      <- c(1950-studyPeriod[2],1950-studyPeriod[1])
xxStep      <- timeStep				  # Interval for fitted x's
bandWidth   <- 25				  # half-window width
kernel.name	<- 'gauss'      # Kernel to use (a function from kernels.r)
# Confidence Intervals parameters
compute.CI  <- 1 
nboot       <- 1000
alpha       <- 0.1


# Number of records

# Bootsrap