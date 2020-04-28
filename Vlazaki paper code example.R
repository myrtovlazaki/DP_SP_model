## Example script to perform Minimum Divergence Estimation (MDE) of parameters of a stochastic system in single-phenotype and a dual-phenotype model
##
## The example given below reproduces the results of the ampicillin-treated group for both the SP- and DP- model during day 1, in "A data-based mathematical modelling study to quantify the effects of ciprofloxacin and ampicillin on the within-host dynamics of Salmonella enterica during treatment and relapse", by Vlazaki et al. (2020)
## Sub-ordinate functions can be found in SPEEDI.R package
##
## Written by: Myrto Vlazaki, Callum McLean, David Price and Olivier Restif*
## Contact: mv382@cam.ac.uk
## Date: 23th of April 2020

# Moments for ampicillin-treated group (Rossi et al., 2017) provided in file:
# Rossi.amp.wt.moments.correct.R

# Moments for ciprofloxacin-treated group (Rossi et al., 2017) provided in file:
# Rossi.cip.wt.moments.correct.R


install.packages("SPEEDI.R")
library(SPEEDI.R)


# Function to calculate the KL divergence between predicted and experimental moments

KL.div.o2p <- 
function(N, obs,pred,subs=1:N){
  # Predicted moments
  mu0 <- pred[subs]
  cov0 <- v.2.m(pred[(N+1):length(pred)], N)[subs,subs]
  # Observed moments
  mu1<- obs[subs]
  cov1 <- v.2.m(obs[(N+1):length(pred)], N)[subs,subs]
  # KL divergence
  KL.div(mu0,cov0,mu1,cov1)
}



KL.div <- 
function(mu0,cov0,mu1,cov1){
  k <- length(mu0)
  # inv.C1 <- solve(cov1)
  # as.numeric(sum(diag(inv.C1 %*% cov0)) + (mu1-mu0) %*% inv.C1 %*% (mu1-mu0) - k + log(det(cov1)/det(cov0)))/2
  as.numeric(sum(diag(apply(cov0, 2, FUN = function(x){return(solve(cov1,x))}))) + (mu1-mu0) %*% solve(cov1, (mu1-mu0)) - k + log(det(cov1)/det(cov0)))/2
}


# Function to calculate divergence between predicted and experimental moments

div.measure  <- 
function(N,obs,pred,subs=1:N, div){
  switch(div,
         KL = KL.div.o2p(N, obs,pred,subs=subs),
         Hell = Hell.dist.o2p(N, obs,pred,subs=subs),
         Maha = Maha.dist.o2p(N, obs,pred,subs=subs),
         Chisq = Chisq.dist.o2p(obs,pred))
}


# Function to facilitate identification and replacement of parameters by name

replace.par 
function(x,rep)
{
  if(length(rep)==0) return(x)
  pos <- sapply(names(rep),function(n) which(names(x)==n))
  replace(x,pos,rep)
}


# Function that defines the structure of a network with N compartments

all.poss <- 
function(N){
  
  # All migrations allowed
  migrations <- matrix(1,N,N)
  
  # (except, of course, migrations from a compartment to itself)
  diag(migrations) <- 0
  
  # Birth + death allowed in all compartments
  kills <- rep(1,N)
  replications <- rep(1,N)
  
  return(list('migrations'=migrations, 'kills'=kills, 'replications'=replications))}




# Function that yields the parameter names relevant to a compartmental model with N compartments

network.params <- 
function(network.structure){
  
  capacitated <- !is.null(network.structure$capacities)
  
  # Split network structure into migrations, kills, replications
  migrations <- network.structure$migrations
  kills <- network.structure$kills
  replications <- network.structure$replications
  
  # Check if network is capacitated, act appropriately
  if (capacitated){
    capacities <- network.structure$capacities
  }
  
  else {capacities <- 0}
  
  # Get number of compartments
  N = length(kills)
  
  # Get total number of parameters
  total.parameters = sum(migrations) + sum(kills) + sum(replications) + sum(capacities)
  
  # Initialise vector of names
  out.names <- vector('character', total.parameters)
  
  # Start position tracker at 1
  current.parameter <- 1
  
  # Add permitted migrations
  # If migration permitted, add its name to the name vector
  # and increment the position tracker
  for (i in 1:N){
    for (j in 1:N){
      if (migrations[i,j] == 1){
        
        out.names[current.parameter] = paste('m', i, '.', j, sep='')
        
        current.parameter <- current.parameter + 1
        
      }
    }
  }
  
  
  # Add permitted kills
  # If kill permitted, add its name to the name vector
  # and increment the position tracker
  for (i in 1:N){
    if (kills[i] == 1){
      
      out.names[current.parameter] = paste('k', i, sep='')
      
      current.parameter = current.parameter + 1
      
    }
  }
  
  
  # Add permitted replications
  # If replication permitted, add its name to the name vector
  # and increment the position tracker
  for (i in 1:N){
    if (replications[i] == 1){
      
      out.names[current.parameter] = paste('r', i, sep='')
      
      current.parameter = current.parameter + 1
      
    }
  }
  
  # Add permitted capacitances
  # If capacitance permitted, add its name to the name vector
  # and increment the position tracker
  if (capacitated){
    for (i in 1:N){
      if (capacities[i] == 1){
        
        out.names[current.parameter] = paste('c', i, sep='')
        
        current.parameter = current.parameter + 1
        
      }
    }
  }
  
  
  # Output named vector of zeros with names equal to names of parameters
  out <- vector('numeric', total.parameters)
  
  names(out) <- out.names
  
  return(out)
  
}


# Function to calculate the moments for a network with all possible inter-compartmental migrations, and intra-organ replication and killing

moment.sol.N.general <- 
function(N,t,par,M0,met="Higham08.b"){
  
  # Retrieves r_i from par
  
  replication <- function(i){
    par.name <- paste('r', i, sep='')
    par.val <- par[par.name]
    
    if (is.na(par.val)) {
      return(0)
    }
    
    else {
      return(par.val)
    }
  }
  
  # Retrieves k_i from par
  
  kill <- function(i){
    par.name <- paste('k', i, sep='')
    par.val <- par[par.name]
    
    if (is.na(par.val)) {return(0)}
    else {return(par.val)}
  }
  
  # Retrieves m_i1,i2 from par
  
  migration <- function(i1, i2){
    par.name <- paste('m', i1, '.', i2, sep='')
    par.val <- par[par.name]
    
    if (is.na(par.val)) {return(0)}
    else {return(par.val)}
    
  }
  
  
  # Calculates total migration out of a compartment
  
  total_migration <- function(i){
    
    out <- 0
    
    for (j in 1:N){
      if (i == j) next
      out <- out + migration(i,j)
    }
    
    return(out)
    
  }
  
  # First moment: M.1(t) = (E[n1],E[n2],...,E[nN]) is solution of M.1'(t) = A * M.1
  A <- matrix (0, N, N)
  
  for (i in 1:N){
    A[i,i] = replication(i) -  kill(i) - total_migration(i)
  }
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      A[i,j] = migration(j,i)
      A[j,i] = migration(i,j)
    }
  }
  
  # Solution for first moment
  M.1 <- expm::expAtv(A,M0[1:N],t)$eAtv
  
  
  # Second moment: M.2(t) = (Var(n1),Var(n2),..,Var(nN),Cov(n1,n2),..,Cov(nN-1,nN)) is solution of M.2' = B*M.1 + C*M.2
  # B and C matrices
  B <- matrix(0, N*(N+1)/2, N)
  C <- matrix(0, N*(N+1)/2, N*(N+1)/2)
  
  # Calculates index for (n_i1.n_i2)
  
  index <- function(i1, i2){
    return(i1*N - choose(i1,2) + i2 - i1)
  }
  
  # Var(ni)
  for (i in 1:N){
    
    B[i,i] <- total_migration(i) + kill(i) + replication(i)
    
    C[i,i] <- 2*A[i,i]
  }
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      
      m_ij <- migration(i,j)
      m_ji <- migration(j,i)
      
      B[i,j] <- m_ji
      B[j,i] <- m_ij
      
      ij = index(i,j)
      
      C[i, ij] <- 2*m_ji
      C[j, ij] <- 2*m_ij
      
    }
  }
  
  # Cov(ni, nj)
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      
      ij <- index(i,j)
      
      m_ij <- migration(i,j)
      m_ji <- migration(j,i)
      
      B[ij, i] <- -m_ij
      B[ij, j] <- -m_ji
      C[ij, i] <- m_ij
      C[ij, j] <- m_ji
      
      C[ij, ij] <- (A[i,i] + A[j,j])
      
      for (k in 1:N){
        if (k != i && k != j){
          
          # mins and maxes used here as the index function requires i1 < i2
          ik = index(min(i,k),max(i,k))
          jk = index(min(j,k),max(j,k))
          
          C[ij, ik] <- migration(k,j)
          C[ij, jk] <- migration(k,i)
        }
      }
    }
  }
  
  # Matrix exponential solution to integral
  
  D <- matrix(0, N*(N+1)/2, N)
  K <- matrix(0, N*(N+3)/2, N*(N+3)/2)
  
  # As a block matrix, K is:
  
  # C B
  # 0 A
  
  K[1:(N*(N+1)/2),1:(N*(N+1)/2)] <- C[1:(N*(N+1)/2),1:(N*(N+1)/2)]
  K[1:(N*(N+1)/2),(N*(N+1)/2+1):(N*(N+3)/2)] <- B[1:(N*(N+1)/2),1:N]
  K[(N*(N+1)/2+1):(N*(N+3)/2),(N*(N+1)/2+1):(N*(N+3)/2)] <- A[1:N,1:N]
  
  ## We take the exponential of K*t
  
  M <- expm::expm(t*K,m=met)
  
  ## D represents the matrix solution of the second moment equation
  
  D[1:(N*(N+1)/2),1:N] <- M[1:(N*(N+1)/2),(N*(N+1)/2+1):(N*(N+3)/2)]
  
  ## Solution for second moment
  
  M.2 <- D %*% M0[1:N]
  M.2 <- as.numeric(expm::expm(t*C,m = met) %*% M0[(N+1):(N*(N+3)/2)] + M.2)
  
  # Solutions
  
  return(c(M.1,M.2))
  
}





# Modified mde.est from SPEEDI.R package to include all possible inter-compartmental migrations

general.mde.est.noper <- 
  function(N, t, obs.moments, network, divergence, par.start, 
           par.fixed=NULL, init.moments.fun=NULL, init.moments.vals=NULL, 
           subs=rep(list(1:N), times=length(t)), constraint.fun = function(x) {
             return(min(x) < 0)
           },optim.routine=minqa::bobyqa, 
           combine="sum", capacitated=FALSE, ...){
    
    # Check if there are more parameters to be estimated than moments from which to estimate them
    if (length(obs.moments) < length(par.start)){
      print('Warning: There are fewer observed moments than parameters to estimate. Estimates may be very inaccurate.')
    }
    
    # Check if initial conditions are to be estimated
    
    calc.i0 <- FALSE
    if (!is.null(init.moments.vals)){
      calc.i0 <- TRUE
      
      # If initial conditions are to be estimated, but no names included,
      # add names to initial values vector
      
      if (is.null(names(init.moments.vals))){
        names(init.moments.vals) <- paste("init.val.",1:length(init.moments.vals), sep="")
      }
    }
    
    # Initialise parameter vector
    all.par <- network.params(network)
    
    # Update fixed parameters to their value
    all.par <- replace.par(all.par, par.fixed)
    
    # Stop if par.start is not named at all
    if(is.null(names(par.start))){
      stop("Please name entries in par.start")
    }
    # Stop if names of par.start are incorrect
    if (!all(names(par.start) %in% names(all.par))){
      stop('Ensure that names of par.start are correct.')
    }
    
    # Make sure obs.moments input is a matrix
    if (length(t)==1){
      obs.moments <- matrix(obs.moments, nrow=1)
    }
    
    # Pre-define matrix to contain moments of predicted system
    mom.i <- matrix(NA, nrow=length(t), ncol=14)
    div.i <- vector("numeric", length=length(t))
    
    
    # Initialise parameter vector (initial values + dynamic parameters)
    par.start <- c(init.moments.vals, par.start)
    
    # Calculate minimum divergence estimate
    mde.out <- do.call(optim.routine, list(par=par.start, fn = function(par.ext){
      
      names(par.ext) <- names(par.start)
      
      # If any parameters are guessed to be negative, return a prohibitively high value
      if(any(par.ext < 0)) return(1E100)
      
      if (calc.i0){
        # Get dynamic parameters (i0 only needed for initial conditions)
        par <- par.ext[-(1:length(init.moments.vals))]
        # Estimate proportion of inoculum
        i0 <- par.ext[1:length(init.moments.vals)]
      } else{
        par <- par.ext
        i0 <- NULL
      }
      # Update parameters being estimated
      par.i <- replace.par(all.par,par)
      
      # Calculate initial moments from initial conditions (e.g. Poisson from inoculum dose)
      init.mom.i <- init.moments.fun(i0)
      
      # Calculate moments + divergence measure for given parameters
      for (i in 1:length(t)){
        
        if (capacitated == TRUE){
          mom.i[i,] <- do.call(moment.sol.N.general.capacitated, list(N, t[i], par.i, init.mom.i))
        }
        
        else{
          mom.i[i,] <- do.call(moment.sol.N.general, list(N=N, t[i], par.i, init.mom.i))
          
          
        }
        
        
        div.i[i] <- div.measure(N = N, obs = obs.moments, pred = mom.i[i,], subs = subs[[i]], div = divergence)
      }
      
      # Combine divergence measures according to given function
      do.call(combine, as.list(div.i))
    })
    )
    
    # Retrieve parameter estimates (initial conditions + dynamic parameters)
    par.estimates <- c(mde.out$par, par.fixed)
    
    # If no names for parameters, add them
    
    if (is.null(names(par.estimates))){
      
      par.names <- vector('character', length(par.estimates))
      
      # If initial conditions are to be calculated, 
      # add names for initial conditions + dynamic parameters separately
      
      if (calc.i0){
        par.names[1:length(init.moments.vals)] <- names(init.moments.vals)
        par.names[-(1:length(init.moments.vals))] <- names(all.par)
      }
      
      # Otherwise just add names for dynamic parameters
      
      else {
        par.names <- names(c(par.start, par.fixed))
      }
      
      # Add names to parameter estimates
      
      names(par.estimates) <- par.names
      
    }
    
    # Retrieve observed divergence
    obs.div <-  mde.out$value
    
    
    # Return estimated parameters and observed divergence
    return(c(par.estimates, "obs.div"= obs.div))
    
  }




# Modified mde.est function in the SPEEDI.R package to combine moments in a dual phenotype 4-compartment model.
general.mde.est <- 
  function(N, t, obs.moments, network, divergence, par.start, 
           par.fixed=NULL, init.moments.fun=NULL, init.moments.vals=NULL, 
           subs=rep(list(1:N), times=length(t)), constraint.fun = function(x) {
             return(min(x) < 0)
           },optim.routine=minqa::bobyqa, 
           combine="sum", capacitated=FALSE, ...){
    
    # Check if there are more parameters to be estimated than moments from which to estimate them
    if (length(obs.moments) < length(par.start)){
      print('Warning: There are fewer observed moments than parameters to estimate. Estimates may be very inaccurate.')
    }
    
    # Check if initial conditions are to be estimated
    
    calc.i0 <- FALSE
    if (!is.null(init.moments.vals)){
      calc.i0 <- TRUE
      
      # If initial conditions are to be estimated, but no names included,
      # add names to initial values vector
      
      if (is.null(names(init.moments.vals))){
        names(init.moments.vals) <- paste("init.val.",1:length(init.moments.vals), sep="")
      }
    }
    
    # Initialise parameter vector
    all.par <- network.params(network)
    
    # Update fixed parameters to their value
    all.par <- replace.par(all.par, par.fixed)
    
    # Stop if par.start is not named at all
    if(is.null(names(par.start))){
      stop("Please name entries in par.start")
    }
    # Stop if names of par.start are incorrect
    if (!all(names(par.start) %in% names(all.par))){
      stop('Ensure that names of par.start are correct.')
    }
    
    # Make sure obs.moments input is a matrix
    if (length(t)==1){
      obs.moments <- matrix(obs.moments, nrow=1)
    }
    
    # Pre-define matrix to contain moments of predicted system
    mom.x <- matrix(nrow=length(t), ncol=44)
    mom.i <- matrix(NA, nrow=length(t), ncol=14)
    div.i <- vector("numeric", length=length(t))
    
    
    # Initialise parameter vector (initial values + dynamic parameters)
    par.start <- c(init.moments.vals, par.start)
    
    # Calculate minimum divergence estimate
    mde.out <- do.call(optim.routine, list(par=par.start, fn = function(par.ext){
      
      names(par.ext) <- names(par.start)
      
      # If any parameters are guessed to be negative, return a prohibitively high value
      if(any(par.ext < 0)) return(1E100)
      
      if (calc.i0){
        # Get dynamic parameters (i0 only needed for initial conditions)
        par <- par.ext[-(1:length(init.moments.vals))]
        # Estimate proportion of inoculum
        i0 <- par.ext[1:length(init.moments.vals)]
      } else{
        par <- par.ext
        i0 <- NULL
      }
      # Update parameters being estimated
      par.i <- replace.par(all.par,par)
      
      # Calculate initial moments from initial conditions (e.g. Poisson from inoculum dose)
      init.mom.i <- init.moments.fun(i0)
      
      # Calculate moments + divergence measure for given parameters
      for (i in 1:length(t)){
        
        if (capacitated == TRUE){
          mom.i[i,] <- do.call(moment.sol.N.general.capacitated, list(N, t[i], par.i, init.mom.i))
        }
        
        else{
          mom.x[i,] <- do.call(moment.sol.N.general, list(N=8, t[i], par.i, init.mom.i))
          # combine 
          # means
          mom.i[i,1] <- mom.x[i,1]+mom.x[i,5]
          mom.i[i,2] <- mom.x[i,2]+mom.x[i,6]
          mom.i[i,3] <- mom.x[i,3]+mom.x[i,7]
          mom.i[i,4] <- mom.x[i,4]+mom.x[i,8]
          
          #var
          mom.i[i,5] <- mom.x[i,9]+mom.x[i,13] +2*mom.x[i,20]
          mom.i[i,6] <- mom.x[i,10]+mom.x[i,14]+2*mom.x[i,27]
          mom.i[i,7] <- mom.x[i,11]+mom.x[i,15]+2*mom.x[i,33]
          mom.i[i,8] <- mom.x[i,12]+mom.x[i,16]+2*mom.x[i,38]
          
          # cov
          mom.i[i,9] <- mom.x[i,17]+mom.x[i,21]+mom.x[i,26]+mom.x[i,39]
          mom.i[i,10] <- mom.x[i,18]+mom.x[i,22]+mom.x[i,31]+mom.x[i,40]
          mom.i[i,11] <- mom.x[i,19]+mom.x[i,23]+mom.x[i,35]+mom.x[i,41]
          mom.i[i,12] <- mom.x[i,24]+mom.x[i,28]+mom.x[i,32]+mom.x[i,42]
          mom.i[i,13] <- mom.x[i,25]+mom.x[i,29]+mom.x[i,36]+mom.x[i,43]
          mom.i[i,14] <- mom.x[i,30]+mom.x[i,34]+mom.x[i,37]+mom.x[i,44]
          
        }
        
        
        div.i[i] <- div.measure(N = 4, obs = obs.moments, pred = mom.i[i,], subs = subs[[i]], div = divergence)
      }
      
      # Combine divergence measures according to given function
      do.call(combine, as.list(div.i))
    })
    )
    
    # Retrieve parameter estimates (initial conditions + dynamic parameters)
    par.estimates <- c(mde.out$par, par.fixed)
    
    # If no names for parameters, add them
    
    if (is.null(names(par.estimates))){
      
      par.names <- vector('character', length(par.estimates))
      
      # If initial conditions are to be calculated, 
      # add names for initial conditions + dynamic parameters separately
      
      if (calc.i0){
        par.names[1:length(init.moments.vals)] <- names(init.moments.vals)
        par.names[-(1:length(init.moments.vals))] <- names(all.par)
      }
      
      # Otherwise just add names for dynamic parameters
      
      else {
        par.names <- names(c(par.start, par.fixed))
      }
      
      # Add names to parameter estimates
      
      names(par.estimates) <- par.names
      
    }
    
    # Retrieve observed divergence
    obs.div <-  mde.out$value
    
    
    # Return estimated parameters and observed divergence
    return(c(par.estimates, "obs.div"= obs.div))
    
  }


######################################################################################
######################################################################################
load("Rossi.amp.wt.moments.correct.RData")
load("Rossi.cip.wt.moments.correct.RData")



######################################################################################
######################################################################################
######################################################################################

# Example of code for parameter inference using the SP model for ampicillin-treated mice during day 1 of treatment
x <- matrix(nrow=200, ncol=21)
# Define a general 4-compartmental network with killing and death in all compartments and all allowed inter-compartmental bacterial migration.
network <- all.poss(4)
t <- c(24)
N <- 4
obs.moments <- rbind(Rossi.amp.wt.moments.correct[2,])
divergence <- "KL"

# Initialisation of the optimisation algorithm from multiple points and selection of the parameter estimate with the smallest KL divergence
for (i in 1:200)
{
  set.seed(i)
  par.start <- c( "k1"=runif(1),
                  "k2"=runif(1), 
                  "k3"=runif(1), "k4"=runif(1), "r2"=runif(1),
                  "r3"=runif(1), "r4"=runif(1), "m4.2"=runif(1))
  
  par.fixed <- c( "m1.3"=0,"m1.2"=0, "m1.4"=0, 
                  "m2.1"=0, "m2.3"=0, "m2.4"=0, "m4.1"=0,
                  "m3.1"=0,  
                  "r1"=0,
                  "m4.3"=0, "m3.2"=0, "m3.4"=0
  )
  
  constraint.fun <- function(x){ return(x["r2"]>1.5 ||
                                          x["r3"]>1.5 ||
                                          x["r4"]>1.5 ||min(x)<0)}
  
  
  init.moments.fun <- function(x){return(c(Rossi.amp.wt.moments.correct[1,]))}
  subs <- list(c(c(1,2,3,4)))
  x[i,] <- general.mde.est.noper(N=N, t=t, obs.moments=obs.moments, network=network, par.start=par.start,
                                 init.moments.fun = init.moments.fun,divergence=divergence, subs=subs,
                                 constraint.fun=constraint.fun, par.fixed=par.fixed,
                                 optim.routine = "optim")
  print(x[i,])
}

est3to4.wt.amp <- x[which.min(x[,21]),1:21]





# Simulations using the same initial conditions and best parameter estimates to obtain the bootstapped samples
est3to4.wt.amp.simulations <- list()
o <- matrix(nrow=80, ncol=5)

for (i in 1:500)
{
  transitios <-  list(
    c(B=-1, L=1, M=0, S=0), #mBL
    c(B=-1, L=0, M=1, S=0), #mBL
    c(B=-1, L=0, M=0, S=1), #mBS
    
    
    c(B=1, L=-1, M=0, S=0), #mLB
    c(B=1, L=0, M=-1, S=0), #mMB
    c(B=1, L=0, M=0, S=-1), #mSB
    c(B=0, L=1, M=0, S=-1), #mSL
    
    
    c(B=-1, L=0, M=0, S=0), #kB
    c(B=0, L=-1, M=0, S=0), #kL
    c(B=0, L=0, M=-1, S=0), #kM
    c(B=0, L=0, M=0, S=-1), #kS
    
    c(B=0, L=1, M=0, S=0), #rL
    c(B=0, L=0, M=1, S=0), #rM
    c(B=0, L=0, M=0, S=1)) #rS
  
  
  
  my.rates <- function(x, params, t) {
    return(c(params$mBL*x["B"],
             params$mBM*x["B"],
             params$mBS*x["B"],
             
             params$mLB*x["L"],
             params$mMB*x["M"],
             params$mSB*x["S"],
             params$mSL*x["S"],
             
             params$kB*x["B"],
             params$kL*x["L"],
             params$kM*x["M"],
             params$kS*x["S"],
             
             
             params$rL*x["L"],
             params$rM*x["M"],
             params$rS*x["S"])
    )
  }
  params = list(
    mBL=0,
    mBM=0,
    mBS=0,
    mLB=0,
    mMB=0,
    mSB=0,
    mSL=est3to4.wt.amp[8],
    kB=est3to4.wt.amp[1],
    kL=est3to4.wt.amp[2],
    kM=est3to4.wt.amp[3],
    kS=est3to4.wt.amp[4],
    rL=est3to4.wt.amp[5],
    rM=est3to4.wt.amp[6],
    rS=est3to4.wt.amp[7])
  
  for (j in 1:80)
  {
    colnames(Rossi.wt.amp_3) <- c("B", "L", "M", "S")
    start.conditions <- round(Rossi.wt.amp_3[sample(nrow(Rossi.wt.amp_3),size=1,replace=T),])
   
    o[j,] <- tail(ssa.adaptivetau(init.values=start.conditions, transitions=transitios,
                                  rateFunc=my.rates, params=params, tf=24), n=1)
  }
  est3to4.wt.amp.simulations[[i]] <- o[,2:5]
}

# Calculation of moments for bootstrapped samples
est3to4.wt.amp.simulations.moments <- matrix(nrow=500, ncol=14)
for (i in 1:500)
{
  est3to4.wt.amp.simulations.moments[i,] <- 
    SPEEDI.R:::network.moments(est3to4.wt.amp.simulations[[i]])
}


################################################################################
################################################################################
################################################################################
# Convert the moments for the 4-compartment SP model into the moments for the 8-compartment DP model
load("Rossi.wt.amp_3.RData")
Rossi.wt.amp_3.het <- cbind(Rossi.wt.amp_3, rep(0,80), rep(0,80), 
                            rep(0,80), rep(0,80))
Rossi.wt.amp_3.het.moments <- SPEEDI.R:::network.moments(Rossi.wt.amp_3.het)



# Example of code for parameter inference using the DP model for ampicillin-treated mice during day 1 of treatment
x <- matrix(nrow=200, ncol=73)
# Define a general 8-compartmental network (4 growing, antibiotic-sensitive populations and 4 non-growing, antibiotic recalcitrant populations) with killing and death in all compartments and all allowed inter-compartmental bacterial migration.
network <- all.poss(8)
t <- c(24)
N <- 8
obs.moments <- rbind(Rossi.amp.wt.moments.correct[2,])
divergence <- "KL"
for (i in 1:200)
{
  set.seed(i)
  par.start <- c( 
    "k1"=runif(1),
    "k2"=runif(1) ,
    "k3"=runif(1),
    "k4"= runif(1),
    "r2"=runif(1),
    "r3"=runif(1),
    "r4"=runif(1), 
    "m4.2"=runif(1)
  )
  
  par.fixed <- c( "m1.3"=0,  "m1.6"=0, "m1.7"=0, "m1.8"=0,     "m2.1"=0,
                  "m2.3"=0, "m2.7"=0, "m2.8"=0, 
                  "m3.1"=0, "m3.2"=0, "m3.4"=0, "m3.5"=0, "m3.6"=0, "m3.8"=0,
                  "m4.1"=0,"m4.3"=0, "m4.5"=0, "m4.6"=0, "m4.7"=0,
                  "m5.1"=0, "m5.2"=0, "m5.3"=0, "m5.4"=0, "m5.6"=0, "m5.7"=0, "m5.8"=0,
                  "m6.1"=0, "m6.2"=0, "m6.3"=0, "m6.4"=0, "m6.5"=0, "m6.7"=0, "m6.8"=0,
                  "m7.1"=0, "m7.2"=0,"m7.4"=0, "m7.5"=0, "m7.6"=0, "m7.8"=0,
                  "m8.1"=0, "m8.2"=0, "m8.3"=0,  "m8.5"=0,  
                  "k5"=0, "k6"=0, "k7"=0,
                  "r1"=0,  "r5"=0, "r6"=0, "r7"=0, "r8"=0,
                  
                  "m1.5"=0, 
                  "m7.3"=0, "m2.5"=0, 
                  "m8.4"=0, 
                  "m8.7"=0, 
                  
                  "k8"=0, "m2.4"=0,     
                  "m1.2"=0, "m1.4"=0 , "m2.6"=0,
                  "m8.6"=0, "m3.7"=0, "m4.8"=0
                  
                  
                  
                  
  )
  
  
  constraint.fun <- function(x){ return(x["r2"]>1.5 ||
                                          x["r3"]>1.5 ||
                                          x["r4"]>1.5 ||min(x)<0)}
  
  
  init.moments.fun <- function(x){return(c(Rossi.wt.amp_3.het.moments))}
  subs <- list(c(1:4))
  x[i,] <- general.mde.est(N=N, t=t, obs.moments=obs.moments, network=network, 
                           par.start=par.start,
                           init.moments.fun = init.moments.fun,
                           divergence=divergence, subs=subs,
                           constraint.fun=constraint.fun, par.fixed=par.fixed,
                           optim.routine = "optim", capacitated = F)
  print(x[i,])
}

est3to4.wt.amp.hetno3.2 <- x[which.min(x[,73]),1:73]



# Simulations using the same initial conditions and best parameter estimates to obtain the bootstapped samples
est3to4.wt.amp.simulations.het <- list()
o <- matrix(nrow=80, ncol=9)

for (i in 1:500)
{
  transitios <-  list(
    
    c(B=-1, L=0, M=0, S=0, Bp=0, Lp=0, Mp=0, Sp=0), #kB
    c(B=0, L=-1, M=0, S=0, Bp=0, Lp=0, Mp=0, Sp=0), #kL
    c(B=0, L=0, M=-1, S=0, Bp=0, Lp=0, Mp=0, Sp=0), #kM
    c(B=0, L=0, M=0, S=-1, Bp=0, Lp=0, Mp=0, Sp=0), #kS
    c(B=0, L=1, M=0, S=0, Bp=0, Lp=0, Mp=0, Sp=0), #rL
    
    c(B=0, L=0, M=1, S=0, Bp=0, Lp=0, Mp=0, Sp=0), #rM
    c(B=0, L=0, M=0, S=1, Bp=0, Lp=0, Mp=0, Sp=0), #rS
    
    c(B=0, L=1, M=0, S=-1, Bp=0, Lp=0, Mp=0, Sp=0), #mSL
    c(B=0, L=0, M=0, S=0, Bp=0, Lp=1, Mp=0, Sp=-1), #mSL
    c(B=0, L=0, M=0, S=-1, Bp=0, Lp=0, Mp=0, Sp=1), #StoS'
    c(B=0, L=0, M=-1, S=0, Bp=0, Lp=0, Mp=1, Sp=0),#MtoM'
    c(B=0, L=0, M=0, S=0, Bp=0, Lp=0, Mp=0, Sp=-1)) #kSp
  
  
  params = list(
    
    kB=est3to4.wt.amp.hetno3.2 [1],
    kL=est3to4.wt.amp.hetno3.2 [2],
    kM=est3to4.wt.amp.hetno3.2 [3],
    kS=est3to4.wt.amp.hetno3.2 [4],
    rL=est3to4.wt.amp.hetno3.2 [5],
    rM=est3to4.wt.amp.hetno3.2 [6],
    rS=est3to4.wt.amp.hetno3.2 [7],
    mSL=est3to4.wt.amp.hetno3.2 [8],
    mSpLp=0,
    StoSp=0,
    MtoMp=0,
    kSp=0
  )
  
  
  my.rates <- function(x, params, t) {
    return(c(params$kB*x["B"],
             params$kL*x["L"],
             params$kM*x["M"],
             params$kS*x["S"],
             params$rL*x["L"],
             params$rM*x["M"],
             params$rS*x["S"],
             params$mSL*x["S"],
             params$mSpLp*x["Sp"],
             params$StoSp*x["S"],
             params$MtoMp*x["M"],
             params$kSp*x["Sp"]))
  }
  
  
  for (j in 1:80)
  {
    colnames(Rossi.wt.amp_3.het) <- c("B", "L", "M", "S", "Bp", "Lp", "Mp", "Sp")
    start.conditions <- round(Rossi.wt.amp_3.het[sample(nrow(Rossi.wt.amp_3.het),size=1,replace=T),])
  
    o[j,] <- tail(ssa.adaptivetau(init.values=start.conditions, transitions=transitios,
                                  rateFunc=my.rates, params=params, tf=24), n=1)
  }
  est3to4.wt.amp.simulations.het[[i]] <- o[,2:9]
}


# Calculation of moments for bootstrapped samples (8-compartments)
est3to4.wt.amp.simulations.het.moments <- matrix(nrow=500, ncol=44)
for (i in 1:500)
{
  est3to4.wt.amp.simulations.het.moments[i,] <- 
    SPEEDI.R:::network.moments(est3to4.wt.amp.simulations.het[[i]])
}

# Combination of moments for bootstrapped samples from 8 to 4 compartments
est3to4.wt.amp.simulations.het.moments.comb <- matrix(nrow=500, ncol=14)
for (i in 1:500)
{
  # combine 
  # means
  est3to4.wt.amp.simulations.het.moments.comb[i,1] <- est3to4.wt.amp.simulations.het.moments[i,1]+est3to4.wt.amp.simulations.het.moments[i,5]
  est3to4.wt.amp.simulations.het.moments.comb[i,2] <- est3to4.wt.amp.simulations.het.moments[i,2]+est3to4.wt.amp.simulations.het.moments[i,6]
  est3to4.wt.amp.simulations.het.moments.comb[i,3] <- est3to4.wt.amp.simulations.het.moments[i,3]+est3to4.wt.amp.simulations.het.moments[i,7]
  est3to4.wt.amp.simulations.het.moments.comb[i,4] <- est3to4.wt.amp.simulations.het.moments[i,4]+est3to4.wt.amp.simulations.het.moments[i,8]
  
  #var
  est3to4.wt.amp.simulations.het.moments.comb[i,5] <- est3to4.wt.amp.simulations.het.moments[i,9]+est3to4.wt.amp.simulations.het.moments[i,13] +2*est3to4.wt.amp.simulations.het.moments[i,20]
  est3to4.wt.amp.simulations.het.moments.comb[i,6] <- est3to4.wt.amp.simulations.het.moments[i,10]+est3to4.wt.amp.simulations.het.moments[i,14]+2*est3to4.wt.amp.simulations.het.moments[i,27]
  est3to4.wt.amp.simulations.het.moments.comb[i,7] <- est3to4.wt.amp.simulations.het.moments[i,11]+est3to4.wt.amp.simulations.het.moments[i,15]+2*est3to4.wt.amp.simulations.het.moments[i,33]
  est3to4.wt.amp.simulations.het.moments.comb[i,8] <- est3to4.wt.amp.simulations.het.moments[i,12]+est3to4.wt.amp.simulations.het.moments[i,16]+2*est3to4.wt.amp.simulations.het.moments[i,38]
  
  # cov
  est3to4.wt.amp.simulations.het.moments.comb[i,9] <- est3to4.wt.amp.simulations.het.moments[i,17]+est3to4.wt.amp.simulations.het.moments[i,21]+est3to4.wt.amp.simulations.het.moments[i,26]+est3to4.wt.amp.simulations.het.moments[i,39]
  est3to4.wt.amp.simulations.het.moments.comb[i,10] <- est3to4.wt.amp.simulations.het.moments[i,18]+est3to4.wt.amp.simulations.het.moments[i,22]+est3to4.wt.amp.simulations.het.moments[i,31]+est3to4.wt.amp.simulations.het.moments[i,40]
  est3to4.wt.amp.simulations.het.moments.comb[i,11] <- est3to4.wt.amp.simulations.het.moments[i,19]+est3to4.wt.amp.simulations.het.moments[i,23]+est3to4.wt.amp.simulations.het.moments[i,35]+est3to4.wt.amp.simulations.het.moments[i,41]
  est3to4.wt.amp.simulations.het.moments.comb[i,12] <- est3to4.wt.amp.simulations.het.moments[i,24]+est3to4.wt.amp.simulations.het.moments[i,28]+est3to4.wt.amp.simulations.het.moments[i,32]+est3to4.wt.amp.simulations.het.moments[i,42]
  est3to4.wt.amp.simulations.het.moments.comb[i,13] <- est3to4.wt.amp.simulations.het.moments[i,25]+est3to4.wt.amp.simulations.het.moments[i,29]+est3to4.wt.amp.simulations.het.moments[i,36]+est3to4.wt.amp.simulations.het.moments[i,43]
  est3to4.wt.amp.simulations.het.moments.comb[i,14] <- est3to4.wt.amp.simulations.het.moments[i,30]+est3to4.wt.amp.simulations.het.moments[i,34]+est3to4.wt.amp.simulations.het.moments[i,37]+est3to4.wt.amp.simulations.het.moments[i,44]
  
}
