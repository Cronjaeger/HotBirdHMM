# Implements sampling and ML-estimation for a HMM
#
# © Mathias C. Cronjäger 12/9/2017
# Developped as part of a project @BIRC Aarhus w. Asger Hobolth, Kasper M. Terkelsen and Jonas Berglund




#--------------------------
# Step 1: Setting the stage
#--------------------------

#-------------------
#  1.1: Import code
#-------------------
source('./R/from-hobolth/TransMat.R')
# Note: Since I didn't write TransMat.R myself, I have not put it online.
# Ask Asger if you want the code (maybie I'll add it later should he permit)


#---------------------------------
#  1.2: Define Auxiliary functions
#---------------------------------


#' Generate discretization timepoints
#'
#' Compute the sequence of discretization timepoints (cf. Asger's note)
#'
#' @param M The number of intervals in the discretization; should be a positive
#'   integer
#' @return A vector of upper limits of each interval. last element will be
#'   \code{Inf}
#' @export
d_seq <- function(M){

    #Validate input
  if(!is.numeric(M) || M < 1) stop('Invalid value of M')

  return(-log(1 - 1:(M) / M))
}

#' Compute emission probabilities for HMM
#'
#' @param M the number of states of the hidden variable X
#' @param theta the rate of mutation between two lineages
#' @return A vector mu s.t. mu[i] = P(Y = 1 | X = i) for i=1...M
#' @export
compute_mu <- function(M,theta){

  #validate input
  if(!is.numeric(M) || M < 1) stop('Invalid value of M')
  if(!is.numeric(theta) || theta < 0) stop('Invalid theta argument!')

  ds_padded <- c(0, d_seq(M)) # d_sequence with 0-padding
  d <- ds_padded[1:(length(ds_padded)-1)] # the d_m-1 sequence (cf. asgers Note)
  D <- ds_padded[2:(length(ds_padded))]   # the d_m sequence

  a <- exp(-theta * d)
  b <- 1 - exp((1+theta) * (d - D))
  c <- (1 + theta) * (1 - exp(d - D))
  result <- 1 - (a * b / c)
  #Note@self(1): If d_m = -ln(1 - m/M) for all m=1...M, exp-terms might simplify.
  #Note@self(2): Think about if there is a more numerically stable order of operations.

  return(result)
}

#' Compute hidden state transition kernel
#'
#' A wrapper function for calling A. Hobolths code for computin the hiden state
#' rate matrix by way of matrix exponentiation (as outlined in his note).
#'
#' @param M the number of hissden states
#' @param rho The recombination rate between two sites
#' @return An M-by-M matrix P, s.t. \code{P[i,j]} = P(X_t = j | X_{t-1} = i)
#' @export
compute_Px <- function(M, rho){
  # Compute the transition matrix of the discretized tree-hight process X
  # Arguments:
  #   M: number of intervals in discretization
  #   rho: recombination rate

  if(!is.numeric(rho) || rho < 0) stop('Invalid value of rho!')

  return(TransMat(d_seq(M), rho))
  }

#' Count changepoints
#'
#' Scan along a sequence and count the number of times an entry does not match
#' the previous entry
#'
#' @param vec A vector to be scanned
#' @return a positive integer
#' @export
count_changepoints <- function(vec){
  return(sum(diff(vec) != 0))
}

#' Count runs
#'
#' Scan along a sequence; return the number of intervals wherein ajacent entries
#' are identical
#'
#' @param vec A vector to be scanned
#' @return The number of runs
#' @export
count_runs <- function(vec){
  return(count_changepoints(vec) + 1)
}



#-------------------------
# Step 2: Define Samplers
#-------------------------


#' Sample hidden states
#'
#' Sample X consequtive states from the hidden markov model
#'
#' @param sample_X0 A prior distribution on X0; sample_X0() should return an
#'   integer in the range 1:(dim(Px)[1]); e.g. “function() sample(M,1)” for a
#'   uniform distribution on 1:M or “function() i” for a point mass at i.
#' @param L the length of the desired sequence of samples
#' @param Px The transition rates of the markov chain X to be sampled from.
#'
#' @return A vector of samples from X of length \code{L}.
#'
#' @export
sample_X <- function(sample_X0, L, Px){

  #validate input
  if(!is.numeric(L) || L < 1) stop('Invalid value of L')

  #initalize X
  X <- vector(mode = 'integer', length = L)

  M <- dim(Px)[1]

  update <- function(x_prev){
    x_new <- sample(x=M, size = 1, prob = Px[x_prev,])
    return(x_new)
  }

  X[1] <- sample_X0()

  i <- 2
  while(i <= L){
    X[i] <- update(X[i-1])
    i <- i+1
  }

  return(X)
}

#' Sample emissions
#'
#' Samples a sequence of binary emissions (Y[i] = 1 is interpreted as,
#' “Heterozygous sites at position i”; whereas Y[i] = 0 is interpreted as
#' signifying homozugosity.)
#'
#' @param X a vector of integers encoding the states of the hidden variables
#' @param theta the rate of mutation between two sequences
#' @param M the number of states of the hidden variables X
#' @return A vector of the same length as \code{X} s.t. Y[i] = 1 with
#'   probability \code{mu[X[i]]}.
#' @export
sample_YgivenX <- function(X,theta,M){
  #Returns a Y-sequence (binary) given an X sequence (in 1:M)
  # corresponds to Y|(X=1)

  #Validate input
  if(M < max(X)) stop('stop: M < max(X)!')

  #initialize
  Y <- vector(mode = 'integer',length=length(X))
  mu_precomputed <- compute_mu(M,theta)

  sampleYi <- function(x){
    p <- mu_precomputed[x]
    return(sample(2,size = 1,prob = c(1-p,p)) - 1)
  }

  Y <- sapply(X,sampleYi)

}

#' Sample full model
#'
#' Samples both hidden states and emissions
#'
#' @param rho, rathe of recombination between ajacent sites; affects transition
#'   rate of hidden states.
#' @param theta Rate of mutation; affects distribution of Y|X.
#' @param M number of hidden states of X.
#' @param L desired length of sampled sequences.
#' @param sample_x0 a function which generates samples the initial state of the
#'   hidden variables X.
#' @return a list with names “X” and “Y”
#'   X is a vector of length \code{L} encoding the hidden states
#'   Y is a vector of length \code{L} encoding the emissions.
#' @export
sample_XandY <- function(rho,theta,M,L,sample_x0 = NULL){

  #validate prior
  if(is.null(sample_x0)) sample_x0 <- function() sample(M,1)
  if(is.numeric(sample_x0)){
    x0 <- sample_x0
    sample_x0() <- function() x0
  }
  if(!is.function(sample_x0)) stop('sample_x0 must be either integer, NULL, or a function')

  Px <- compute_Px(M,rho)

  X <- sample_X(sample_x0,L,Px)

  Y <- sample_YgivenX(X,theta,M)

  return(list(X=X,Y=Y))
  }

#------------------------------
# Step 3: Implement ML estimator
#------------------------------

#' Watterson's estimator
#'
#' Compute Watterson's estimator of the site-wise mutation rate
#'
#' @param s number of segregating sites observed.
#' @param n number of sequences sampled.
#' @param L number of sites (segregating + non-segregating) observed in sample.
#' @return A non-negative floating point number
#' @export
theta_watterson <- function(s,n,L){
  # Compute Watterson's estimator of the mutation rate
  # Arguments:
  #  s: number of segregating sites
  #  n: number of sequences
  #  L: sequence length (i.e. total number of sites)
  return((s / L) / sum(1 / 1:(n-1)))
}


#' Forward algorithm
#'
#' Computes alpha_T[i] = P(X_T = i, Y_1:T = y_1:T) for i=1...M, where M is the
#' number of hidden states, using the so-called forward-algorithm. This
#' implementation is less numerically stable than
#' \code{forward_algorithm_logspace}, which is otherwise identical, but works in
#' logspace instead.
#'
#' @param Y sequence of observed emissions
#' @param trnasition_probs an M-by-M matrix of transition probabilities
#' @param emission_probs a vector of emission probabilities. The ith entry is
#'   interpreted as the probability of Y=1 given X=i.
#' @param prior a vector encoding the a prior distribution on X0; not X[1], but
#' the stete preceding it, i.e. there is no observed emission of the variable X0
#' @return a vector alpha with entries P(X_T = i, Y_1:T = y_1:T) for i=1...M.
#' @export
forward_algorithm <- function(Y, transition_probs, emission_probs, prior){
  alpha0 <- prior
  mu1 <- emission_probs
  mu0 <- 1 - emission_probs

  alpha_update <- function(alpha, y){
    if(y==1) return(c((alpha %*% transition_probs) * mu1))
    if(y==0) return(c((alpha %*% transition_probs) * mu1))
    else stop('y neither 0 or 1; this should not happen.')
  }

  #compute alpha_T = (P(X_T = i, Y_1:T = y_1:T), i = 1...M)
  alpha <- Reduce(f = alpha_update, x = Y, init = alpha0)

  alpha <- c(alpha) #we're done with matrix arithmetics; turn into a vector.

  return(alpha)
}

#' Forward algorithm
#'
#' Computes ln( P(X_T = i, Y_1:T = y_1:T) ) for i=1...M, where M is the
#' number of hidden states, using the so-called forward-algorithm.
#' (cf. e.g. \url{https://en.wikipedia.org/wiki/Forward_algorithm})
#' This implementation is more stable than \code{algorithm_logspace}, which is
#' otherwise identical, but which does not work in logspace.
#'
#' @param Y sequence of observed emissions
#' @param trnasition_probs an M-by-M matrix of transition probabilities
#' @param emission_probs a vector of emission probabilities. The ith entry is
#'   interpreted as the probability of Y=1 given X=i.
#' @param prior a vector encoding the a prior distribution on X0; not X[1], but
#' the stete preceding it, i.e. there is no observed emission of the variable X0
#' @return a vector alpha with entries ln( P(X_T = i, Y_1:T = y_1:T) )
#'   for i=1...M.
#' @export
forward_algorithm_logspace <- function(Y, transition_probs, emission_probs, prior){

  ltrans <- log(transition_probs)
  lemiss1 <- log(emission_probs)
  lemiss0 <- log(1 - emission_probs)

  lalpha0 <- log(prior)

  M <- length(emission_probs)

  #we define the update-step
  lalpha_update <- function(lalpha,y){
    lalpha_prev <- lalpha
    j <- 1
    while(j <= M){

      #compute log(p_1j * a_1), ..., log(p_Mj * a_M)
      log_summands <-ltrans[,j]+lalpha_prev

      # compute sum_i(p_ij * a_i ) = P(X_t = j | X_(t-1) = i, Y_{1:t-1}=y_{1:t-1})
      prob_sum <- sum(exp(log_summands))

      # compute log(a) = log(p(Y = y | X=j) * log(prob_sum))
      lalpha[j] <- ifelse(y == 1,
                          { lemiss1[j] + log(prob_sum) },
                          { lemiss0[j] + log(prob_sum)})

      j <- j+1
    }
    return(lalpha)
  }

  result <- Reduce(f = lalpha_update, x = Y, init = lalpha0)

  return(result)
}

#' Marginal likelihood of Y
#'
#' Marginalizes out the hidden variables X using the forward algorithm to
#' compute the likelihood of an observed sequence of emissions.
#'
#' @param Y observed sequence of emissions (1: heterozygous site;
#'   0: homozygous site).
#' @param M size of hidden statespace; an integer.
#' @param rho the recombinartion rate between two ajacent sites.
#' @param theta mutation rate between two sequences at any given position.
#' @param Px transition kernel of hidden variables. May be passed instead of
#'    rho.
#' @param mu a vector of emission probabilities. The ith entry is
#'   interpreted as the probability of Y=1 given X=i. May be passed instead of
#'   theta
#' @param x0_prior initial distribution of the hidden markov chain. Defaults to
#'   a uniform distribution.
#' @return a single number, corresponding to L(theta,rho|Y)
#' @export
likelihood_Y <- function(Y, M, rho= NULL, theta = NULL, Px = NULL, mu=NULL, prior_x0 = NULL){
  return(exp(log_likelihood_Y(Y,M, rho, theta, Px, mu, prior_x0)))
}

#' Marginal log-likelihood of Y
#'
#' Marginalizes out the hidden variables X using the forward algorithm to
#' compute the log-likelihood of an observed sequence of emissions.
#'
#' @param Y observed sequence of emissions (1: heterozygous site;
#'   0: homozygous site).
#' @param M size of hidden statespace; an integer.
#' @param rho the recombinartion rate between two ajacent sites.
#' @param theta mutation rate between two sequences at any given position.
#' @param Px transition kernel of hidden variables. May be passed instead of
#'    rho.
#' @param mu a vector of emission probabilities. The ith entry is
#'   interpreted as the probability of Y=1 given X=i. May be passed instead of
#'   theta
#' @param x0_prior initial distribution of the hidden markov chain. Defaults to
#'   a uniform distribution.
#' @return a single number, corresponding to l(theta,rho|Y)
#' @export
log_likelihood_Y <- function(Y, M, rho= NULL, theta = NULL, Px = NULL, mu=NULL, prior_x0 = NULL){

  if(is.null(Px)){
    if(is.null(rho)) stop('You must pass either Px or rho as an argument')
    Px <- compute_Px(M,rho)
  } else {
    if(!is.null(rho)) warning('Px passed as argument; ignoring rho-argument')
  }

  if(is.null(mu)){
    if(is.null(theta)) stop('You must pass either mu or theta as an argument')
    mu <- compute_mu(M,theta)
  } else {
    if(!is.null(theta)) warning('mu passed as argument; ignoring theta-argument')
  }

  if(is.null(prior_x0)) prior_x0 <- rep(1,M)/M

  lalpha <- forward_algorithm_logspace(Y,Px,mu,prior_x0) # log alpha_values
  #(lalpha[i]: log(P(X_T = i, Y_1 = y_1, ... , Y_T = y_T)))

  return(log(sum(exp(lalpha))))
}

#' Joint likelihood of (X,Y)
#'
#' Compute the "full" likelihood of both emissions (Y) and hidden states (X).
#'
#' @param X Sequence of hidden states
#' @param Y sequence of emissions (1: heterozygous site; 0: homozygous site).
#' @param M size of hidden statespace; an integer.
#' @param rho the recombinartion rate between two ajacent sites.
#' @param theta mutation rate between two sequences at any given position.
#' @param Px transition kernel of hidden variables. May be passed instead of
#'    rho.
#' @param mu a vector of emission probabilities. The ith entry is
#'   interpreted as the probability of Y=1 given X=i. May be passed instead of
#'   theta
#' @param x0_prior initial distribution of the hidden markov chain. Defaults to
#'   a uniform distribution.
#' @return a single number, corresponding to L(theta,rho|Y)
log_likelihood_XY <- function(X, Y, M, rho= NULL, theta = NULL, Px = NULL, mu=NULL, prior_x0 = NULL){

  if (is.null(Px)){
    if(is.null(rho)) stop('You must pass either Px or rho as an argument')
    Px <- compute_Px(M,rho)
  } else {
    if (!is.null(rho)) warning('Px passed as argument; ignoring rho-argument')
  }

  if(is.null(mu)){
    if(is.null(theta)) stop('You must pass either mu or theta as an argument')
    mu <- compute_mu(M,theta)
  } else {
    if(!is.null(theta)) warning('mu passed as argument; ignoring theta-argument')
  }

  if(is.null(prior_x0)) prior_x0 <- rep(1,M) / M

  lPx <- log(Px)
  lmu1 <- log(mu)
  lmu0 <- log(1-mu)
  L <- length(X)

  acc <- log(sum(prior_x0 * Px[,X[1]])) + ifelse(Y[1] == 1, lmu1[X[1]], lmu1[X[1]] )
  i <- 2
  while(i <= L){
    stop('not implemented')
    #TODO: FINISH THIS
    i <- i+1
  }

}

#' Maximum likelihood of recombination rate
#'
#' Approximate the maximum likelihood estimate of the recombination rate.
#' This is done using nlm() for newtonian maximization of the likelihood
#' function.
#' @param Y A binary sequence classifying all observed sites as either
#'   homozygous (encoded Y[i] == 0) or heterozygous (encoded Y[i] == 1)
#' @param M the number of hidden states to be used in the HMM
#' @param rho0 initial guess used as starting point of gradiant-algorithm.
#' @param theta The mutation rate to be used (defaults to Watterson's
#'   estimator).
#' @return a maximum likelihood estimate of the recombination rate between two
#'   sites.
#' @export
compute_ML_rho <- function(Y,M,theta = NULL,rho0 = 1e-8){

  if(is.null(theta)) theta <- theta_watterson(sum(Y), 2, length(Y))

  mu <- compute_mu(M,theta)

  objective <- function(rho){
    -log_likelihood_Y(Y,M,rho = rho, mu = mu)
  }

  guess_scale <- 1/M * 1/length(Y)

  solver_result <- nlm(objective, rho0, typsize = rho0, fscale = guess_scale, steptol = 1e-6, stepmax = 1e3,print.level = 1)
  #solver_result <- nlm(objective,p = rho0)


  if(solver_result$code == 3) warning('ast global step failed to locate a point lower than estimate. Either estimate is an approximate local minimum of the function, or steptol is too small (i.e. we should consider allowing larger steps ')

  if(solver_result$code == 4) warning('Iteration limit exceeded (computing compute_ML_rho); take result with a grain of salt')

  if(solver_result$code == 5) warning('Maximum stepsize exceeded 5 consecutive times; consider either increasing stepsize or starting with a better initial guess')

  return(solver_result$estimate)
}
