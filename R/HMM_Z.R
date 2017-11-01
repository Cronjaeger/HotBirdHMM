# Implements recombination estimation based on “incompatability-distance”
#
# (C) Mathias C. Cronjäger 6/10/2017
# Developped as part of a project @BIRC Aarhus w. Asger Hobolth, Kasper M.
# Terkelsen and Jonas Berglund
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#' Perform 4 (or 3) gamete test
#'
#' Returns TRUE if vec1 and vec2 pass the four gamete test
#' (cf. e.g. Willson 1965 https://doi.org/10.2307/2411550)
#'
#' If the ancestral type is known, and passed as an argument, the three-gamete
#' test will be perfomed.
#'
#' @param vec1 the first vector to be examined
#' @param vec2 the seccond vector to be examined
#' @param ancestral_types (optional) indication of the ancestral type of each vector. Must be either NULL (indicating that ancestor is unknown and that a 4 gammete test should be performed), or a vector of length 2: the first entry being the ancestral type of vec1 and the seccond being the ancestral type of vec2.
#' @return a boolean indicating if two characters are compatible or not.
#' @export
gamete_test <- function(vec1, vec2, ancestral_types = NULL){
  if (length(vec1) != length(vec2)) stop('Arguments must have same length')
  if (length(unique(vec1)) > 2 | length(unique(vec1)) >2) warning('Four gammete test only valid for input vectors that each have ony two types!')

  #convert inputs to boolean vectors with 0 in first entry
  if(length(ancestral_types) == 0 && is.null(ancestral_types)){
    v1 <- vec1 != vec1[1]
    v2 <- vec2 != vec2[1]
  } else {
    if(length(ancestral_types) != 2) stop('Ancestor-type must be Null or a vector of length 2')
    v1 <- vec1 != ancestral_types[1]
    v2 <- vec2 != ancestral_types[2]
  }

  a = any( v1 &  v2) # does 00 occur?
  b = any( v1 & !v2) # does 01 occur?
  c = any(!v1 &  v2) # does 10 occur?
  d = any(!v1 & !v2) # does 11 occur?

  return(!(a & b & c & d))
}

#' Test if a character is informative
#'
#' A character is considered informative if the derived type has greater
#' frequency than 1; if the ancestral state is unknown this is equivalent to the
#' at least two nucleotides that appearing more than once.
#'
#' @param vec the observed character
#' @param ancestral_state the state of the most recent common ancestor, if known
#' defalults to NULL which signifies an unknown ancestral state
#' @return a boolean value: TRUE if the character is informative; oftherwise FALSE
is_informative_char <- function(vec, ancestral_state = NULL){

  if (!is.null(ancestral_state)){
    vec <- c(ancestral_state,vec)
  }

  char_states <- unique(vec)

  sorted_counts <- sort(sapply(char_states, function(x) sum(vec == x)), decreasing = TRUE)

  return (length(sorted_counts) > 1 && sorted_counts[2] > 1)

  # #The nucleotide of the first haplotype is our reference nucleotide
  # x <- vec[1]
  # same_as_x <- vec == x
  #
  # #test if all other nucleotides are the same (i.e. test if the character is non-segregating)
  # if (all(same_as_x)) return(FALSE)
  #
  # if sum(same_as_x) < 2

}

#' Compute the probability of at least one mutation ocurring
#'
#' The probability that no mutations occur in a coalescent with n sequences and mutation at constant rate \theta/2 is given by
#' p_no_mut = (1 / (1+theta)) * (2 / (2+theta)) * ... * (n-1 / (n-1+theta))
#' This function returns 1 - p_no_mut.
#'
#' @param n the number of sequences. Should be >= 2
#' @param theta the rate of mutation. Should be a non-negative float.
#' @result a float in the range (0,1)
probMutationOccurs <- function(n, theta){
  if(n < 2) stop('n should be >= 2')
  if(theta < 0) stop('theta should be non-negative.')
  return(1 - prod((function(i,x) i/(x+i))(1:(n-1),theta)))
}


#' Compute distance to nearest incompatibility
#'
#' Given a matrix where each row is a sequence/haplotype; and each collum coresponds to an alligned position: we for each position i=1...ncols determine the value of Z(i) = max{ 0 <= k <= ncols-i | coumn i compatible with all collumns in the index set i...i+k} (if the argument \code{dir} is \code{'left'} which is the default). If the argument \code{dir} is \code{'right'}, we instead determine Z'(i) max{ 0 <= k < i | coumn i compatible with all collumns in the index set i-k...i}. The general idea is that these values would be expected to behave like (cencored) samples from a geometric distribution with a mean inversely proportional to the recombination rate.
#'
#' @param SeqMatrix a matrix with rows as haplotypes and columns as characters ad individual loci.
#' @param dir the direction to scan in when finding incompatabilities. should be either 'left' or 'right'.
#' @param cencor_with_NA a boolean indicating if we should handle cencored
#' observations. Doing so will essentially lead to non-informative sites and
#' informative sites which are compatible with all subsequent sited being treated
#' equivalently (i.e. ignored). If set to false, *all* informative sites will be
#' assigned an integer value (which in turn means that analysis should take this
#' into account).
#'
#' @return a vector of (Z_1, Z_2, ... , Z_ncols) values computed as outlined above.
#' @export
compute_Z <- function(SeqMatrix, dir = 'left', cencor_with_NA = TRUE){
  L <- ncol(SeqMatrix)

  Z_vals <- vector(mode = 'numeric', length = L)

  if (dir == 'left'){
    z_step <- +1
  } else if (dir == 'right'){
    z_step <- -1
  } else stop('dir should be either \'left\' or \'right\' ' )

  i <- 1
  while(i <= L){
    if (is_informative_char(SeqMatrix[, i])){
      incompatibility_found <- FALSE
      Z_signed = z_step # initialize by taking a single step
      while(!incompatibility_found && 1 <= i + Z_signed  && i + Z_signed <= L){
        if (gamete_test(SeqMatrix[, i], SeqMatrix[, i + Z_signed], ancestral_types = c(FALSE, FALSE))){
          Z_signed <- Z_signed + z_step
        } else {
          #if(cencor_with_NA) Z_signed <- NA
          incompatibility_found <- TRUE
        }
      }

      #undo the effect of the step which lead us to an incompatiblilty
      Z_signed <- Z_signed - z_step

      # Change sign (nessecary when scanning right)
      Z_val <- Z_signed * z_step

      #If we scanned all the way to the boundary without finding an
      #incompatability, we may want to cencor.
      if(cencor_with_NA && !incompatibility_found){
        Z_val <- NA
      }

      Z_vals[i] <- Z_val
    } else {
      Z_vals[i] <- NA
    }
    i <- i+1
  }
  return(Z_vals)
}

#' Forward-backward recombination_rate determination
#'
#' Compute marginal distributions [i.e. compute P(R_t, Z_1:T) for t=1:T] of the
#' R-chain in the following HMM:
#'
#'     R_0 -> R_1 -> R_2 -> ... -> R_L
#'             |      |             |
#'             V      V             V
#'            Z_1    Z_2    ...    Z_L
#'
#' Where R is a MC with two states rho_low and rho_hot and transition matrix
#' P = 1-a  a
#'      b  1-b
#' (typically rho_hot > rho_low; and a<<b<<1 hold)
#' Z can either be NA (indicating that the chracter at a position is
#' non-informative (i.e. carries no mutation on an internal tree-branch)), or an
#' integer, indicating how many steps we have to take left/right without hitting
#' an incompatible character.
#'
#' @param SeqMatrix the matrix of observed sequences including non-segregating
#'   sites. Should be encoded as either a 0-1 matrix (1 entries indicating
#'   derived type and 0 ancestal) or as a boolean TRUE-FALSE matrix.
#' @param r_normal the recombination rate between two ajacent sites outside
#'   hotspots (on the timescale of the coalescent, i.e. a rate of 1 indicates
#'   that recombination-events occur at the same rate as a pairwise
#'   coalescence-event).
#' @param r_hot the recombination rate between two ajacent sites inside hotspots
#'   (on the timescale of the coalescent, i.e. a rate of 1 indicates that
#'   recombination-events occur at the same rate as a pairwise
#'   coalescence-event).
#' @param rate_norm2hot the rate at which the hidden state changes from normal
#'   recombination to hot. This should be inversely proportional to the expected
#'   distance between right-boundaries of hotspots.
#' @param rate_hot2norm the rate at which the hidden state changes from hot to
#'   normal recombination. This should be inversely proportional to the expected
#'   length of a hotspot.
#' @param theta (double) the presumed site-wise mutation rate. Defaults to wattersons estimate
#'   (n_seg_sites / ( sites * (n-1)st harmonic number)) if not supplied.
#'   not supplied.
#' @param r_factor Our modell presumes that the time untill an incompatability
#'   is observed will be geometrically distributed with a mean dependant on recombination and mutation rate. Since not all recombinations change the topoplogy, and not all mutations (once the topology has changed) cause incompatabilities, we scale the mean of this distribution by an arbitrary factor to account for this.
#'   It is an estimate of rate(topology_changes between two
#'   sites)/rate(recombination events in the ARG). Giving some thought to how to
#'   pick this factor may yield improvements.
#' @return A matrix with two rows and a column per site. Each column corresponds
#'   to the marginal distribution of the hidden state at that position (under
#'   the model outlined above and conditioned on the observed haplotypes).
#' @export
compute_marginal_recomb_state <- function(SeqMatrix,
                                          r_normal = 4e-4,
                                          r_hot = 4e-3,
                                          rate_norm2hot = 1e-5,
                                          rate_hot2norm = 1e-3,
                                          theta = NA,
                                          r_factor = 1.0){

  #trnasition matrix: state 1 is normal; state 2 is hot
  P <- matrix(data = c(1 - rate_norm2hot, rate_hot2norm, rate_norm2hot, 1- rate_norm2hot), nrow = 2, ncol = 2)

  #we chose the invariant distribution of P as the initial distribution
  alpha0 <- c(rate_hot2norm,rate_norm2hot) / (rate_norm2hot + rate_hot2norm) #

  #TODO: Do we have a better guess for r_factor, and will changing it improve our methodology?
  #r_factor <- 1.0

  #We scan both forward and backward
  emissions <- matrix(
    data = c(compute_Z(SeqMatrix, dir = 'left'),
             compute_Z(SeqMatrix, dir = 'right')),
    ncol = 2
  )
  #Z_forward  <- compute_Z(SeqMatrix, dir = 'left' )
  #Z_backward <- compute_Z(SeqMatrix, dir = 'right')

  #compute parameters for distriburtion of Z-values
  n <- nrow(SeqMatrix)
  sites <- ncol(SeqMatrix)

  #probability that a recombination occurs between two sites
  p_normal <- probMutationOccurs(n = n, theta = r_normal)
  p_hot <- probMutationOccurs(n = n, theta = r_hot)

  #probability that a mutation occurs at a new site
  if(is.na(theta)){# compute Watterons estimate of theta if nessecary.
    s = sum( apply(SeqMatrix, MARGIN = 2, FUN = any) )
    H <- sum(1/1:(n-1))
    theta <- (s / sites) / H
  }
  p_mutaion <- probMutationOccurs(n = n, theta = theta)

  # compute parameter for Z ~ geo(p_Z)
  p_Z_normal <- r_factor * p_normal * p_mutaion / (p_normal + p_mutaion)
  p_Z_hot <- r_factor * p_hot * p_mutaion / (p_hot + p_mutaion)

  density_Z_given_R <- function(y,x,t){
    result <- 1.0
    for(z in y){
      if (is.na(z)) next # if a value is missing, we do nothing
      p_Z <- ifelse(x==1, p_Z_normal, p_Z_hot)
      result <- result*dgeom(z,p_Z)
    }
    return(result)
  }

  log_likelihoods <- forwad_backward_logspace(emissions = emissions,
                                              alpha0 = alpha0,
                                              transition_probs = P,
                                              emission_density = density_Z_given_R)

  marginal_probs <- apply(X=log_likelihoods, MARGIN = 2, FUN = log_probs_2_normalised_probs)

  return(marginal_probs)
}

#' Compute marginal distribution, given marginal log_likelihoods
#'
#' Given log(p1), log(p2), ... , log(pN), return
#' p1 / (p1 + ... + pN), p2 / (p1 + ... + pN), ..., pN / (p1 + ... + pN)
#'
#' @param log_vector a vector of log-values x_1, ..., x_N
#' @return a vector where the kth entry is exp(x_k) / (exp(x_1) + ... + exp(x_N))
log_probs_2_normalised_probs <- function(log_vector){
  return( exp( log_vector - log_sum_exp(log_vector) ) )
}

#' Forward backward algoritm in logspace
#'
#' Given
#'  - a chain X with states indexed by 1...N and transition matrix P,
#'  - emissions Y_1:T with conditional marginal distributios
#'    P( Y_t = j | X_t = i ) = mu_i(j)
#'  - Ovserved emissions of the form Y_{1:T} = y_{1:T} (partial observations (i.e. not for all time-points) are allowed)
#' Compute log( P(X_t = i, Y_{1:T} = y_{1:T}) ) for t = 1...T using the
#' forward-backward algorithm.
#'
#' @param alpha0 the prior distribution of X
#' @param emissions The sequence of observations. Either a vector of length T,
#'   or an array with a first axis of length T; where T denotes the number of
#'   observations.
#' @param transition_probs An N-by-N stichastic matrix encoding the forward
#'   transition probabilities of the hidden chain X.
#' @param emission_density a function f(y,x,t) returning the density P(Y_t = y |
#'   X_t = x) (when Y is discrete) or lim_h->0 P(Y_t = y±h | X_t = x) / (2h)
#'   (when Y is continuous). Missing values should be handled by assigning 'Na'
#'   observations the weight 1 in discrete cases (in non-iscrete cases missing
#'   values nust be handled differently).
#' @return a vector with entries log( P(X_t = i, Y_{1:T} = y_{1:T}) ) for t =
#'   1...T
#' @export
# @param emission_probs A N-by-M matrix (where M is the number of states of Y)
# such that the (i,j)th entry corresponds to P(Y_t = j | X_t = i).
forwad_backward_logspace <- function(emissions, alpha0, transition_probs, emission_density){
  lalpha0 <- log(alpha0)
  lP <- log(transition_probs)
  #lmu <- log(emission_probs)
  N_states <- dim(lP)[1]

  #if emissions are a vector, the number of emissions is the length,
  if (is.vector(emissions)){
    N_emissions <- length(emissions)
    old_emissions_dim <- dim(emissions)
    dim(emissions) <- c(N_emissions,1)
  } else if (is.array(emissions)){
    N_emissions <- dim(emissions)[1]
  } else stop('emissions misty be passed as either a vector or an array')

  smoothed_values <- matrix(c(0), nrow = N_states, ncol = N_emissions)

  #compute forward predictions
  #lalphas <- matrix(c(0), nrow = N_states, ncol = N_emissions)

  forward_step <- function(lalpha, emission, t = 0){
    lalpha_prev <- lalpha
    j <- 1
    while(j <= N_states){

      #compute log(p_1j * a_1), ..., log(p_Mj * a_M)
      log_summands <-lP[,j]+lalpha_prev

      # compute sum_i(p_ij * a_i ) = P(X_t = j | X_(t-1) = i, Y_{1:t-1}=y_{1:t-1})
      lprob_sum <- log_sum_exp(log_summands)

      # compute log(a) = log(p(Y = y | X=j) * log(prob_sum)).
      # an emission of NA indicates that no observation occurred
      #if (!is.na(emission)) lalpha[j] <- lmu[j,emission] + lprob_sum
      # if (!is.na(emission)) lalpha[j] <- log(emission_density(y = emission, x = j, t = t)) + lprob_sum
      # else lalpha[j] <- lprob_sum
      lalpha[j] <- log(emission_density(y = emission, x = j, t = t)) + lprob_sum

      j <- j+1
    }
    return(lalpha)
  }

  backward_step <- function(emission,lbeta_now,t = 0){
    lbeta_next <- lbeta_now
    i <- 1
    while (i <= N_states){
      log_summands <- lP[i, ] + lbeta_next

      # if (!is.na(emission)){
      #   #log_summands <- log_summands + lmu[ , emission]
      #   log_summands <- log_summands +
      #     emission_density(y = emission, x = 1:N_states, t = t)
      # }
      log_summands <- log_summands + emission_density(y = emission, x = 1:N_states, t = t)
      lbeta_now[i] <- log_sum_exp(log_summands)
      i <- i+1
    }
    return(lbeta_now)
  }

  #Forward steps
  #lalphas <- Reduce(f = prediction_step, x = emissions, accumulate = TRUE, init = lalpha0)

  #Backward steps
  #lbetas  <- Reduce(f = backward_step,   x = emissions, accumulate = TRUE, init = rep(0, N_states), right = TRUE)

  #Perform forward and backward steps simultaneously
  lbetas <- matrix(c(0), nrow = N_states, ncol = N_emissions)
  lalphas <- matrix(c(0), nrow = N_states, ncol = N_emissions)

  lalphas[ , 1] <- forward_step(lalpha0,emissions[1, ])
  lbetas[ , N_emissions] <- backward_step(emissions[1, ] , rep(0, N_states))
  i <- 1
  while (i < N_emissions){
    lalphas[, 1 + i] <- forward_step(lalphas[, i], emissions[1 + i, ], t = i+1)
    i <- i+1
  }
  i <- N_emissions - 1
  while (0 < i){
    lbetas[ , i] <- backward_step(emissions[i, ], lbetas[ , i + 1], t = i)
    i <- i-1
  }

  #Smoothing
  smoothed_values <- lalphas + lbetas

  return(smoothed_values)
}
