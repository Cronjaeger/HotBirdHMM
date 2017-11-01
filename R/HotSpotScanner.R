# Scan for hotspots
#
# (C) Mathias C. Cronjäger 12/9/2017
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




#--------------------------
# Step 1: Setting the stage
#--------------------------

#-------------------
#  1.1: Import code
#-------------------
#source('./R/HMM.R')

#---------------------
# Auxiliary functions
#---------------------

#' Get starting indices for contiguous runs
#'
#' @param vec a numerical vector
#' @return a vector of all indicees i such that vec[i] != vec[i-1]. The first index will always be 1.
#' @examples
#' > v <- c(1,1,1,1,0,0,1)
#' > get_run_indices(v)
#' [1] 1 5 7
#' @export
find_start_of_runs <- function(vec){
  return(c(1,(1 + 1:length(vec))[diff(vec) != 0]))
}

#' Compute equidistant (up to rounding) breakpoints
#'
#' takes the contiguous integers n...N and finds the breakpoints if they are to
#' be divided into a given number of cintiguous intervals of integers of equal
#' (up to rounding) length.
#'
#' @param upper the upper bound
#' @param lower the lower bound
#' @param breaks the numbr of intervals into which {upper...lower} is to be
#' divided.
#' @return a vector of length \code{breaks - 1} containing the first integer
#' after (or equal to) each breakpoint.
#' @export
#' @examples
#' > compute_equidistant_breaks(upper = 99, lower = 0,3)
#' [1] 34 67
#' > compute_equidistant_breaks(upper = 99, lower = 0,5)
#' [1] 20 40 60 80
#' > compute_equidistant_breaks(upper = 100, lower = 1,5)
#' [1] 21 41 61 81
compute_equidistant_breaks <- function(upper, lower = 1, breaks = 2){
  if (breaks == 1) return(c())
  return((find_start_of_runs(cut(lower:upper, breaks, labels = FALSE)) - 1 + lower)[-1])
}

#' Convert data in MS format to raw sequences
#'
#' In MS sequences are encoded only by specifying seg. sites (in a matrix where
#' rows are sequences and collumns are segregating sites) along with a vector
#' assigning each collumn a position p in the interval [0,1]. This function
#' converts to a matrix which includes non-segregating sites (the total numebr
#' of sites has to be manually specified, and is presumed to be known to the
#' user). Each segregating site is the ninserted at position round(p*sites).
#'
#' @param sites the number of total (segregating and non-segregating) sites. The
#'   output will have this number of columns.
#' @param position vector a vector of floats between 0 and 1 indicating the
#'   scaled location of each segregating site.
#' @param segSiteMatrix A 0-1 matrix where each row is an individual, and each
#'   column is a site. An entry of 1 in the (i,j)th position indicates that the
#'   ith individual is segregating at the jth segregating site.
#' @return A logical (to conserve memmory) matrix with \code{sites} columns and
#'   the same numebr of rows as \code{segSiteMatrix}.
#' @export
ms_out2SeqArray <- function(sites,position_vector,segSiteMatrix){
  n = nrow(segSiteMatrix)
  SeqArray <- matrix(data = logical(length = n * sites), nrow = n, ncol = sites)
  for (seg_index in 1:length(position_vector)){
    real_index <- round(position_vector[seg_index]*sites)
    SeqArray[,real_index] <- as.logical(segSiteMatrix[,seg_index])
  }
  return(SeqArray)
}

#----------------------
# Main functions
#----------------------

#' Sample hidden states with hotspots
#'
#' @param rho_normal recombination rate outside hotspots
#' @param rho_hot recombination rate in hotspots
#' @param hotspot_indicator A binary vector hotspot_indicator[i]==TRUE iff the position $i$ lies within a hotspot.
#' @param M the number of hidden states to be used
#' @param x0 the initial value; if not specified, x0 will be uniformly sampled from 1...M
#' @return A vector (X1, ... ,Xn) where n = length(hotspot_indicator) where Xi | X(i-1) is sampled according to recombination rate rho_high if hotspot_indicator[i] = TRUE and according to recombination rate rho_low if hotspot_indicator[i] == FALSE.
#' @export
sample_X_with_hotspots <- function(rho_normal,rho_hot,hotspot_indicator,M,x0 = NULL){

  if (is.null(x0)) x0 <- sample(1:M, 1)

  Px_normal <- compute_Px(M,rho_normal)
  Px_hot <- compute_Px(M,rho_hot)

  L <- length(hotspot_indicator)

  X <- vector(mode = 'integer', length = L)

  if (L == 0) return(X) # A special case.

  update <- function(x_prev,is_hot){
    if (is_hot){
      prob <- Px_hot[x_prev, ]
    } else {
      prob <- Px_normal[x_prev, ]
    }
    x_new <- sample(x=M, size = 1, prob = prob)
    return(x_new)
  }

  X[1] <- update(x0,hotspot_indicator[1])

  i <- 2
  while(i <= L){
    X[i] <- update(X[i-1],hotspot_indicator[i])
    i <- i+1
  }

  return(X)
}

#' Sample both hidden states and emissions
#'
#' @param theta the mutation rate
#' @param rho_normal recombination rate outside hotspots
#' @param rho_hot recombination rate in hotspots
#' @param hotspot_indicator A binary vector hotspot_indicator[i]==TRUE iff the position $i$ lies within a hotspot.
#' @param M the number of hidden states to be used
#' @param x0 the initial value; if not specified, x0 will be uniformly sampled from 1...M
#' @return A data-frame with two collumns named 'X' and 'Y'. The number of rows is equal to the length of \code{hotspot_indicator}. The first row is the first (X,Y)-pair, the seccond row is the second, and so on.
#' @export
sample_XandY_with_hotspots <- function(theta, rho_normal, rho_hot, hotspot_indicator, M, x0 = NULL){
  X <- sample_X_with_hotspots(rho_normal, rho_hot, hotspot_indicator, M, x0)
  Y <- sample_YgivenX(X, theta, M) # from HMM.R
  return(data.frame(X=X,Y=Y))
}

#' Compute recombination rate in contiguous windows
#'
#' @param Y_seq the sequence of contiguously observed 0-1 values indicating if
#' a site is heterozygous (y == 1) or homozygous (y == 0).
#' @param breakpoints a strictly increasing vector containing the first index in
#' each window, e.g. c(26, 51, 76) if Y_seq contains 100 observations that are
#' to be divided into four windows of equal size (i.e. 1:25, 26:50, 51:75, and
#' 76:100).
#' @param M the number of hidden states to be used in HMM governing rounded
#' tree-heights.
#' @param theta the mutation rate to be presumed. If none is supplied,
#' Waterson's estimator will be used.
#' @return a vector of length breakpoints + 1 containig the maximum likelihood estimate within each window. The first window always starts at index 1, and the last window always ends on index length(Y_seq).
#' @export
compute_windowed_ML_rho <- function(Y_seq, breakpoints, M, theta = NULL) {

  L <- length(Y_seq)

  if (any( breakpoints  < 1 || breakpoints  > L )){
    warning('Some breakpoints out of range 1...length(Y_seq); these will be ignored.')
  }

  breakpoints <- Filter(f = function(i) 1 < i && i < L, breakpoints)
  breakpoints <- c(1,breakpoints,L)
  if (!(all(diff(breakpoints) > 0))) stop('breakpoints must be a strictly increasing sequence')

  # if (breakpoints[1] != 1) breakpoints <- c(1,breakpoints)
  # if (breakpoints[length(breakpoints)] != length(Y)) breakpoints <- c(breakpoints,length[Y])

  if (is.null(theta)) theta <- theta_watterson(s = sum(Y_seq == 1), n = 2, L = length(Y_seq))

  rhos <- vector(mode = 'numeric', length = length(breakpoints) - 1)
  bounds_upper <- vector(mode = 'integer', length = length(breakpoints) - 1)
  bounds_lower  <-vector(mode = 'integer', length = length(breakpoints) - 1)
  i <- 1
  while(i < length(breakpoints)){
    left_endpoint <- breakpoints[i]
    right_endpoint <- breakpoints[i+1] - 1
    ml_rho <- compute_ML_rho(Y = Y_seq[left_endpoint:right_endpoint], M = M, theta = theta)
    rhos[i] <- ml_rho
    bounds_lower[i] <- left_endpoint
    bounds_upper[i] <- right_endpoint
    i <- i+1
  }

  result = data.frame(lower = bounds_lower, upper = bounds_upper, rho = rhos)

  return(result)

  # return(rhos)
}


locate_hotspots <- function(Y_sequences){
  stop('NOT IMPLEMENTED')
}

