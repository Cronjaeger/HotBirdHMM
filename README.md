# HotBirdHMM

A project aiming to develop novel methods for finding recombination hotspots
with a novel HMM-based approach.

Joint work w. Berklund, Terkelsen and Hobolth  @BIRC (Aarhus Uni)

## Overview of code:

 * **R/TransMat.R** (written by Asger Hobolth) contains aux. functions for
  computing transition-rates of the hidden states in our HMM. This is done via
  matrix exponentiation using the '[expm](https://cran.r-project.org/web/packages/expm/index.html)' library.

 * **R/HMM.R** Code for sampling from and computing likelihoods for our model of
  pairwise sequence heterozygosity.

## Project Roadmap

 1. Implement basics of HMM:
    1. Sampler [done]
    2. Likelihood computation (via forward algorithm) [done]

 2. Work on maximum likelihood estimator
    1. Implement ML estimation of recombination rate [done]
    2. Test performance and power of method [initial testing done]
    3. Explore improvements to ML algorithm, if necessary (e.g. use a Baum–Welch/EM- type algorithm instead of Newtonian method currently employed) [currently not under consideration]

 3. Search for hotspots
    1. Implement scanning algorithm [next step]
    2. Test on simulated data
    3. Apply to real data
