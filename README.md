# HotBird-HMM

A project aiming to develop novel methods for finding recombination hotspots
with a novel HMM-based approach.

Joint work w. Berklund, Terkelsen and Hobolth  @BIRC (Aarhus Uni)

## Overview of code:

 * **R/HMM.R** Code for sampling from and computing likelihoods for our model of
  pairwise sequence heterozygosity. Depends on code written by Asger Hobolth
  which has not been uploaded to this repositiory (as It's not mine). Ask asger,
  or implement something yourself if you want this to work.

## Project Roadmap

 1. Implement basics of HMM:
    1. Sampler [done]
    2. Likelihood computation (via forward algorithm) [done]

 2. Work on maximum likelihood estimator
    1. Implement ML estimation of recombination rate [done, but needs improvement]
    2. Test performance and power of method
    3. Expore improvements to ML algorithm, if nessecary (e.g. use a Baum–Welch/EM- type algorithm instead of newtonian method currently employed)

 3. Search for hotspots
    1. Implement scanning algorithm
    2. Test on simulated data
    3. Apply to real data