# GLSP-count: Global-local shrinkage priors for analyzing sequence of counts
This package implements a new shrinkage method based on global-local shrinkage priors for Poisson likelihood via Markov Chain Monte Carlo algorithm, as proposed by the following paper.

Hamura, Y., Irie, K. and Sugasawa, S. (2019). On Global-local Shrinkage Priors for Count Data. to appear in *Bayesian Analysis*.(https://arxiv.org/abs/1907.01333)

The repository includes the following files.

* GLSP-function.R: Script implementing the proposed global-local shrinkage prior method 
* data-generation.R: Script generating the simulated dataset
* Example.RData: Simulated dataset 
* demo.R: Script for an application of the proposed method to the simulated dataset
