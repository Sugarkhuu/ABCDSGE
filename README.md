# ABCDGSE

Release 1 has the files for replicating the results for the DSGE model in the paper "Bayesian Indirect
Inference and the ABC of GMM" by Creel, Gao, Hong and Kristensen, http://arxiv.org/abs/1512.07385 The file DSGE.pdf is an extract, which
gives a description and the results.


The current version is updated to use parameter dependent moments and neural net statistics. This will be described in
the next version of the paper, which should be ready in June.


The estimator is an Approximate Bayesian Computing estimator,
computed using importance sampling and local linear nonparametric
regression. 

Requires Open MPI, Julia (with MPI and Distances packages), Octave (with MPI
package) and other supporting code available at https://github.com/mcreel/Econometrics

To replicate the results, execute "sh MasterScript" from the bash prompt.

For questions, write to michael.creel@uab.es
