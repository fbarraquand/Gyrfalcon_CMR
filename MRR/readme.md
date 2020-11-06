The paper presents three models: A, B and C

* Model A (without covariates) in `Stan_juvad`. 

* Model B (with the ptarmigan prey as covariate only) in `Stan_juvad_wcovariates`. Logistic link to juvenile survival, with 1 covariate. Several versions, depending on both the formulation of the logistic (and constraints on parameters) and whether the covariate is normalised. v0: factorized + not-normalised + gamma>0. v1:  v0: factorized + not-normalised + gamma unconstrained. 

* Model C (with the three variables, adding two weather covariates) in `Stan_juvad_wcovariates`. Logistic link to juvenile survival, with 3 covariates.

The other subrepositories contain:

* Code for Kéry & Schaub's Mark-Recapture-Recovery model in `Stan_example`, translated into Stan by Hiroki Ito and Bob Carpenter. 

* `Jags`: similar models in JAGS, yielding similar results. However, speed was slow (due to computation of all the latent states through the Gibbs sampler) and convergence unsatisfactory (the latent states did make sense though). Hence the move to Stan, with use of the Forward algorithm. 


