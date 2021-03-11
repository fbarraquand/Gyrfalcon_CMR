     
We now use more informative priors and use the beta distribution for parameters between zero and one
* mean_s[1] ~ beta(2,4) instead of uniform(0, 1); // Prior for mean juvenile survival 
* mean_s[2] ~ beta(4,2) instead of uniform(0, 1); // Prior for mean adult survival
* mean_eta ~ beta(1,5) instead of uniform(0,0.5); // Prior for mean carcass survival prob
* mean_r ~ beta(1,5) instead of uniform(0, 0.5); // Prior for mean recovery
* mean_p ~ beta(1,5) instead of uniform(0, 0.5); //Prior for mean recapture

For logistic models on juvenile survival probability
* mu_juvsurv ~ normal(0,1); // instead of 0.5
* beta[1] ~ normal(0,1); // We try not to change everything at the same time, so we keep the same weakly informative priors there
* beta[2] ~ normal(0,1);
* beta[3] ~ normal(0,1);
