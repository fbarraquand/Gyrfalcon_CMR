
data {
  int<lower=0> nind;
  int<lower=0> n_occasions;
  int<lower=1,upper=3> y[nind, n_occasions];
  int<lower=0> f[nind];
  int<lower=1,upper=2> xj[nind,n_occasions-1]; //stage matrix, 1 if juvenile 2 if adult
  real<lower=0> prey_abund[n_occasions]; 
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
}

parameters {
  real<lower=0,upper=1> mean_s2; // Mean adulty survival 
  real<lower=0,upper=1> mean_eta;  // Mean fidelity
  real<lower=0,upper=1> mean_r;    // Mean recovery
  real<lower=0,upper=1> mean_p;    // Mean recapture
  real<lower=0> gamma_prey;      // slope logistic -- not necessarily positive, but shouldn't we assume so biologically? 
  //For this parameterization (gamma,mu) the non-positive version does not converge well -> positive constraint. 
  real<lower=0> mu_prey;         // Mean for the logistic
}

transformed parameters {
  matrix<lower=0,upper=1>[n_occ_minus_1,2] s; // Survival probability
  vector<lower=0,upper=1>[n_occ_minus_1] eta; // Survival of carcass
  vector<lower=0,upper=1>[n_occ_minus_1] r; // Recovery probability
  vector<lower=0,upper=1>[n_occ_minus_1] p; // Recapture/resighting probability
  simplex[3] ps[3, nind, n_occ_minus_1];   // ps[j, i, t - 1, k] i: indiv, t: year, j: previous state, k: next state
  simplex[3] po[3, nind, n_occ_minus_1];


  // Constraints
  for (t in 1:n_occ_minus_1) {
    s[t,1] = inv_logit(gamma_prey*(prey_abund[t]-mu_prey));   //previously mean_s[1], variable here;
    s[t,2] = mean_s2;
    eta[t] = mean_eta;
    r[t] = mean_r;
    p[t] = mean_p;
  }


    // Define state-transition and observation matrices 	
    for (i in 1:nind){

    // Define probabilities of state S(t+1) given S(t)
    for (t in 1:n_occ_minus_1){ // starts at 1 instead of f[i] to define the whole matrices
    ps[1,i,t,1] = s[t,xj[i,t]]; 		// xj[i,t], the stage matrix, enters here. 
    ps[1,i,t,2] = 1.0-s[t,xj[i,t]];  //shouldn't it be t-1? -> no that's OK because we use ps[1,i,t-1,2] in the equation later on. 
    ps[1,i,t,3] = 0.0;
    ps[2,i,t,1] = 0.0;
    ps[2,i,t,2] = eta[t];
    ps[2,i,t,3] = 1.0-eta[t];
    ps[3,i,t,1] = 0.0;
    ps[3,i,t,2] = 0.0;
    ps[3,i,t,3] = 1.0;
    
    // Define probabilities of O(t+1) given S(t+1) //notice the timelag, cf below. 
    po[1,i,t,1] = p[t];
    po[1,i,t,2] = 0.0;
    po[1,i,t,3] = 1-p[t];
    po[2,i,t,1] = 0.0;
    po[2,i,t,2] = r[t];
    po[2,i,t,3] = 1.0-r[t];
    po[3,i,t,1] = 0.0;
    po[3,i,t,2] = 0.0;
    po[3,i,t,3] = 1.0;

    } //t
    } //i
  
}

model {
     real acc[3];
     vector[3] gamma[n_occasions];

     // Uniform priors are implicitly defined.
     mean_s2 ~ uniform(0, 1); // Prior for mean adult survival
     mean_eta ~ uniform(0,0.5);  // Prior for mean carcass survival prob
     mean_r ~ uniform(0, 0.5);   // Prior for mean recovery
     mean_p ~ uniform(0, 0.5);   //Prior for mean recapture
    // Logistic model on juvenile survival probability
     gamma_prey ~ exponential(100); //better for >0 constraint that previous normal(0,1); mean = 0.01
     mu_prey ~ normal(6,2); // more vague prior

     // Likelihood
     // Forward algorithm derived from Stan Modeling Language
     // User's Guide and Reference Manual

     for (i in 1:nind)
     {

        if ((f[i] > 0) && (f[i]<n_occasions)) //some individuals have data from the last year. 
        { //check the definition of f -> do we have the same? Yes. 
           for (k in 1:3)
                gamma[f[i], k] = (k == y[i, f[i]]); //Puts gamma to the observed state probability at first capture. Perhaps check that. 

      for (t in (f[i] + 1):n_occasions) {
        for (k in 1:3) {
          for (j in 1:3)
            acc[j] = gamma[t - 1, j] * ps[j, i, t-1, k]
                    * po[k, i, t-1, y[i, t]]; 
                    // beware of that mismatch in time -- means that p[t] is in fact one lag after t in time? So far as these are constant no pb though. 
          gamma[t, k] = sum(acc)+ 0.000001;// add small number in case we have a numerical zero problem. 
         }
      }
      target += log(sum(gamma[n_occasions]));
     }
  }

}


