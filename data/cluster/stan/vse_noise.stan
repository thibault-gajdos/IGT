// from https://github.com/CCS-Lab/hBayesDM/blob/develop/commons/stan_files/igt_pvl_delta.stan
// modified to include learning noise

/* functions { */
/*   real normal_bounded_rng(real lb, real ub) { */
/*     real p_lb = normal_cdf(lb, 0, 1); */
/*     real p_ub = normal_cdf(ub, 0, 1); */
/*     real x = uniform_rng(lb,ub); */
/*     real y = Phi(x); */
/*     return(y); */
/*   } */
/* } */

data {
  int<lower=1> T;
  int choice[T];
  real gain[T];
  real loss[T];
}
transformed data {
  vector[4] init_explore;
  vector[4] init_exploit;
  real noise[T];
  
  init_exploit  = rep_vector(0.0, 4);
  init_explore  = rep_vector(1.0, 4); 

  //generate noise 
  for (t in 1:T){
    //noise[i,t] =  normal_bounded_rng(-1,1);
    noise[t] =  normal_rng(0,1);
  }
}


parameters {
// Declare all parameters as vectors for vectorizing
 
 // Subject-level raw parameters (for Matt trick)
  real alpha_pr; // risk aversion
  real cons_pr; //consitency  
  real gamma_pr; // exploration rate
  real delta_pr; // decay rate
  real phi; // exploration bonus
  real zeta_pr; // noise 
  
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=0, upper=2> alpha;
  real<lower=0, upper=5> cons;
  real<lower=0, upper=1> gamma;
  real<lower=0, upper=1> delta;
  real<lower=0, upper=2> zeta;


  alpha  = Phi_approx(alpha_pr)*2;
  cons   = Phi_approx(cons_pr) * 5;
  gamma  = Phi_approx(gamma_pr);
  delta  = Phi_approx(delta_pr);
  zeta  = Phi_approx(zeta_pr) * 2;
}


model {
  // Define values
  vector[4] exploit;
  vector[4] explore;
  real curUtil;     // utility of curFb
  real absUtil;     // abs utility of curFb
  real theta;       // theta = 3^c - 1
  
  // individual parameters
  alpha_pr  ~ normal(0, 1);
  cons_pr   ~ normal(0, 1);
  gamma_pr ~ normal(0, 1);
  delta_pr ~ normal(0, 1);
  zeta_pr ~ normal(0, 1);
  phi ~ normal(0,1);

  // Initialize values
  theta = pow(3, cons) -1;
  exploit = init_exploit; // initial exploit values
  explore = init_explore; // initial explore values
    
  for (t in 1:T) {
      // softmax choice
      choice[t] ~ categorical_logit(theta * (exploit+explore));

      // curent utility
      curUtil = pow(gain[t], alpha) - pow(loss[t], alpha);
    
      // update exploit
      exploit *= delta;
      exploit[choice[t]] += curUtil + zeta * fabs(curUtil) * noise[t];

      // update explore
      for (k in 1:4){
	if (k == choice[t]) {
	  explore[k] =  0;  
	} else {
	  explore[k] += gamma * (phi - explore[k]);
	}	  
      }     
  }
}



generated quantities {  
  // For log likelihood calculation
  real log_lik;
  
  // For posterior predictive check
  real y_pred[T];

  // Set all posterior predictions to -1 (avoids NULL values)
  for (t in 1:T) {
    y_pred[t] = -1;
  }

  
  { // local section, this saves time and space
    // Define values
    vector[4] exploit;
    vector[4] explore;
    real curUtil;     // utility of curFb
    real absUtil;     // abs utility of curFb
    real theta;       // theta = 3^c - 1
    
    // Initialize values
    log_lik = 0;
    theta      = pow(3, cons) -1;
    exploit = init_exploit; // initial exploit values
    explore = init_explore; // inirial explore values
    
    for (t in 1:T) {
      // softmax choice
      log_lik += categorical_logit_lpmf(choice[t] | theta * (exploit + explore));
	  
      // generate posterior prediction for current trial
      y_pred[t] = categorical_rng(softmax(theta * (exploit + explore)));
      
      // curent utility
       curUtil = pow(gain[t], alpha) - pow(loss[t], alpha);
      
      // update exploit
      exploit *= delta;
      exploit[choice[t]] += curUtil + fabs(curUtil) * noise[t] * zeta;
	  
      // update explore
      for (k in 1:4) {
	if (k == choice[t]) {
	  explore[k] =  0;  
	} else {
	  explore[k] += gamma * (phi - explore[k]);
	    }	  
      }
    }
  }
}



