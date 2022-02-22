// from https://github.com/CCS-Lab/hBayesDM/blob/develop/commons/stan_files/igt_pvl_delta.stan
// modified to include learning noise


data {
  int<lower=1> T;
  int choice[T];
  real outcome[T];
  real gain[T];
  real loss[T];
}
transformed data {
  vector[4] initV;
  real noise[T];
  
  initV  = rep_vector(0.0, 4); //initialize expected value to 0

  //generate noise 
  for (t in 1:T){
      noise[t] =  normal_rng(0,1);
  }
}

parameters {
  //Rraw parameters (for Matt trick)
  real A_pr;
  real alpha_pr;
  real cons_pr;
  real zeta_pr;
  real lambda_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  real<lower=0, upper=1>  A;
  real<lower=0, upper=1>  alpha;
  real<lower=0, upper=5>  cons;
  real<lower=0, upper=10>  lambda;
  real<lower=0, upper=2>  zeta;

  A      = Phi_approx(A_pr);
  alpha  = Phi_approx(alpha_pr);
  cons   = Phi_approx(cons_pr) * 5;
  zeta   = Phi_approx(zeta_pr) * 2;
  lambda   = Phi_approx(lambda_pr) * 10;
}

model {
  // Define values
  vector[4] ev;
  real curUtil;     // utility of curFb
  real theta;       // theta = 3^c - 1
  real epsilon;  //abs value pred error
  
// individual parameters
  A_pr      ~ normal(0, 1);
  alpha_pr  ~ normal(0, 1);
  cons_pr   ~ normal(0, 1);
  lambda_pr   ~ normal(0, 1);
  zeta_pr   ~ normal(0, 1);

  
  
  // Initialize values
  theta = pow(3, cons) -1;
  ev = initV; // initial ev values
  
  for (t in 1:T) {
    // softmax choice
    choice[t] ~ categorical_logit(theta * ev);

     curUtil = pow(gain[t], alpha) - lambda * pow(loss[t], alpha);

    //compute abs value pred error
    if (curUtil - ev[choice[t]] >= 0){
      epsilon = curUtil - ev[choice[t]];
    }else{
      epsilon =  ev[choice[t]] - curUtil;
    }
    
    // update EV
    ev[choice[t]] += A * (curUtil - ev[choice[t]]) +  noise[t] * epsilon * zeta; //learning
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
      vector[4] ev;
      real curUtil;     // utility of curFb
      real theta;       // theta = 3^c - 1
      real epsilon;  //abs value pred error

      // Initialize values
      log_lik = 0;
      theta      = pow(3, cons) -1;
      ev         = initV; // initial ev values

      for (t in 1:T) {
        // softmax choice
        log_lik += categorical_logit_lpmf(choice[t] | theta * ev);

        // generate posterior prediction for current trial
        y_pred[t] = categorical_rng(softmax(theta * ev));

        curUtil = pow(gain[t], alpha) - lambda * pow(loss[t], alpha);

        //absolute value pred error
	if (curUtil - ev[choice[t]] >= 0){
	  epsilon = curUtil - ev[choice[t]];
	}else{
	  epsilon =  ev[choice[t]] - curUtil;
      }

	// choice
	ev[choice[t]] += A * (curUtil - ev[choice[t]])+ noise[t] * epsilon * zeta; //learning
      }
  }
}


