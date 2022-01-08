// from https://github.com/CCS-Lab/hBayesDM/blob/develop/commons/stan_files/igt_pvl_delta.stan
// modified to include learning noise

data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  int choice[N, T];
  real outcome[N, T];
}
transformed data {
  vector[4] initV;
  real noise[N, T];
  
  initV  = rep_vector(0.0, 4); //initialize expected value to 0

  // generate noise -- one for each chain
  for (i in 1:N){
    for (t in 1:T){
      noise[i,t] =  normal_rng(0, 1);
	}
  }
}

parameters {
// Declare all parameters as vectors for vectorizing
  // Hyper(group)-parameters
  vector[5] mu_pr;
  vector<lower=0>[5] sigma;

  // Subject-level raw parameters (for Matt trick)
  vector[N] A_pr;
  vector[N] alpha_pr;
  vector[N] cons_pr;
  vector[N] lambda_pr;
  vector[N] zeta_pr;
}

transformed parameters {
  // Transform subject-level raw parameters
  vector<lower=0, upper=1>[N]  A;
  vector<lower=0, upper=2>[N]  alpha;
  vector<lower=0, upper=5>[N]  cons;
  vector<lower=0, upper=10>[N] lambda;
  vector<lower=0, upper=1.5>[N] zeta;

  for (i in 1:N) {
    A[i]      = Phi_approx(mu_pr[1] + sigma[1] * A_pr[i]);
    alpha[i]  = Phi_approx(mu_pr[2] + sigma[2] * alpha_pr[i]) * 2;
    cons[i]   = Phi_approx(mu_pr[3] + sigma[3] * cons_pr[i]) * 5;
    lambda[i] = Phi_approx(mu_pr[4] + sigma[4] * lambda_pr[i]) * 10;
    zeta[i]  = Phi_approx(mu_pr[5] + sigma[5] * zeta_pr[i])*1.5;
  }
}

model {
// Hyperparameters
  mu_pr  ~ normal(0, 1);
  sigma ~ normal(0, 0.2);

  // individual parameters
  A_pr      ~ normal(0, 1);
  alpha_pr  ~ normal(0, 1);
  cons_pr   ~ normal(0, 1);
  lambda_pr ~ normal(0, 1);
  zeta_pr ~ normal(0, 1);

  for (i in 1:N) {
    // Define values
    vector[4] ev;
    real curUtil;     // utility of curFb
    real theta;       // theta = 3^c - 1
    real epsilon;  //abs value pred error

    // Initialize values
    theta = pow(3, cons[i]) -1;
    ev = initV; // initial ev values
    
    for (t in 1:Tsubj[i]) {
      // softmax choice
      choice[i, t] ~ categorical_logit(theta * ev);

      if (outcome[i, t] >= 0) {  // x(t) >= 0
        curUtil = pow(outcome[i, t], alpha[i]);
      } else {                  // x(t) < 0
        curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
      }

      // compute abs value pred error
      if (curUtil - ev[choice[i, t]] >= 0){
	epsilon = curUtil - ev[choice[i, t]];
      }else{
	epsilon =  ev[choice[i, t]] - curUtil;
      }
    
      // update EV
      ev[choice[i, t]] += A[i] * (curUtil - ev[choice[i, t]]) + epsilon*zeta[i]*noise[i,t];
    }
  }
}

generated quantities {
  // For group level parameters
  real<lower=0, upper=1>  mu_A;
  real<lower=0, upper=2>  mu_alpha;
  real<lower=0, upper=5>  mu_cons;
  real<lower=0, upper=10> mu_lambda;
  real<lower=0, upper=1> mu_zeta;
  

  // For log likelihood calculation
  real log_lik[N];

  // For posterior predictive check
  real y_pred[N, T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }

  mu_A      = Phi_approx(mu_pr[1]);
  mu_alpha  = Phi_approx(mu_pr[2]) * 2;
  mu_cons   = Phi_approx(mu_pr[3]) * 5;
  mu_lambda = Phi_approx(mu_pr[4]) * 10;
  mu_zeta = Phi_approx(mu_pr[5]);

  { // local section, this saves time and space
    for (i in 1:N) {
      // Define values
      vector[4] ev;
      real curUtil;     // utility of curFb
      real theta;       // theta = 3^c - 1
      real epsilon;  //abs value pred error

      // Initialize values
      log_lik[i] = 0;
      theta      = pow(3, cons[i]) -1;
      ev         = initV; // initial ev values

      for (t in 1:Tsubj[i]) {
        // softmax choice
        log_lik[i] += categorical_logit_lpmf(choice[i, t] | theta * ev);

        // generate posterior prediction for current trial
        y_pred[i, t] = categorical_rng(softmax(theta * ev));

        if (outcome[i, t] >= 0) {  // x(t) >= 0
          curUtil = pow(outcome[i, t], alpha[i]);
        } else {                  // x(t) < 0
          curUtil = -1 * lambda[i] * pow(-1 * outcome[i, t], alpha[i]);
        }

        // absolute value pred error
	if (curUtil - ev[choice[i, t]] >= 0){
	epsilon = curUtil - ev[choice[i, t]];
	}else{
	  epsilon =  ev[choice[i, t]] - curUtil;
      }

	// choice
	ev[choice[i, t]] += A[i] * (curUtil - ev[choice[i, t]]) + epsilon*zeta[i]*noise[i,t];
      }
    }
  }
}

