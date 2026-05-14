data {
  int<lower=1> N;

  int<lower=0> treatment_events[N];
  int<lower=0> treatment_total[N];

  int<lower=0> control_events[N];
  int<lower=0> control_total[N];

  real logRR_prior_mu;
  real<lower=0> logRR_prior_sigma;

  real<lower=0> tau_prior_sigma;
}

parameters {
  real theta;
  real<lower=0> tau;

  vector[N] study_effect;
  real mu;
}

transformed parameters {

  vector[N] logRR;
  vector<lower=0, upper=1>[N] p_treat;
  vector<lower=0, upper=1>[N] p_ctrl;

  for (i in 1:N) {

    logRR[i] = theta + study_effect[i];

    p_treat[i] = inv_logit(mu + logRR[i] / 2);
    p_ctrl[i]  = inv_logit(mu - logRR[i] / 2);
  }
}

model {

  // treatment effect prior (IMPORTANT FIX)
  theta ~ normal(logRR_prior_mu, logRR_prior_sigma);

  mu ~ normal(0, 10);

  tau ~ normal(0, tau_prior_sigma);

  study_effect ~ normal(0, tau);

  for (i in 1:N) {
    treatment_events[i] ~ binomial(treatment_total[i], p_treat[i]);
    control_events[i] ~ binomial(control_total[i], p_ctrl[i]);
  }
}