functions {
  // Work on the log concentration scale, but handle zero explicitly.
  real response_mean(real x, real bottom, real top,
                     real log_ec50, real hill) {
    if (x == 0) return top;
    return bottom + (top - bottom)
                    * inv_logit(-hill * (log(x) - log_ec50));
  }
}
data {
  int<lower=1> N;
  int<lower=1> J;
  array[N] int<lower=1, upper=J> compound;
  vector<lower=0>[N] concentration;
  vector[N] y;
  int<lower=1> M;
  array[M] int<lower=1, upper=J> compound_grid;
  vector<lower=0>[M] concentration_grid;
}
parameters {
  vector[J] bottom;
  vector[J] log_gap;
  vector[J] log_ec50;
  vector[J] log_hill;
  vector<lower=0>[J] sigma;
}
transformed parameters {
  vector[J] top = bottom + exp(log_gap);
  vector<lower=0>[J] ec50 = exp(log_ec50);
  vector<lower=0>[J] hill = exp(log_hill);
  vector[N] mu;
  for (n in 1:N) {
    int j = compound[n];
    mu[n] = response_mean(concentration[n], bottom[j], top[j],
                          log_ec50[j], hill[j]);
  }
}
model {
  // Explicit priors; log_gap ensures top > bottom.
  // These priors are choices for this reanalysis, not the paper's priors.
  target += normal_lpdf(bottom | 0, 0.2);
  target += normal_lpdf(log_gap | 0, 0.25);
  target += normal_lpdf(log_ec50 | log(300), 1);
  target += normal_lpdf(log_hill | log(1.5), 0.8);
  target += normal_lpdf(sigma | 0, 0.15);
  
  for (n in 1:N)
    target += normal_lpdf(y[n] | mu[n], sigma[compound[n]]);
}

generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;
  vector[M] mu_grid;
  vector[M] y_grid_rep;
  for (n in 1:N) {
    y_rep[n] = normal_rng(mu[n], sigma[compound[n]]);
    log_lik[n] = normal_lpdf(y[n] | mu[n], sigma[compound[n]]);
  }
  for (m in 1:M) {
    int j = compound_grid[m];
    mu_grid[m] = response_mean(concentration_grid[m], bottom[j], top[j],
                               log_ec50[j], hill[j]);
    y_grid_rep[m] = normal_rng(mu_grid[m], sigma[j]);
  }
}
