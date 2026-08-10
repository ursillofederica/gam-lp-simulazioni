// lp_superficie_lkj.stan: superficie f(s,h) + Sigma LIBERA (LKJ).
// parametrizzazione centrata.

data {
  int<lower=1> TT; int<lower=1> K; int<lower=1> J; int<lower=1> H1;
  int<lower=1> N; int<lower=1> K_ctrl;
  matrix[N, J*K] TP; matrix[TT, H1] y_mat; matrix[TT, K_ctrl] X_ctrl;
  matrix[H1, J*K] TP_pos; matrix[H1, J*K] TP_neg; matrix[H1, J*K] TP_zero;
  real<lower=0> sigma_theta;
}

parameters {
  vector[J*K] gamma; real<lower=0,upper=1> delta;
  real<lower=0> sigma_h; real<lower=0> sd_shock;
  cholesky_factor_corr[H1] L_Omega; vector<lower=0>[H1] tau;
  vector[K_ctrl] theta;
}

transformed parameters {
  matrix[TT, H1] mu_mat;
  { vector[N] mv = TP*gamma; vector[TT] cp = X_ctrl*theta;
    for (t in 1:TT) for (h in 1:H1) mu_mat[t,h] = mv[(h-1)*TT+t] + cp[t]; }
  matrix[H1,H1] L_Sigma = diag_pre_multiply(tau, L_Omega);
}

model {
  for (t in 1:TT) y_mat[t]' ~ multi_normal_cholesky(mu_mat[t]', L_Sigma);
  delta ~ beta(2,2); sigma_h ~ gamma(2,0.1);
  L_Omega ~ lkj_corr_cholesky(2); tau ~ student_t(3,0,1);
  sd_shock ~ lognormal(0,0.7); theta ~ normal(0,sigma_theta);
  for (j in 1:J) {
    gamma[(j-1)*K+1] ~ normal(0, sigma_h/sqrt(1-square(delta)));
    for (k in 2:K) gamma[(j-1)*K+k] ~ normal(delta*gamma[(j-1)*K+(k-1)], sigma_h);
  }

  for (k in 1:K) for (j in 3:J)
    target += normal_lpdf(gamma[(j-1)*K+k]-2*gamma[(j-2)*K+k]+gamma[(j-3)*K+k] | 0, sd_shock);
}

generated quantities {
  vector[H1] irf_pos = TP_pos*gamma - TP_zero*gamma;
  vector[H1] irf_neg = TP_neg*gamma - TP_zero*gamma;
  vector[TT] log_lik;
  for (t in 1:TT) log_lik[t] = multi_normal_cholesky_lpdf(y_mat[t]' | mu_mat[t]', L_Sigma);
}
