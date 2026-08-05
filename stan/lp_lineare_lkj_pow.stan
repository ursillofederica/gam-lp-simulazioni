// ============================================================================
// lp_lineare_lkj.stan — Cap. 5, modello lineare + Sigma LIBERA (LKJ)
//
// DERIVATO da cap5/stan/lp_lineare_strutturata.stan: UNICA differenza il
// blocco di Sigma, sostituito con la parametrizzazione libera del cap. 4
// (eq. sigma_lkj): Sigma = diag(tau) * R * diag(tau),
// R ~ LKJ(eta = 2), tau_h ~ t3(0,1)+.
// In Stan: fattore di Cholesky della correlazione (L_Omega) e scale per
// orizzonte (tau); L_Sigma = diag_pre_multiply(tau, L_Omega).
// Tutto il resto (NC su gamma, media, controlli per orizzonte, prior,
// log_lik per la ripesatura) e' identico al file strutturato.
// ============================================================================

data {
  int<lower=1> TT;             // n. righe del sistema (T)
  int<lower=1> K;              // dimensione base su h
  int<lower=1> H1;             // n. orizzonti (H+1)
  int<lower=1> N;              // = TT * H1
  int<lower=1> K_ctrl;         // n. controlli (P1: 5 ritardi; P2: 1)
  matrix[N, K] TP;             // design: shock_t * psi_k(h)
  matrix[TT, H1] y_mat;
  matrix[TT, K_ctrl] X_ctrl;
  matrix[H1, K] TP_s1;         // base a shock unitario: irf = TP_s1 * gamma
  real<lower=0> sigma_theta;
  real<lower=0,upper=1> c_pow;   // potenza della pseudo-verosimiglianza
}

parameters {
  vector[K] gamma_raw;             // NC: standard normal iid
  real<lower=0,upper=1> delta;     // decadimento AR(1) lungo k
  real log_sigma_h;                // log-reparam scala di gamma
  cholesky_factor_corr[H1] L_Omega; // correlazione tra orizzonti (libera)
  vector<lower=0>[H1] tau;         // scale per orizzonte
  matrix[K_ctrl, H1] theta;        // controlli LIBERI per orizzonte
}

transformed parameters {
  real<lower=0> sigma_h = exp(log_sigma_h);

  // NC reconstruction di gamma (AR(1) lungo k)
  vector[K] gamma;
  {
    real scale_first = sigma_h / sqrt(1 - square(delta));
    gamma[1] = scale_first * gamma_raw[1];
    for (k in 2:K) {
      gamma[k] = delta * gamma[k-1] + sigma_h * gamma_raw[k];
    }
  }

  // Media: parte shock (TP * gamma) + controlli per orizzonte
  matrix[TT, H1] mu_mat;
  {
    vector[N] mu_vec = TP * gamma;
    matrix[TT, H1] ctrl_part = X_ctrl * theta;   // [TT, H1]
    for (t in 1:TT) {
      for (h in 1:H1) {
        mu_mat[t, h] = mu_vec[(h - 1) * TT + t] + ctrl_part[t, h];
      }
    }
  }

  // Sigma libera: L_Sigma = diag(tau) * L_Omega
  matrix[H1, H1] L_Sigma = diag_pre_multiply(tau, L_Omega);
}

model {
  // Verosimiglianza (pseudo-): righe indipendenti per costruzione del criterio
  for (t in 1:TT) {
    target += c_pow * multi_normal_cholesky_lpdf(y_mat[t]' | mu_mat[t]', L_Sigma);
  }

  // Prior — cap. 4.4/4.5
  gamma_raw ~ std_normal();
  delta ~ beta(2, 2);
  log_sigma_h ~ normal(log(0.05), 1);
  L_Omega ~ lkj_corr_cholesky(2);
  tau ~ student_t(3, 0, 1);
  to_vector(theta) ~ normal(0, sigma_theta);
}

generated quantities {
  // IRF al contrasto unitario
  vector[H1] irf = TP_s1 * gamma;

  // log g per riga: sum(log_lik) = log g(y; theta) per la ripesatura (§ 5.1)
  vector[TT] log_lik;
  for (t in 1:TT) {
    log_lik[t] = multi_normal_cholesky_lpdf(y_mat[t]' | mu_mat[t]', L_Sigma);
  }

  // Posterior predictive
  matrix[TT, H1] y_rep_mat;
  for (t in 1:TT) {
    y_rep_mat[t] = (multi_normal_cholesky_rng(mu_mat[t]', L_Sigma))';
  }
}
