// ============================================================================
// lp_lineare_strutturata_pow.stan — variante POWER POSTERIOR di
// lp_lineare_strutturata.stan: unica differenza la verosimiglianza elevata a c_pow
// (target += c_pow * lpdf), per le celle a c fisso e come base della
// ripesatura per importance sampling. Tutto il resto e' identico al
// file base.
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
  real<lower=0> sigma;             // scala base errori
  real<lower=0,upper=0.95> phi;    // autocorrelazione Sigma strutturata
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

  // Sigma strutturata sqrt(ij) — identica all'originale
  matrix[H1, H1] Sigma_h;
  for (i in 1:H1) {
    for (j in 1:H1) {
      Sigma_h[i, j] = square(sigma) * sqrt(1.0 * i * j) / (1 - square(phi))
                      * pow(phi, 1.0 * abs(i - j));
    }
  }
  matrix[H1, H1] L_Sigma = cholesky_decompose(Sigma_h);
}

model {
  // Verosimiglianza (pseudo-): righe indipendenti per costruzione del criterio
  for (t in 1:TT) {
    target += c_pow * multi_normal_cholesky_lpdf(y_mat[t]' | mu_mat[t]', L_Sigma);
  }

  // Prior — cap. 4.5
  gamma_raw ~ std_normal();
  delta ~ beta(2, 2);
  phi ~ beta(2, 4);                       // ⚠️ troncata: vedi header
  log_sigma_h ~ normal(log(0.05), 1);
  sigma ~ student_t(3, 0, 1);
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
