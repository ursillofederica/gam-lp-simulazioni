# ==============================================================================
# 01_genera.R — Cap. 5: generazione dei dataset simulati (protocollo § 5.1)
#
# Nota numerazione: capitoli e paragrafi citati nel codice seguono la
# stesura interna; nel PDF della tesi scalano di uno (il capitolo dei
# risultati, qui cap. 5, nel PDF e' il 4: l'introduzione non e' numerata).
#
# Contratto: ogni genera_<nome>(r) restituisce
#   list(y_mat [T x H1], X_ctrl [T x K_ctrl], s [T], irf_vera, info)
# Seed: SEED_BASE[dgp] + r. Ogni replica salvata in output/dati/<dgp>/rep_<r>.rds
# e MAI rigenerata se il file esiste (congelamento, § 5.1).
# ==============================================================================

# ---- Costanti del protocollo (§ 5.1) ----------------------------------------
TT      <- 100                 # righe del sistema
H       <- 16                  # orizzonte massimo
H1      <- H + 1               # numero orizzonti (h = 0..16)
BURN    <- 1000                # burn-in traiettorie
R_MAIN  <- 100                 # repliche per processo
R_PILOT <- 20                  # repliche cancello pilota

SEED_BASE <- c(pilota = 1000, riferimento = 2000,
               persistente = 3000, nonlineare = 4000)

DIR_OUT <- file.path("output", "dati")

# ---- Sigma0 strutturata (pilota): sigma0 = 0.5, phi0 = 0.6 ------------------
# Entrambe le parametrizzazioni la contengono -> il pilota valida entrambe.
costruisci_sigma_strutturata <- function(sigma0, phi0, H1) {
    idx <- 1:H1
    outer(idx, idx, function(i, j)
        sigma0^2 * sqrt(i * j) / (1 - phi0^2) * phi0^abs(i - j))
}
SIGMA0 <- costruisci_sigma_strutturata(sigma0 = 0.5, phi0 = 0.6, H1 = H1)

# ---- Helper: righe sovrapposte da una traiettoria ---------------------------
# Da y (lunghezza >= inizio + TT + H) costruisce y_mat[t, h+1] = y_{inizio+t+h}.
costruisci_righe <- function(y, inizio) {
    y_mat <- matrix(NA_real_, TT, H1)
    for (t in 1:TT) y_mat[t, ] <- y[(inizio + t - 1) + 0:H]
    y_mat
}

# ==============================================================================
# PILOTA — righe indipendenti ben specificate (fase 1 + errori da SIGMA0)
# Ruolo: criterio = verosimiglianza vera -> copertura nominale o e' un bug.
# ==============================================================================
genera_pilota <- function(r) {
    set.seed(SEED_BASE["pilota"] + r)
    h_grid <- 0:H
    theta0 <- 0.5 * h_grid * exp(-0.2 * h_grid)      # gobba di fase 1
    s      <- rnorm(TT)
    V      <- matrix(rnorm(TT * H1), TT, H1) %*% chol(SIGMA0)
    y_mat  <- outer(s, theta0) + V
    # X_ctrl = sola intercetta (alpha_h della spec LP, cap. 2):
    # il DGP ha media zero ma il modello non ha intercetta propria.
    X_ctrl <- matrix(1, TT, 1)
    list(y_mat = y_mat, X_ctrl = X_ctrl, s = s,
         irf_vera = theta0,
         info = list(dgp = "pilota", r = r,
                     seed = SEED_BASE["pilota"] + r,
                     sigma0 = 0.5, phi0 = 0.6))
}

# ==============================================================================
# PROCESSO 1 — riferimento, spirito HMP (eq:dgp_riferimento)
#   y_t = sum_{p=1..5} phi_p y_{t-p} + eps_t + alpha*T^(-1/2) sum_j a_j eps_{t-j}
#   eps_t = s_t + u_t;  alpha = 2;  a_j ~ N(0,1) iid PER REPLICA (j <= 10)
#   phi_p: FISSI, da calibrazione AR(5) su crescita output USA (phi_calibrati.rds)
#   Controlli: intercetta + 5 ritardi di y (liberi per orizzonte, nel modello)
#   IRF vera: ricorsione psi (eq:irf_ricorsione), replica per replica
# ------------------------------------------------------------------------------
# Calibrazione: GDPC1 log-diff x100, 1960-2019Q4 (00_calibrazione_riferimento.R)
# ==============================================================================
genera_riferimento <- function(r) {
    set.seed(SEED_BASE["riferimento"] + r)
    phi   <- readRDS(file.path("output", "phi_calibrati.rds"))$phi   # FISSI
    alpha <- 2
    peso  <- alpha * TT^(-1/2)          # il peso MA: alpha * T^(-1/2)
    a     <- rnorm(10)                  # coefficienti MA: PER REPLICA

    # Traiettoria: burn-in + 5 (ritardi controlli) + TT + H
    ntot <- BURN + 5 + TT + H
    s_full <- rnorm(ntot)
    u_full <- rnorm(ntot)
    eps    <- s_full + u_full           # eps_t = s_t + u_t
    y <- numeric(ntot)
    for (t in 11:ntot) {                # parte da 11: servono 5 lag y, 10 lag eps
        y[t] <- sum(phi * y[t - (1:5)]) + eps[t] + peso * sum(a * eps[t - (1:10)])
    }

    # Righe del sistema: la riga t usa la traiettoria da 'inizio + t - 1'
    inizio <- BURN + 6
    y_mat  <- costruisci_righe(y, inizio)
    s      <- s_full[inizio + 0:(TT - 1)]

    # Controlli: intercetta + 5 ritardi di y (non standardizzati: scale ~ 1;
    # la scala del prior si fissa nel driver)
    ritardi <- sapply(1:5, function(p) y[(inizio - p) + 0:(TT - 1)])
    X_ctrl  <- cbind(1, ritardi)

    # IRF vera: ricorsione psi (eq:irf_ricorsione), dipende dagli a_j pescati
    psi <- numeric(H1)                  # psi[h+1] = psi_h
    psi[1] <- 1
    for (h in 1:H) {
        ar_part <- sum(phi[1:min(h, 5)] * psi[h + 1 - (1:min(h, 5))])
        ma_part <- if (h <= 10) peso * a[h] else 0
        psi[h + 1] <- ar_part + ma_part
    }

    list(y_mat = y_mat, X_ctrl = X_ctrl, s = s,
         irf_vera = psi,
         info = list(dgp = "riferimento", r = r,
                     seed = SEED_BASE["riferimento"] + r,
                     phi = phi, a = a, alpha = alpha))
}

# ==============================================================================
# PROCESSO 2 — persistente, Esempio 2.1 (eq:dgp_persistente)
#   y_t = rho y_{t-1} + eps_t;  eps_t = s_t + u_t;  rho = 0.9
#   Controlli: intercetta + y_{t-1};  IRF vera: theta_h = rho^h
# ==============================================================================
genera_persistente <- function(r) {
    set.seed(SEED_BASE["persistente"] + r)
    rho  <- 0.9
    ntot <- BURN + 1 + TT + H           # un solo ritardo per i controlli

    s_full <- rnorm(ntot)
    u_full <- rnorm(ntot)
    eps    <- s_full + u_full           # eps_t = s_t + u_t

    y <- numeric(ntot)                  # allocato: parte da 0, il burn-in dimentica
    for (t in 2:ntot) {
        y[t] <- rho * y[t-1] + eps[t]
    }

    inizio <- BURN + 2
    y_mat  <- costruisci_righe(y, inizio)
    s      <- s_full[inizio + 0:(TT-1)]
    X_ctrl <- cbind(1, y[(inizio - 1) + 0:(TT-1)])   # intercetta + y_{t-1}

    irf_vera <- rho^(0:H)               # Esempio 2.1: theta_h = rho^h

    list(y_mat = y_mat, X_ctrl = X_ctrl, s = s,
         irf_vera = irf_vera,
         info = list(dgp = "persistente", r = r,
                     seed = SEED_BASE["persistente"] + r,
                     rho = rho))
}

# ==============================================================================
# PROCESSO 3 — non lineare, GHKP Sez. 6.1 (eq:dgp_nonlineare)
#   x_t = eps1_t;  y_t = 0.5 y_{t-1} + 0.5 x_t + 0.3 x_{t-1}
#                        - 0.4 f(x_t) - 0.3 f(x_{t-1}) + eps2_t;  f = pmax(x, 0)
#   IRF vera (eq:irf_vera_nonlin), contrasti delta = +2 e -2
# ------------------------------------------------------------------------------
# Controlli: intercetta + y_{t-1} + x_{t-1}.
# ==============================================================================
irf_ghkp <- function(delta) {
    irf <- numeric(H1)
    irf[1] <- 0.5 * delta - 0.4 * max(delta, 0)
    irf[2] <- 0.5 * irf[1] + 0.3 * delta - 0.3 * max(delta, 0)
    for (h in 2:H) {
        irf[h+1] <- 0.5 * irf[h]
    }
    irf                                  # ultima espressione = valore restituito
}

genera_nonlineare <- function(r) {
    set.seed(SEED_BASE["nonlineare"] + r)
    ntot <- BURN + 1 + TT + H

    x    <- rnorm(ntot)                  # x_t = eps1_t: lo shock osservato
    eps2 <- rnorm(ntot)
    fx   <- pmax(x, 0)                   # f = max{x, 0}, elemento per elemento

    y <- numeric(ntot)
    for (t in 2:ntot) {
        y[t] <- 0.5 * y[t-1] + 0.5 * x[t] +
                0.3 * x[t-1] - 0.4 * fx[t] - 0.3 * fx[t-1] + eps2[t]
    }

    inizio <- BURN + 2
    y_mat  <- costruisci_righe(y, inizio)
    s      <- x[inizio + 0:(TT-1)]
    X_ctrl <- cbind(1,
                    y[(inizio - 1) + 0:(TT-1)],
                    x[(inizio - 1) + 0:(TT-1)])   # intercetta + y_{t-1} + x_{t-1}

    irf_vera <- list(pos = irf_ghkp(2), neg = irf_ghkp(-2))

    list(y_mat = y_mat, X_ctrl = X_ctrl, s = s,
         irf_vera = irf_vera,
         info = list(dgp = "nonlineare", r = r,
                     seed = SEED_BASE["nonlineare"] + r))
}

# ==============================================================================
# Generazione e congelamento su disco
# ==============================================================================
genera_e_salva <- function(dgp, n_rep) {
    gen <- switch(dgp,
                  pilota      = genera_pilota,
                  riferimento = genera_riferimento,
                  persistente = genera_persistente,
                  nonlineare  = genera_nonlineare,
                  tanh        = genera_tanh)
    dir_dgp <- file.path(DIR_OUT, dgp)
    dir.create(dir_dgp, recursive = TRUE, showWarnings = FALSE)
    for (r in 1:n_rep) {
        percorso <- file.path(dir_dgp, sprintf("rep_%03d.rds", r))
        if (!file.exists(percorso)) saveRDS(gen(r), percorso)   # congelati
    }
    message(dgp, ": ", n_rep, " repliche in ", dir_dgp)
}

# ---- Esecuzione (decommentare quando i generatori sono pronti) --------------
# genera_e_salva("pilota",      R_PILOT)
# genera_e_salva("riferimento", R_MAIN)
# genera_e_salva("persistente", R_MAIN)
# genera_e_salva("nonlineare",  R_MAIN)

# ==============================================================================
# TANH — DGP liscio a segnale forte (banco di recupero, § 5.4 fase 2)
# f(s,h) = 1.5*h*exp(-0.2h)*tanh(0.5*s); righe INDIPENDENTI (isola il recupero
# della media), errori da SIGMA0. tanh e' DISPARI -> contrasti di TAGLIA
# delta in {1,2} (la saturazione), non di segno.
# ==============================================================================
SEED_BASE["tanh"] <- 5000
genera_tanh <- function(r) {
    set.seed(SEED_BASE["tanh"] + r)
    ff <- function(sh, h) 1.5 * h * exp(-0.2 * h) * tanh(0.5 * sh)
    h_grid <- 0:H
    s  <- rnorm(TT)
    mu <- sapply(s, function(ss) ff(ss, h_grid))   # H1 x TT
    # Errori IID come fase 2 (sd = sqrt(0.3)): il tanh e' banco di RECUPERO,
    # la struttura d'errore e' irrilevante e IID e' il setup validato in fase 2.
    V  <- matrix(rnorm(TT * H1, 0, sqrt(0.3)), TT, H1)
    y_mat <- t(mu) + V
    X_ctrl <- matrix(1, TT, 1)                     # righe indipendenti: sola intercetta
    list(y_mat = y_mat, X_ctrl = X_ctrl, s = s,
         irf_vera = list(d1 = ff(1, h_grid) - ff(0, h_grid),
                         d2 = ff(2, h_grid) - ff(0, h_grid)),
         info = list(dgp = "tanh", r = r, seed = SEED_BASE["tanh"] + r))
}
