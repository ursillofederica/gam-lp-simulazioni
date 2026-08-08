# 02_fitta.R: stima dei modelli sulle repliche congelate
#
# Legge le repliche da output/dati/<dgp>/, costruisce basi e stan_data
# (contratto di 01_genera.R), stima con cmdstanr e salva
# in output/fit/<dgp>/<modello>/rep_<r>.rds:
#   draw della IRF, log g per draw (per la ripesatura), diagnostiche.
#
# Specifiche: 4 catene, 1000 warmup + 2000 sampling, adapt_delta 0.95.
# Base su h: B-spline cubica K = 8, nodi interni equispaziati

library(cmdstanr)
library(splines)

source("01_genera.R")     

# Base sull'orizzonte: K = 8, cubica, nodi interni equispaziati 
K_H     <- 8
h_grid  <- 0:H
n_nodi  <- K_H - 4                                   # cubica con intercetta: df = nodi + 4
nodi_h  <- seq(0, H, length.out = n_nodi + 2)[-c(1, n_nodi + 2)]  # interni, equispaziati
B_H     <- bs(h_grid, knots = nodi_h, degree = 3, intercept = TRUE)
stopifnot(ncol(B_H) == K_H)

# stan_data 
# Ordinamento del vettore impilato: h esterno, t interno (come nei .stan:
# mu_vec[(h-1)*TT + t]). TP = kronecker(B_H, s)
costruisci_stan_data <- function(d, sigma_theta = 1) {
    TP <- kronecker(as.matrix(B_H), matrix(d$s, ncol = 1))
    list(TT = TT, K = K_H, H1 = H1, N = TT * H1,
         K_ctrl = ncol(d$X_ctrl),
         TP = TP, y_mat = d$y_mat, X_ctrl = d$X_ctrl,
         TP_s1 = as.matrix(B_H),         
         sigma_theta = sigma_theta)
}

# Base sullo shock (superficie, Processo 3): J = 6
J_S     <- 6
n_nodi_s <- J_S - 4 
nodi_s  <- seq(-4, 4, length.out = n_nodi_s + 2)[-c(1, n_nodi_s + 2)]

# stan_data per la superficie
# TP[(h-1)*TT + t, (j-1)*K + k] = B_s[t, j] * B_H[h, k]: colonna per colonna
# via prodotto di Kronecker vettoriale (h esterno, t interno come nei .stan).
costruisci_stan_data_superficie <- function(d, sigma_theta = 1, deltas = c(2, -2)) {
    # deltas = c(a, b): irf_pos = contrasto a delta a, irf_neg a delta b.
    # GHKP: c(2, -2) segno; tanh (saturazione, dispari): c(2, 1) taglia.
    B_S     <- bs(d$s, knots = nodi_s, degree = 3, intercept = TRUE,
                  Boundary.knots = c(-4, 4))
    b_pos  <- predict(B_S, deltas[1]); b_neg <- predict(B_S, deltas[2]); b_zero <- predict(B_S, 0)
    TP <- matrix(0, TT * H1, J_S * K_H)
    TP_pos <- TP_neg <- TP_zero <- matrix(0, H1, J_S * K_H)
    for (j in 1:J_S) for (k in 1:K_H) {
        col <- (j - 1) * K_H + k
        TP[, col]      <- kronecker(B_H[, k], B_S[, j])
        TP_pos[, col]  <- b_pos[j]  * B_H[, k]
        TP_neg[, col]  <- b_neg[j]  * B_H[, k]
        TP_zero[, col] <- b_zero[j] * B_H[, k]
    }
    list(TT = TT, K = K_H, J = J_S, H1 = H1, N = TT * H1,
         K_ctrl = ncol(d$X_ctrl),
         TP = TP, y_mat = d$y_mat, X_ctrl = d$X_ctrl,
         TP_pos = TP_pos, TP_neg = TP_neg, TP_zero = TP_zero,
         sigma_theta = sigma_theta)
}

MODELLI <- list(
    strutturata            = cmdstan_model(file.path("stan", "lp_lineare_strutturata.stan")),
    lkj                    = cmdstan_model(file.path("stan", "lp_lineare_lkj.stan")),
    superficie_strutturata = cmdstan_model(file.path("stan", "lp_superficie.stan")),
    superficie_lkj         = cmdstan_model(file.path("stan", "lp_superficie_lkj.stan"))
)

# Fit di una replica 
fitta_replica <- function(dgp, r, modello,
                          adapt_delta = 0.95, seed_stan = 999) {
    d  <- readRDS(file.path(DIR_OUT, dgp, sprintf("rep_%03d.rds", r)))
    superficie <- startsWith(modello, "superficie")
    deltas <- if (dgp == "tanh") c(2, 1) else c(2, -2)  
    sd_stan <- if (superficie) costruisci_stan_data_superficie(d, deltas = deltas)
               else            costruisci_stan_data(d)

    fit <- MODELLI[[modello]]$sample(
        data = sd_stan, chains = 4, parallel_chains = 4,
        iter_warmup = 1000, iter_sampling = 2000,
        adapt_delta = adapt_delta, seed = seed_stan + r,
        init = 0.1,   # init per evitare un Cholesky quasi singolare all'avvio (catena bloccata in partenza)
        refresh = 0, show_messages = FALSE)

    if (superficie) {
        irf_draws <- list(pos = fit$draws("irf_pos", format = "matrix"),
                          neg = fit$draws("irf_neg", format = "matrix"))
        riass <- fit$summary(variables = c("irf_pos", "irf_neg"))
    } else {
        irf_draws <- fit$draws("irf", format = "matrix")      # S x H1
        riass     <- fit$summary(variables = c("irf"))
    }
    log_g <- rowSums(fit$draws("log_lik", format = "matrix"))  # S: log g per draw
    diagn <- list(
        rhat_max    = max(riass$rhat, na.rm = TRUE),
        divergenze  = sum(fit$diagnostic_summary()$num_divergent),
        ess_bulk_min = min(riass$ess_bulk, na.rm = TRUE),
        rifit = FALSE)

    # Politica repliche fallite: un rifit con adapt_delta alzato;
    if (diagn$rhat_max > 1.01 && adapt_delta < 0.99) {
        out2 <- fitta_replica(dgp, r, modello,
                              adapt_delta = 0.99, seed_stan = seed_stan + 500)
        out2$diagn$rifit <- TRUE
        saveRDS(out2, file.path("output", "fit", dgp, modello,
                                sprintf("rep_%03d.rds", r)))   
        return(invisible(out2))
    }

    out <- list(dgp = dgp, r = r, modello = modello,
                irf_draws = irf_draws, log_g = log_g,
                irf_vera = d$irf_vera, diagn = diagn)
    dir_fit <- file.path("output", "fit", dgp, modello)
    dir.create(dir_fit, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file.path(dir_fit, sprintf("rep_%03d.rds", r)))
    invisible(out)
}

# Ciclo su un blocco di repliche
fitta_blocco <- function(dgp, modello, repliche) {
    for (r in repliche) {
        percorso <- file.path("output", "fit", dgp, modello,
                              sprintf("rep_%03d.rds", r))
        if (file.exists(percorso)) next
        out <- fitta_replica(dgp, r, modello)
        cat(sprintf("%s | %s | rep %03d | rhat %.3f | div %d\n",
                    dgp, modello, r, out$diagn$rhat_max, out$diagn$divergenze))
    }
}

# Fit diretti alla potenza c
MODELLI_POW <- new.env()
modello_pow <- function(modello) {
    if (is.null(MODELLI_POW[[modello]]))
        MODELLI_POW[[modello]] <- cmdstan_model(
            file.path("stan", paste0("lp_lineare_", modello, "_pow.stan")))
    MODELLI_POW[[modello]]
}

fitta_replica_pow <- function(dgp, r, modello, c_pow,
                              adapt_delta = 0.95, seed_stan = 7000) {
    d <- readRDS(file.path(DIR_OUT, dgp, sprintf("rep_%03d.rds", r)))
    sd_stan <- costruisci_stan_data(d)
    sd_stan$c_pow <- c_pow
    fit <- modello_pow(modello)$sample(
        data = sd_stan, chains = 4, parallel_chains = 4,
        iter_warmup = 1000, iter_sampling = 2000,
        adapt_delta = adapt_delta, seed = seed_stan + r,
        init = 0.1,  
        refresh = 0, show_messages = FALSE)
    irf_draws <- fit$draws("irf", format = "matrix")
    riass <- fit$summary(variables = "irf")
    out <- list(dgp = dgp, r = r, modello = modello, c_pow = c_pow,
                irf_draws = irf_draws,
                log_g = rowSums(fit$draws("log_lik", format = "matrix")),
                irf_vera = d$irf_vera,
                diagn = list(rhat_max = max(riass$rhat, na.rm = TRUE),
                             divergenze = sum(fit$diagnostic_summary()$num_divergent),
                             ess_bulk_min = min(riass$ess_bulk, na.rm = TRUE)))
    dir_fit <- file.path("output", "fit", dgp,
                         sprintf("%s_c%02d", modello, round(100 * c_pow)))
    dir.create(dir_fit, recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, file.path(dir_fit, sprintf("rep_%03d.rds", r)))
    invisible(out)
}

fitta_blocco_pow <- function(dgp, modello, repliche, c_pow) {
    for (r in repliche) {
        percorso <- file.path("output", "fit", dgp,
                              sprintf("%s_c%02d", modello, round(100 * c_pow)),
                              sprintf("rep_%03d.rds", r))
        if (file.exists(percorso)) next
        out <- fitta_replica_pow(dgp, r, modello, c_pow)
        cat(sprintf("%s | %s c=%.1f | rep %03d | rhat %.3f | div %d\n",
                    dgp, modello, c_pow, r,
                    out$diagn$rhat_max, out$diagn$divergenze))
    }
}

