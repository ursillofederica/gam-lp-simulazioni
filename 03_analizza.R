# 03_analizza.R: ripesatura (PSIS) e metriche di copertura
#
# Ripesatura per importance sampling via loo::psis (Vehtari, Gelman, Gabry;
# PSIS: Vehtari et al. 2024, JMLR). Bande via Hmisc::wtd.quantile (Harrell).
# Il calcolo delicato (stabilizzazione, lisciatura delle code, quantile
# ponderato) è nei pacchetti; codice originale = solo le definizioni di
# copertura/bias/ampiezza (Morris, White & Crowther 2019, Tab. 6).

suppressMessages({ library(loo); library(Hmisc) })

GRIGLIA_C <- c(0.7, 0.8, 0.9, 1)
LIVELLO   <- c(0.05, 0.95)                 # bande equal-tailed 90%

# Pesi per un dato c: PSIS su (c-base_c)*log_g -> pesi + Pareto-k + ESS
# base_c = potenza a cui sono stati campionati i fit (1 per le celle
# ordinarie; 0.8/0.9/... per le celle a c fisso, dedotta da c_pow).
pesi_psis <- function(log_g, c, base_c = 1) {
    if (c == base_c) {                     # base: nessuna ripesatura, pesi uniformi
        S <- length(log_g)
        return(list(w = rep(1/S, S), k = NA_real_, ess = S))
    }
    lr <- (c - base_c) * log_g             # log importance ratios
    ps <- suppressWarnings(loo::psis(matrix(lr, ncol = 1), r_eff = 1))
    w  <- as.vector(weights(ps, log = FALSE))
    list(w = w, k = loo::pareto_k_values(ps), ess = 1 / sum(w^2))
}

# Metriche di una serie di draw contro un vero (le 3 def. del 5.1)
metriche_serie <- function(draws, vera, w) {
    H1 <- ncol(draws); lo <- md <- hi <- numeric(H1)
    for (h in 1:H1) {
        q <- Hmisc::wtd.quantile(draws[, h], weights = w,
                                 probs = c(LIVELLO[1], 0.5, LIVELLO[2]), normwt = TRUE)
        lo[h] <- q[1]; md[h] <- q[2]; hi[h] <- q[3]
    }
    list(dentro   = vera >= lo & vera <= hi,    # copertura (Morris et al. 2019)
         bias     = md - vera,                  # bias
         ampiezza = hi - lo)                    # ampiezza
}

# Metriche di una replica a un dato c
# Un contrasto (E1/E2): irf_vera vettore, irf_draws matrice.
# Due contrasti (GHKP, tanh): irf_vera lista pos/neg. Il modello lineare produce
# il solo contrasto unitario, si scala per delta; la superficie dà già pos/neg.
metriche_replica <- function(irf_draws, irf_vera, log_g, c, deltas = NULL,
                             base_c = 1) {
    pw <- pesi_psis(log_g, c, base_c)
    if (is.list(irf_vera)) {                          # due contrasti
        if (is.list(irf_draws)) {                     # superficie: gia' pos/neg
            dp <- irf_draws$pos; dn <- irf_draws$neg
        } else {                                      # lineare: scala l'unitario
            dp <- deltas[1] * irf_draws; dn <- deltas[2] * irf_draws
        }
        return(list(pos = metriche_serie(dp, irf_vera$pos, pw$w),
                    neg = metriche_serie(dn, irf_vera$neg, pw$w),
                    k = pw$k, ess_rel = pw$ess / length(pw$w)))
    }
    c(metriche_serie(irf_draws, irf_vera, pw$w),      # un contrasto
      list(k = pw$k, ess_rel = pw$ess / length(pw$w)))
}

# Delta dei contrasti per DGP (per scalare il lineare a due contrasti) 
DELTAS <- list(nonlineare = c(2, -2), tanh = c(2, 1))  

# Aggrega le metriche di una serie (pos/neg o singola) sulle repliche
aggrega_serie <- function(serie_list) {  
    list(copertura_h = rowMeans(sapply(serie_list, `[[`, "dentro")),
         bias_h      = rowMeans(sapply(serie_list, `[[`, "bias")),
         ampiezza_h  = apply(sapply(serie_list, `[[`, "ampiezza"), 1, median))
}

# Aggregazione su una cella (dgp x modello), tutta la griglia di c 
analizza_cella <- function(dgp, modello, c_grid = GRIGLIA_C) {
    files <- list.files(file.path("output", "fit", dgp, modello),
                        pattern = "^rep_", full.names = TRUE)
    stopifnot(length(files) > 0)
    fits  <- lapply(files, readRDS)
    due   <- is.list(fits[[1]]$irf_vera)             
    deltas <- DELTAS[[dgp]]
    # base_c dedotta dai fit stessi (c_pow salvato da fitta_replica_pow):
    # le celle campionate a c fisso vengono ripesate dalla loro base vera
    base_c <- if (!is.null(fits[[1]]$c_pow)) fits[[1]]$c_pow else 1
    out <- list()
    for (c in c_grid) {
        m <- lapply(fits, function(f)
            metriche_replica(f$irf_draws, f$irf_vera, f$log_g, c, deltas,
                             base_c = base_c))
        k  <- suppressWarnings(max(sapply(m, `[[`, "k"), na.rm = TRUE))
        er <- median(sapply(m, `[[`, "ess_rel"))
        if (due) {                                    # aggrega pos e neg separati
            out[[as.character(c)]] <- list(
                pos = aggrega_serie(lapply(m, `[[`, "pos")),
                neg = aggrega_serie(lapply(m, `[[`, "neg")),
                pareto_k_max = k, ess_rel_med = er, R = length(fits),
                base_c = base_c)
        } else {
            out[[as.character(c)]] <- c(
                aggrega_serie(m),
                list(pareto_k_max = k, ess_rel_med = er, R = length(fits),
                     base_c = base_c))
        }
    }
    out
}

stampa_cella <- function(res, dgp, modello) {
    cat("====", dgp, "x", modello, "| R =", res[[1]]$R,
        "| base c =", res[[1]]$base_c, "====\n")
    for (c in names(res)) {
        r <- res[[c]]
        kflag <- if (as.numeric(c) == r$base_c) " [fit diretti, non ripesati]"
                 else if (is.finite(r$pareto_k_max) && r$pareto_k_max > 0.7)
                     " [k>0.7: fit diretto]" else ""
        cat(sprintf("c = %-4s | copertura %.3f | ampiezza %.3f | Pareto-k max %.2f%s\n",
                    c, mean(r$copertura_h), mean(r$ampiezza_h),
                    r$pareto_k_max, kflag))
    }
}
