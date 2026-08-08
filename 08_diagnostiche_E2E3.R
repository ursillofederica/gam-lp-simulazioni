# Diagnostiche e calibrazione estesa per E2 (persistente) ed E3 (nonlineare):
# Output: output/diagnostiche_E2E3.rds. Eseguire dalla cartella degli script.
source("03_analizza.R")

INIT_REFIT <- list(
  "persistente/lkj"     = c(8, 22, 29, 39),
  "persistente/lkj_c80" = c(5, 30, 33, 38, 42),
  "persistente/lkj_c90" = c(5, 13, 30, 38, 39, 41, 42, 50),
  "nonlineare/lkj"      = 9)

CELLE <- list(
  E2 = list(dgp = "persistente",
            celle = c("lkj","strutturata","lkj_c80","lkj_c90",
                      "strutturata_c80","lkj_c70","strutturata_c70")),
  E3 = list(dgp = "nonlineare",
            celle = c("lkj","strutturata","lkj_c80","lkj_c70","strutturata_c70")))
GRIGLIA <- c(0.7, 0.8, 0.9, 1)
BK <- c("output/fit/_backup_div_21lug", "output/fit/_backup_div_22lug")

calibra <- function(fits, cc, base_c, deltas = NULL) {
  due <- is.list(fits[[1]]$irf_vera)
  per_serie <- function(get_draws, get_vera) {
    per_rep <- lapply(fits, function(f) {
      w <- pesi_psis(f$log_g, cc, base_c)$w
      q <- apply(get_draws(f), 2, function(col)
        Hmisc::wtd.quantile(col, weights = w, probs = c(0.05, 0.5, 0.95), normwt = TRUE))
      v <- get_vera(f)
      list(dentro = v >= q[1, ] & v <= q[3, ], err = q[2, ] - v, amp = q[3, ] - q[1, ])
    })
    dentro <- sapply(per_rep, `[[`, "dentro"); err <- sapply(per_rep, `[[`, "err")
    amp <- sapply(per_rep, `[[`, "amp"); R <- length(per_rep)
    list(cop_h = rowMeans(dentro),
         mcse_h = sqrt(rowMeans(dentro) * (1 - rowMeans(dentro)) / R),
         bias_h = rowMeans(err), rmse_h = sqrt(rowMeans(err^2)),
         amp_h = apply(amp, 1, median))
  }
  if (due) list(pos = per_serie(function(f) deltas[1] * f$irf_draws,
                                function(f) f$irf_vera$pos),
                neg = per_serie(function(f) deltas[2] * f$irf_draws,
                                function(f) f$irf_vera$neg))
  else per_serie(function(f) f$irf_draws, function(f) f$irf_vera)
}

out <- list()
for (esp in names(CELLE)) {
  dgp <- CELLE[[esp]]$dgp
  for (m in CELLE[[esp]]$celle) {
    fits <- lapply(list.files(file.path("output/fit", dgp, m),
                              pattern = "^rep_", full.names = TRUE), readRDS)
    R <- length(fits)
    base_c <- if (!is.null(fits[[1]]$c_pow)) fits[[1]]$c_pow else 1
    deltas <- if (dgp == "nonlineare") c(2, -2) else NULL

    rhat <- sapply(fits, function(f) f$diagn$rhat_max)
    div  <- sapply(fits, function(f) f$diagn$divergenze)
    essb <- sapply(fits, function(f) f$diagn$ess_bulk_min)
    bkf <- unlist(lapply(BK, list.files, pattern = paste0("^", dgp, "_", m, "_rep"),
                         full.names = TRUE))
    cert <- data.frame(R = R, base_c = base_c,
                       rifit_init = length(INIT_REFIT[[paste(dgp, m, sep = "/")]]),
                       rifit_div = sum(sapply(fits, function(f) isTRUE(f$diagn$rifit))),
                       div_max_pre = if (length(bkf)) max(sapply(bkf, function(p) readRDS(p)$diagn$divergenze)) else NA,
                       rhat_max = max(rhat), div_max_post = max(div),
                       ess_min = round(min(essb)), quota_div0 = mean(div == 0))

    rip <- do.call(rbind, lapply(setdiff(GRIGLIA, base_c), function(cc) {
      k <- sapply(fits, function(f) pesi_psis(f$log_g, cc, base_c)$k)
      er <- sapply(fits, function(f) { p <- pesi_psis(f$log_g, cc, base_c); p$ess / length(p$w) })
      data.frame(c = cc, k_med = median(k), k_max = max(k),
                 quota_k_alto = mean(k > 0.7), ess_rel_med = median(er))
    }))
    k_rep <- sapply(setdiff(GRIGLIA, base_c), function(cc)
      sapply(fits, function(f) pesi_psis(f$log_g, cc, base_c)$k))
    colnames(k_rep) <- setdiff(GRIGLIA, base_c)

    cal <- lapply(setNames(GRIGLIA, GRIGLIA), function(cc)
      calibra(fits, cc, base_c, deltas))

    out[[esp]][[m]] <- list(cert = cert,
                            sd_logg_med = median(sapply(fits, function(f) sd(f$log_g))),
                            ripesatura = rip, k_repliche = k_rep,
                            repliche = data.frame(rhat = rhat, div = div, ess = essb),
                            calibrazione = cal)
  }
}
saveRDS(out, "output/diagnostiche_E2E3.rds")
