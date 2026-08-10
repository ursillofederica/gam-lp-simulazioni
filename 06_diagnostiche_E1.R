# Diagnostiche E1: certificazione di convergenza, affidabilita' della
# ripesatura, calibrazione estesa (errore MC esatto, RMSE, copertura
# congiunta). Riusa pesi_psis e le librerie di 03_analizza.R (loo, Hmisc).
# Eseguire dalla cartella degli script.
source("03_analizza.R")

INIT_REFIT <- list(lkj = c(8, 61, 65, 69), strutturata = integer(0))

out <- list()
for (m in c("lkj", "strutturata")) {
  fits <- lapply(list.files(file.path("output/fit/riferimento", m),
                            pattern = "^rep_", full.names = TRUE), readRDS)
  R <- length(fits)

  # convergenza
  rhat <- sapply(fits, function(f) f$diagn$rhat_max)
  div  <- sapply(fits, function(f) f$diagn$divergenze)
  essb <- sapply(fits, function(f) f$diagn$ess_bulk_min)
  n_div_refit <- sum(sapply(fits, function(f) isTRUE(f$diagn$rifit)))
  bk <- list.files("output/fit/_backup_div_21lug",
                   pattern = paste0("^riferimento_", m), full.names = TRUE)
  div_pre_max <- if (length(bk)) max(sapply(bk, function(p) readRDS(p)$diagn$divergenze)) else NA
  cert <- data.frame(R = R, rifit_init = length(INIT_REFIT[[m]]),
                     rifit_div = n_div_refit, div_max_pre = div_pre_max,
                     rhat_max = max(rhat), div_max_post = max(div),
                     ess_min = round(min(essb)), quota_div0 = mean(div == 0))

  # ripesatura
  sd_logg <- sapply(fits, function(f) sd(f$log_g))
  rip <- do.call(rbind, lapply(c(0.9, 0.8, 0.7), function(cc) {
    pw <- lapply(fits, function(f) pesi_psis(f$log_g, cc))
    k  <- sapply(pw, `[[`, "k"); er <- sapply(pw, function(p) p$ess / length(p$w))
    data.frame(c = cc, k_med = median(k), k_max = max(k),
               quota_k_alto = mean(k > 0.7), ess_rel_med = median(er))
  }))

  # calibrazone estesa
  cal <- lapply(c("1" = 1, "0.9" = 0.9), function(cc) {
    per_rep <- lapply(fits, function(f) {
      w <- pesi_psis(f$log_g, cc)$w
      q <- apply(f$irf_draws, 2, function(col)
        Hmisc::wtd.quantile(col, weights = w, probs = c(0.05, 0.5, 0.95), normwt = TRUE))
      list(dentro = f$irf_vera >= q[1, ] & f$irf_vera <= q[3, ],
           err = q[2, ] - f$irf_vera, amp = q[3, ] - q[1, ])
    })
    dentro <- sapply(per_rep, `[[`, "dentro")            # 17 x R
    err    <- sapply(per_rep, `[[`, "err")
    amp    <- sapply(per_rep, `[[`, "amp")
    list(cop_h  = rowMeans(dentro),
         mcse_h = sqrt(rowMeans(dentro) * (1 - rowMeans(dentro)) / R),
         rmse_h = sqrt(rowMeans(err^2)),
         bias_h = rowMeans(err),
         amp_h  = apply(amp, 1, median),
         cop_congiunta = mean(colSums(dentro) == nrow(dentro)))
  })
    
  k_per_c <- sapply(c("0.9" = 0.9, "0.8" = 0.8, "0.7" = 0.7), function(cc)
    sapply(fits, function(f) pesi_psis(f$log_g, cc)$k))
  out[[m]] <- list(cert = cert, sd_logg_med = median(sd_logg),
                   ripesatura = rip, calibrazione = cal,
                   repliche = data.frame(rhat = rhat, div = div, ess = essb,
                                         sd_logg = sd_logg),
                   k_repliche = k_per_c)

  cat("====", m, "====\n"); print(cert, row.names = FALSE)
  cat("sd(log g) mediana:", round(median(sd_logg), 2), "\n")
  print(rip, row.names = FALSE)
  cat(sprintf("c=1  : cop media %.3f | congiunta %.2f | RMSE medio %.3f\n",
              mean(out[[m]]$calibrazione$`1`$cop_h),
              out[[m]]$calibrazione$`1`$cop_congiunta,
              mean(out[[m]]$calibrazione$`1`$rmse_h)))
  cat(sprintf("c=0.9: cop media %.3f | congiunta %.2f | RMSE medio %.3f\n",
              mean(out[[m]]$calibrazione$`0.9`$cop_h),
              out[[m]]$calibrazione$`0.9`$cop_congiunta,
              mean(out[[m]]$calibrazione$`0.9`$rmse_h)))
}
saveRDS(out, "output/diagnostiche_E1.rds")
