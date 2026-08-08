# Rifit 23/7: repliche con rhat > 1.01 sopravvissute a tutti i giri precedenti.
# Trovate dalla certificazione E2/E3 (08_diagnostiche_E2E3.R): entrambe in
# persistente/lkj, entrambe fit ORIGINALI del 19-20/7, cioe' PRIMA del fix
# init = 0.1 del 21/7. Non erano state prese perche':
#   - la scansione del 21/7 cercava la firma sd(log g) ~ 1e10-1e12 (qui 3.6e3);
#   - il giro divergenze filtra su div > 10 (qui div = 0).
# rep 067: catena 2 bloccata a log g -11256 contro -3011 delle altre tre
#          (rhat 1.578, ess_bulk 7). rep 083: rhat 1.015, ess_bulk 338.
# Protocollo (§ 5.1): rifit, mai esclusione. Backup prima di sovrascrivere.
# eseguire dalla cartella degli script
source("02_fitta.R")

bk <- "output/fit/_backup_rhat_23lug"; dir.create(bk, showWarnings = FALSE)
DA_RIFARE <- list(list("persistente", "lkj", 67), list("persistente", "lkj", 83))

for (sp in DA_RIFARE) {
  dgp <- sp[[1]]; cella <- sp[[2]]; r <- sp[[3]]
  fp <- file.path("output/fit", dgp, cella, sprintf("rep_%03d.rds", r))
  vecchio <- readRDS(fp)
  file.copy(fp, file.path(bk, sprintf("%s_%s_rep_%03d.rds", dgp, cella, r)),
            overwrite = TRUE)
  cat(sprintf("\n== %s/%s rep %03d | PRIMA: rhat %.3f, ess %.0f, div %d\n",
              dgp, cella, r, vecchio$diagn$rhat_max,
              vecchio$diagn$ess_bulk_min, vecchio$diagn$divergenze))
  out <- fitta_replica(dgp, r, cella, adapt_delta = 0.95, seed_stan = 1499)
  out$diagn$rifit <- TRUE
  saveRDS(out, fp)
  S <- length(out$log_g); per <- S / 4
  meds <- sapply(1:4, function(i) mean(out$log_g[((i - 1) * per + 1):(i * per)]))
  cat(sprintf("   DOPO: rhat %.3f, ess %.0f, div %d | log g per catena: %s\n",
              out$diagn$rhat_max, out$diagn$ess_bulk_min, out$diagn$divergenze,
              paste(sprintf("%.1f", meds), collapse = " ")))
}
cat("\nRIFIT RHAT 23/7 COMPLETATI\n")
