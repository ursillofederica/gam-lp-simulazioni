# Giro rifit divergenze del 22/7: protocollo div>10 -> adapt_delta 0.99,
# esteso a TUTTE le celle MC del disegno (14). Salta le repliche gia'
# rifittate (se persiste: riportata) e la superficie (esclusa per
# decisione del 21/7: reperto, non replica MC). Niente unlink.
# eseguire dalla cartella degli script
source("02_fitta.R")
bk <- "output/fit/_backup_div_22lug"; dir.create(bk, showWarnings = FALSE)
SPEC <- list(
  list("riferimento","lkj",NA),        list("riferimento","strutturata",NA),
  list("persistente","lkj",NA),        list("persistente","strutturata",NA),
  list("persistente","lkj_c80",0.8),   list("persistente","lkj_c90",0.9),
  list("persistente","strutturata_c80",0.8),
  list("persistente","lkj_c70",0.7),   list("persistente","strutturata_c70",0.7),
  list("nonlineare","lkj",NA),         list("nonlineare","strutturata",NA),
  list("nonlineare","lkj_c80",0.8),    list("nonlineare","lkj_c70",0.7),
  list("nonlineare","strutturata_c70",0.7))
n <- 0; persistenti <- character()
for (sp in SPEC) {
  dgp <- sp[[1]]; cella <- sp[[2]]; pow <- sp[[3]]
  for (fp in list.files(file.path("output/fit", dgp, cella), pattern="^rep_", full.names=TRUE)) {
    f <- readRDS(fp)
    if (f$diagn$divergenze <= 10) next
    if (isTRUE(f$diagn$rifit)) {
      persistenti <- c(persistenti, sprintf("%s/%s rep %03d (div %d)", dgp, cella, f$r, f$diagn$divergenze)); next
    }
    file.copy(fp, file.path(bk, sprintf("%s_%s_rep_%03d.rds", dgp, cella, f$r)))
    out <- if (is.na(pow)) fitta_replica(dgp, f$r, cella, adapt_delta = 0.99, seed_stan = 1499)
           else fitta_replica_pow(dgp, f$r, sub("_c[0-9]+$","",cella), pow, adapt_delta = 0.99, seed_stan = 7500)
    out$diagn$rifit <- TRUE
    saveRDS(out, fp); n <- n + 1
    cat(sprintf("%s/%s rep %03d | rhat %.3f | div %d\n", dgp, cella, f$r, out$diagn$rhat_max, out$diagn$divergenze))
  }
}
cat(sprintf("RIFIT DIV 22/7 COMPLETATI: %d repliche\n", n))
if (length(persistenti)) { cat("PERSISTENTI (riportate):\n"); cat(paste(" -", persistenti, collapse="\n"), "\n") } else cat("Nessuna persistente.\n")
