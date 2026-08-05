# Rifit repliche con divergenze > 10: adapt_delta 0.99, seed protocollo (+500).
# Regola (decisione 21/7): "div > 10 -> rifit a 0.99; se persiste, riportata".
# Scansione dinamica dello stato su disco; RIPRENDIBILE: le repliche gia'
# rifittate (diagn$rifit == TRUE) non vengono rifittate una seconda volta,
# qualunque sia il loro numero di divergenze -> se persiste, finisce nel
# rapporto finale, come da regola. (Versione 21/7 sera: il primo giro, da
# scratchpad, e' stato interrotto a riferimento/strutturata rep_051.)
setwd("/Users/federicaursillo/Tesi/01_codice/GAM-LP/cap5")
source("02_fitta.R")
bk <- "output/fit/_backup_div_21lug"; dir.create(bk, showWarnings = FALSE)
SOGLIA <- 10
rifit_uno <- function(dgp, cella, r, pow = NA) {
    percorso <- sprintf("output/fit/%s/%s/rep_%03d.rds", dgp, cella, r)
    file.copy(percorso, file.path(bk, sprintf("%s_%s_rep_%03d.rds", dgp, cella, r)))
    # niente unlink qui: saveRDS sovrascrive; cancellare prima del rifit
    # lascia un buco su disco se il run viene interrotto (successo il 21/7)
    out <- if (is.na(pow)) fitta_replica(dgp, r, cella, adapt_delta = 0.99, seed_stan = 1499)
           else fitta_replica_pow(dgp, r, sub("_c[0-9]+$", "", cella), pow,
                                  adapt_delta = 0.99, seed_stan = 7500)
    out$diagn$rifit <- TRUE
    saveRDS(out, percorso)
    cat(sprintf("%s/%s rep %03d | rhat %.3f | div %d\n",
                dgp, cella, r, out$diagn$rhat_max, out$diagn$divergenze))
    out$diagn$divergenze
}
tot_prima <- 0; tot_dopo <- 0; n_rifit <- 0; persistenti <- character()
for (spec in list(list("riferimento","lkj",NA), list("riferimento","strutturata",NA),
                  list("persistente","lkj",NA), list("persistente","strutturata",NA),
                  list("persistente","strutturata_c80",0.8))) {
    dgp <- spec[[1]]; cella <- spec[[2]]; pow <- spec[[3]]
    files <- list.files(file.path("output/fit", dgp, cella), pattern="^rep_", full.names=TRUE)
    for (fp in files) {
        f <- readRDS(fp)
        if (f$diagn$divergenze > SOGLIA) {
            if (isTRUE(f$diagn$rifit)) {
                persistenti <- c(persistenti,
                    sprintf("%s/%s rep %03d (div %d)", dgp, cella, f$r, f$diagn$divergenze))
                next
            }
            tot_prima <- tot_prima + f$diagn$divergenze; n_rifit <- n_rifit + 1
            tot_dopo <- tot_dopo + rifit_uno(dgp, cella, f$r, pow)
        }
    }
}
cat(sprintf("RIFIT DIV COMPLETATI: %d repliche | div totali prima %d -> dopo %d\n",
            n_rifit, tot_prima, tot_dopo))
if (length(persistenti)) {
    cat("DIVERGENZE PERSISTENTI DOPO RIFIT (da riportare, non rifittare oltre):\n")
    cat(paste(" -", persistenti, collapse = "\n"), "\n")
} else cat("Nessuna divergenza persistente dopo rifit.\n")
