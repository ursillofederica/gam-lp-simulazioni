# Driver: analizza tutte le celle dello studio (E1-E3, celle base e a c
# fisso) con 03_analizza.R e salva il riepilogo unico usato da figure
# e tabelle. Eseguire dalla cartella degli script.
source("03_analizza.R")

celle <- list(
    c("riferimento", "lkj"),
    c("riferimento", "strutturata"),
    c("persistente", "lkj"),
    c("persistente", "strutturata"),
    c("persistente", "lkj_c80"),
    c("persistente", "lkj_c90"),
    c("persistente", "strutturata_c80"),
    c("persistente", "lkj_c70"),
    c("persistente", "strutturata_c70"),
    c("nonlineare", "lkj"),
    c("nonlineare", "strutturata"),
    c("nonlineare", "lkj_c80"),
    c("nonlineare", "lkj_c70"),
    c("nonlineare", "strutturata_c70")
)

risultati <- list()
for (cella in celle) {
    dgp <- cella[1]; mod <- cella[2]
    cat("\n>>> analizzo", dgp, "x", mod, "...\n")
    res <- tryCatch(analizza_cella(dgp, mod),
                    error = function(e) { cat("ERRORE:", conditionMessage(e), "\n"); NULL })
    if (!is.null(res)) {
        risultati[[paste(dgp, mod, sep = "_")]] <- res
        stampa_cella(res, dgp, mod)
    }
}

saveRDS(risultati, "output/risultati_completi_22lug.rds")
cat("\nSalvato output/risultati_completi_22lug.rds\n")
