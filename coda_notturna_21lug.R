# Coda notturna 21/7 (decisione di Federica, sera): dopo P3,
# 1) superficie su GHKP, 1 replica, da archiviare (diagnosi § 5.4;
#    il GAM frequentista di controllo si fa domani, scelte di modello
#    da discutere);
# 2) celle c-fisso: completamenti a R=100 delle tre celle esistenti
#    + celle nuove c70 per entrambi i modelli (protocollo k<=0.7:
#    a c=0.7 il Pareto-k supera soglia anche per la strutturata).
# Tutto riprendibile: fitta_blocco* saltano le repliche gia' su disco.
# eseguire dalla cartella degli script
source("02_fitta.R")

cat("== 1. Superficie su GHKP (rep 001, archivio diagnosi) ==\n")
fitta_blocco("nonlineare", "superficie_strutturata", 1)

cat("== 2a. Completamenti a R=100 delle celle esistenti ==\n")
fitta_blocco_pow("persistente", "lkj",         51:100, 0.8)
fitta_blocco_pow("persistente", "lkj",         51:100, 0.9)
fitta_blocco_pow("persistente", "strutturata", 51:100, 0.8)

cat("== 2b. Celle nuove c=0.7, R=100 ==\n")
fitta_blocco_pow("persistente", "lkj",         1:100, 0.7)
fitta_blocco_pow("persistente", "strutturata", 1:100, 0.7)

cat("CODA NOTTURNA COMPLETATA\n")
