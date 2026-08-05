# Completa P3 (nonlineare GHKP), modello lineare, repliche 51:100
# per lkj e strutturata (necessario per il § 5.4). fitta_blocco salta
# le repliche già su disco: riprendibile per costruzione.
setwd("/Users/federicaursillo/Tesi/01_codice/GAM-LP/cap5")
source("02_fitta.R")
fitta_blocco("nonlineare", "lkj", 51:100)
fitta_blocco("nonlineare", "strutturata", 51:100)
cat("NONLINEARE 51:100 COMPLETATO\n")
