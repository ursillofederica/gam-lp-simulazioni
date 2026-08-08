# Celle a c fisso per E3/lkj (approvate da Federica, 21/7 sera):
# la ripesatura lkj su nonlineare regge solo a c=0.9 (k: 0.50 a 0.9,
# 0.98 a 0.8, 1.49 a 0.7) -> fit diretti a 0.8 e 0.7, R=100.
# eseguire dalla cartella degli script
source("02_fitta.R")
fitta_blocco_pow("nonlineare", "lkj", 1:100, 0.8)
fitta_blocco_pow("nonlineare", "lkj", 1:100, 0.7)
# strutturata a 0.7: k=0.98 anche per lei (approvato 21/7 sera)
fitta_blocco_pow("nonlineare", "strutturata", 1:100, 0.7)
cat("CELLE NONLINEARE C-FISSO COMPLETATE\n")
