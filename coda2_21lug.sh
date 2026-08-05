#!/bin/zsh
# Secondo vagone: parte quando la coda notturna ha finito.
cd /Users/federicaursillo/Tesi/01_codice/GAM-LP/cap5
while ! grep -q 'CODA NOTTURNA COMPLETATA' output/coda_notturna_21lug.log 2>/dev/null; do sleep 300; done
Rscript nonlineare_c_fisso.R >> output/nonlineare_c_fisso_21lug.log 2>&1
