#!/bin/zsh
# Attende la fine di P3 (o lo riprende se il run e' morto), poi lancia
# la coda notturna. Staccato dalla sessione (nohup): sopravvive ai /login.
cd /Users/federicaursillo/Tesi/01_codice/GAM-LP/cap5
LOG_P3=output/run_nonlineare_51_100_21lug.log
while true; do
    grep -q 'NONLINEARE 51:100 COMPLETATO' "$LOG_P3" && break
    # se il log e' fermo da >20 min il run P3 e' morto: riprendilo (e' riprendibile)
    if [ -n "$(find "$LOG_P3" -mmin +20 2>/dev/null)" ]; then
        echo "[coda] P3 fermo da >20 min: riprendo run_nonlineare_51_100.R" >> "$LOG_P3"
        Rscript run_nonlineare_51_100.R >> "$LOG_P3" 2>&1
    fi
    sleep 120
done
Rscript coda_notturna_21lug.R >> output/coda_notturna_21lug.log 2>&1
