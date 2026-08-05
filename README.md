# GAM-LP: codice dello studio di simulazione

Codice della parte di simulazione della tesi sulle local projections bayesiane
con learning rate (GAM-LP). Il repository permette di riprodurre le analisi e le
figure del capitolo dei risultati a partire dai fit gia' calcolati (distribuiti
come asset della release, v. sotto), di variarle, o di rigenerare tutto da zero.

Nota sulla numerazione: la cartella storica del progetto si chiama `cap5`, ma
nel PDF della tesi il capitolo dei risultati e' il capitolo 4 (l'introduzione
non e' numerata).

## Struttura

Pipeline numerata (da eseguire con working directory nella radice del repo):

| Script | Cosa fa |
|---|---|
| `00_calibrazione_riferimento.R` | Calibrazione AR(5) su GDPC1 (FRED). Gli output congelati (`output/GDPC1_raw.csv`, `output/phi_calibrati.rds`) sono inclusi nel repo: non serve rieseguirla, ed e' anzi sconsigliato (FRED revisiona le serie). |
| `01_genera.R` | Generatori dei dataset simulati, seed deterministici (`SEED_BASE[dgp] + r`). |
| `02_fitta.R` | Fit Stan via cmdstanr (`fitta_replica`, `fitta_replica_pow`). Attenzione: i nomi dei modelli `_pow` sono costruiti con `paste0`, un grep semplice sui `.stan` non li trova. |
| `03_analizza.R` | Coperture e ripesatura per importance sampling: base dedotta dal `c_pow` salvato nei fit, diagnostica Pareto-k via `loo::psis`, criterio k <= 0,7. |
| `04_driver_analisi.R` | Driver dell'analisi completa; produce `output/risultati_completi_22lug.rds`. |
| `05`-`07` | Figure, certificazione diagnostica e tabelle di E1. |
| `08` | Certificazione diagnostica E2/E3. |
| `09_gam_controllo.R` | GAM frequentista di controllo per E3 (mgcv, REML). |
| `10`-`12` | Figure e tabelle E2/E3, reperti E3. |
| `13`-`15` | Reperto tanh: contrasti e superficie 3D. |

Script operativi dei giri di fit e rifit, inclusi come documentazione del
protocollo dichiarato in tesi (init = 0,1; divergenze > 10 -> rifit con
adapt_delta 0,99; rhat > 1,01 -> rifit): `rifit_div.R`, `rifit_div2.R`,
`rifit_rhat_23lug.R`, `rifit_superficie_23lug.R`, `rifit_tanh_23lug.R`,
`nonlineare_c_fisso.R`, `run_nonlineare_51_100.R`, `coda_notturna_21lug.{R,sh}`,
`coda2_21lug.sh`.

`stan/` contiene i sei modelli: lineare x {LKJ, strutturata} x {base, power
posterior} piu' le due superfici.

## Requisiti

R 4.5.1 con: cmdstanr 0.9.0 (e CmdStan 2.38.0 installato), posterior 1.6.1,
loo 2.8.0, mgcv 1.9.3, ggplot2 3.5.2, gridExtra 2.3, Hmisc 5.2.4, splines
(base R). Le versioni indicate sono quelle usate per i risultati in tesi.

## Riprodurre le analisi senza rifittare (consigliato)

I fit delle repliche (~1,5 GB, gia' in forma compatta) sono negli asset della
release `fit-v1`:

```sh
./scarica_fit.sh          # richiede GitHub CLI autenticata
```

oppure scaricare a mano i `.tar` dalla pagina della release e scompattarli in
`output/` (`tar xf <file>.tar -C output`). Ogni `rep_<r>.rds` contiene:

- `irf_draws`: 8000 draw x 17 orizzonti della IRF a posteriori;
- `log_g`: log-verosimiglianza per draw (serve alla ripesatura in c);
- `irf_vera`: la IRF vera della replica;
- `diagn`: `rhat_max`, `divergenze`, `ess_bulk_min`, `rifit`.

Con i fit su disco, l'ordine e': `04_driver_analisi.R` (analisi completa),
`06_diagnostiche_E1.R` e `08_diagnostiche_E2E3.R` (certificazione), poi gli
script di figure e tabelle. Fuori dal progetto di tesi, figure e tabelle
vengono scritte nelle cartelle locali `figure/` e `tabelle/`.

L'asset `risultati_derivati.tar` contiene anche gli `.rds` di analisi gia'
calcolati (risultati completi, diagnostiche, reperti, GAM di controllo), utili
per un confronto diretto senza rieseguire nulla.

## Variazioni

I punti di ingresso naturali: la griglia dei learning rate e le celle lette in
`04_driver_analisi.R`; il criterio di affidabilita' della ripesatura (soglia
Pareto-k) in `03_analizza.R`; la regola di lettura cella diretta/ripesata nella
lista FONTE di `10_figure_E2E3.R` e `11_tabelle_E2E3.R` (unico posto da toccare
se cambia).

## Rigenerare tutto da zero

`01_genera.R` rigenera i dataset in modo deterministico dai seed;
`02_fitta.R` rifitta (4 catene per replica, ~1400 fit: ore di calcolo su piu'
core). Il protocollo di rifit e' negli script operativi elencati sopra.
