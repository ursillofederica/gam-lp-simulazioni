# GAM-LP: simulation study code

Code for the simulation part of a thesis on Bayesian local projections with a
learning rate (GAM-LP). The repository makes it possible to reproduce the
analyses and figures of the results chapter starting from the precomputed fits
(distributed as release assets, see below), to vary them, or to regenerate
everything from scratch.

## Structure

Numbered pipeline (to be run with the working directory at the repo root):

| Script | What it does |
|---|---|
| `00_calibrazione_riferimento.R` | AR(5) calibration on GDPC1 (FRED). The frozen outputs (`output/GDPC1_raw.csv`, `output/phi_calibrati.rds`) are included in the repo: there is no need to rerun it, and doing so is in fact discouraged (FRED revises its series). |
| `01_genera.R` | Generators of the simulated datasets, one per DGP, with deterministic seeds (`SEED_BASE[dgp] + r`). |
| `02_fitta.R` | Stan fits via cmdstanr (`fitta_replica`, `fitta_replica_pow`), plus the spline bases and model constants reused downstream. Caveat: the `_pow` model names are built with `paste0`, so a plain grep over the `.stan` files will not find them. |
| `03_analizza.R` | Coverages and importance-sampling reweighting across learning rates: base inferred from the `c_pow` stored in the fits, Pareto-k diagnostics via `loo::psis`, reliability criterion k <= 0.7. |
| `04_driver_analisi.R` | Driver for the full analysis over all experiments and cells; produces `output/risultati_completi_22lug.rds`. |
| `05_figure_E1.R` | E1 figures (ggplot2): divergences and Pareto-k per replica, coverage by horizon, interval width + RMSE. |
| `06_diagnostiche_E1.R` | E1 certification: convergence, reliability of the reweighting, extended calibration (exact MC error, RMSE, joint coverage). Produces `output/diagnostiche_E1.rds`. |
| `07_tabelle_E1.R` | E1 LaTeX tables: protocol, diagnostics, reweighting reliability, cell summary, per-horizon profile. |
| `08_diagnostiche_E2E3.R` | Same certification generalized to E2/E3: cells with a fitted base `c`, two contrasts in E3. Produces `output/diagnostiche_E2E3.rds`. |
| `09_gam_controllo.R` | Frequentist GAM control model for E3 (mgcv, REML). Produces `output/gam_controllo_E3.rds`. |
| `10_figure_E2E3.R` | E2/E3 figures, same style as `05`. At each `c`, the direct cell is used if it exists, otherwise the reweighting with Pareto-k <= 0.7 (the FONTE list, shared with `11`). |
| `11_tabelle_E2E3.R` | E2/E3 LaTeX tables (diagnostics, cross-validation, reweighting), same direct-cell/reweighted reading rule as `10`. |
| `12_reperti_E3.R` | E3 findings outside the Monte Carlo cells (one replica per design): Bayesian surface on the GHKP-style DGP, control GAM, strong-signal tanh recovery. Produces `output/reperti_E3.rds` (light summaries), figures and a table. |
| `13_figura_tanh_contrasti.R` | tanh recovery, figure 1: both contrasts (delta = 1 and 2) in one panel, posterior median, 90% band and true value. No refit, reads `reperti_E3.rds`. |
| `14_superficie_tanh_3d.R` | Extracts the fitted surface f(s,h) on a grid of s: reruns the tanh fit once (same seed as the stored finding) only to recover the spline coefficients; does not overwrite the stored fit. Produces `output/superficie_tanh_griglia.rds`. |
| `15_figura_superficie_3d.R` | tanh recovery, figure 2: 3D surface f(s,h), true vs estimated (base-R `persp`, two panels). |

Operational scripts from the fit and refit rounds, included as documentation of
the protocol declared in the thesis (init = 0.1; divergences > 10 -> refit with
adapt_delta 0.99; rhat > 1.01 -> refit): `rifit_div.R`, `rifit_div2.R`,
`rifit_rhat_23lug.R`, `rifit_superficie_23lug.R`, `rifit_tanh_23lug.R`,
`nonlineare_c_fisso.R`, `run_nonlineare_51_100.R`, `coda_notturna_21lug.{R,sh}`,
`coda2_21lug.sh`.

`stan/` contains the six models: linear x {LKJ, structured} x {base, power
posterior} plus the two surface models.

## Requirements

R 4.5.1 with: cmdstanr 0.9.0 (and CmdStan 2.38.0 installed), posterior 1.6.1,
loo 2.8.0, mgcv 1.9.3, ggplot2 3.5.2, gridExtra 2.3, Hmisc 5.2.4, splines
(base R). The versions listed are the ones used for the results in the thesis.

## Reproducing the analyses without refitting (recommended)

The replica fits (~1.5 GB, already in compact form) are in the assets of the
`fit-v1` release:

```sh
./scarica_fit.sh          # requires an authenticated GitHub CLI
```

or download the `.tar` files manually from the release page and extract them
into `output/` (`tar xf <file>.tar -C output`). Each `rep_<r>.rds` contains:

- `irf_draws`: 8000 draws x 17 horizons of the posterior IRF;
- `log_g`: per-draw log-likelihood (needed for the reweighting in c);
- `irf_vera`: the true IRF of the replica;
- `diagn`: `rhat_max`, `divergenze`, `ess_bulk_min`, `rifit`.

With the fits on disk, the order is: `04_driver_analisi.R` (full analysis),
`06_diagnostiche_E1.R` and `08_diagnostiche_E2E3.R` (certification), then the
figure and table scripts. Outside the thesis project, figures and tables are
written to the local `figure/` and `tabelle/` folders.

The `risultati_derivati.tar` asset also contains the precomputed analysis
`.rds` files (full results, diagnostics, findings, control GAM), useful for a
direct comparison without rerunning anything.

## Variations

The natural entry points: the learning-rate grid and the cells read in
`04_driver_analisi.R`; the reweighting reliability criterion (Pareto-k
threshold) in `03_analizza.R`; the direct-cell/reweighted-cell reading rule in
the FONTE list of `10_figure_E2E3.R` and `11_tabelle_E2E3.R` (the only place
to touch if it changes).

## Regenerating everything from scratch

`01_genera.R` regenerates the datasets deterministically from the seeds;
`02_fitta.R` refits (4 chains per replica, ~1400 fits: hours of computation on
several cores). The refit protocol is in the operational scripts listed above.

<!-- TODO (Federica): AI-use statement per the UniPD guidelines, e.g.
"Parts of this codebase and this README were drafted with the assistance of
an AI tool (Claude) and subsequently reviewed and reworked by the author."
Adapt the wording after reworking the README. -->
