# Verifica di riproducibilita' del reperto "superficie spenta su GHKP" (23/7).
# MOTIVO: oggi il fit tanh su disco (rhat 1,225, 95 divergenze) dava IRF
# ESATTAMENTE a zero, e rifittandolo con il codice corrente e' diventato
# pulito (rhat 1,003, 0 divergenze) RECUPERANDO la non linearita' (17/17
# dentro banda). Cioe': la firma "IRF a zero + rhat ~1,2 + molte divergenze"
# NON e' necessariamente il modello che rifiuta di riportare la non
# linearita'; puo' essere un fit rotto.
# Il fit superficie/GHKP ha esattamente quella firma (rhat 1,238, 49 div).
# Il diario lo giustifica come artefatto del rank-normalized su quantita'
# degeneri, ed e' un argomento plausibile, ma dopo il caso tanh va VERIFICATO
# e non assunto: qui si rifa' con un seed diverso e si confronta.
# Se torna a zero -> il reperto e' confermato e ora e' riproducibile.
# Se recupera la soglia -> il § 5.4 va riscritto.
# eseguire dalla cartella degli script
source("02_fitta.R")

bk <- "output/fit/_backup_superficie_23lug"; dir.create(bk, showWarnings = FALSE)
fp <- "output/fit/nonlineare/superficie_strutturata/rep_001.rds"
file.copy(fp, file.path(bk, "nonlineare_superficie_strutturata_rep_001.rds"),
          overwrite = TRUE)

vecchio <- readRDS(fp)
med <- function(dr) apply(dr, 2, median)
q90 <- function(dr) apply(dr, 2, quantile, probs = c(0.05, 0.95))
cat("PRIMA: rhat", round(vecchio$diagn$rhat_max, 3),
    "| div", vecchio$diagn$divergenze,
    "| max |mediana| pos", signif(max(abs(med(vecchio$irf_draws$pos))), 3),
    "| neg", signif(max(abs(med(vecchio$irf_draws$neg))), 3), "\n")

fitta_replica("nonlineare", 1, "superficie_strutturata",
              adapt_delta = 0.95, seed_stan = 1499)
out <- readRDS(fp)

mp <- med(out$irf_draws$pos); mn <- med(out$irf_draws$neg)
bp <- q90(out$irf_draws$pos); bn <- q90(out$irf_draws$neg)
vp <- out$irf_vera$pos; vn <- out$irf_vera$neg
nc <- attr(out$irf_draws$pos, "nchains"); per <- nrow(out$irf_draws$pos) / nc
per_cat <- sapply(1:nc, function(i)
  max(abs(apply(out$irf_draws$neg[((i - 1) * per + 1):(i * per), , drop = FALSE],
                2, median))))

cat("\nDOPO: rhat", round(out$diagn$rhat_max, 3), "| div", out$diagn$divergenze,
    "| ess", round(out$diagn$ess_bulk_min), "\n")
cat("max |mediana| pos:", signif(max(abs(mp)), 3),
    "| neg:", signif(max(abs(mn)), 3), "\n")
cat("max |mediana| per catena (neg):",
    paste(signif(per_cat, 3), collapse = " "), "\n")
cat("mediana neg h0-4:", paste(sprintf("%.3f", mn[1:5]), collapse = " "), "\n")
cat("vera    neg h0-4:", paste(sprintf("%.3f", vn[1:5]), collapse = " "), "\n")
cat("mediana pos h0-4:", paste(sprintf("%.3f", mp[1:5]), collapse = " "), "\n")
cat("vera    pos h0-4:", paste(sprintf("%.3f", vp[1:5]), collapse = " "), "\n")
cat("dentro banda 90%: pos", sum(vp >= bp[1, ] & vp <= bp[2, ]), "/17",
    "| neg", sum(vn >= bn[1, ] & vn <= bn[2, ]), "/17\n")
cat("ampiezza mediana banda: pos", round(median(bp[2, ] - bp[1, ]), 4),
    "| neg", round(median(bn[2, ] - bn[1, ]), 4), "\n")
cat("\nESITO:",
    if (max(abs(mn)) < 0.01) "SUPERFICIE ANCORA SPENTA -> reperto confermato"
    else "SUPERFICIE NON PIU' SPENTA -> il § 5.4 va rivisto", "\n")
cat("VERIFICA SUPERFICIE 23/7 COMPLETATA\n")
