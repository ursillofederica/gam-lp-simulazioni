# Rifit del reperto tanh (23/7). MOTIVO: il fit su disco
# (output/fit/tanh/superficie_strutturata/rep_001.rds, 20/7 ore 14:17) ha le
# IRF mediane ESATTAMENTE a zero, mentre il riassunto salvato
# (output/tanh_recupero_rep1.rds, stesso giorno ore 14:57) contiene bande che
# seguono la verita' (16/17 e 17/17 dentro). I due sono incompatibili: il fit
# che ha prodotto il riassunto NON e' mai stato salvato. DGP, s e irf_vera
# coincidono nei tre file, quindi non e' un disallineamento di dati.
# In piu' il fit su disco e' del 20/7, cioe' PRECEDENTE al fix init = 0,1 del
# 21/7: le IRF a zero possono essere l'ennesima catena degenere.
# Qui si rifa' il fit con il codice corrente e si confronta col riassunto.
# eseguire dalla cartella degli script
source("02_fitta.R")

bk <- "output/fit/_backup_tanh_23lug"; dir.create(bk, showWarnings = FALSE)
fp <- "output/fit/tanh/superficie_strutturata/rep_001.rds"
file.copy(fp, file.path(bk, "tanh_superficie_strutturata_rep_001.rds"), overwrite = TRUE)

vecchio <- readRDS(fp)
cat("PRIMA: rhat", round(vecchio$diagn$rhat_max, 3),
    "| div", vecchio$diagn$divergenze,
    "| max |mediana IRF|",
    signif(max(abs(apply(vecchio$irf_draws$pos, 2, median)),
               abs(apply(vecchio$irf_draws$neg, 2, median))), 3), "\n")

out <- fitta_replica("tanh", 1, "superficie_strutturata",
                     adapt_delta = 0.95, seed_stan = 1499)
out <- readRDS(fp)   # fitta_replica salva da se' (anche in caso di auto-rifit)

med <- function(dr) apply(dr, 2, median)
q90 <- function(dr) apply(dr, 2, quantile, probs = c(0.05, 0.95))
m2 <- med(out$irf_draws$pos); m1 <- med(out$irf_draws$neg)   # deltas = c(2, 1)
b2 <- q90(out$irf_draws$pos); b1 <- q90(out$irf_draws$neg)
v1 <- out$irf_vera$d1; v2 <- out$irf_vera$d2

cat("\nDOPO: rhat", round(out$diagn$rhat_max, 3),
    "| div", out$diagn$divergenze, "| ess", round(out$diagn$ess_bulk_min), "\n")
cat("mediana delta=1 h0-4:", paste(sprintf("%.3f", m1[1:5]), collapse = " "), "\n")
cat("vera    delta=1 h0-4:", paste(sprintf("%.3f", v1[1:5]), collapse = " "), "\n")
cat("mediana delta=2 h0-4:", paste(sprintf("%.3f", m2[1:5]), collapse = " "), "\n")
cat("vera    delta=2 h0-4:", paste(sprintf("%.3f", v2[1:5]), collapse = " "), "\n")
cat("dentro banda 90%: delta=1", sum(v1 >= b1[1, ] & v1 <= b1[2, ]), "/17",
    "| delta=2", sum(v2 >= b2[1, ] & v2 <= b2[2, ]), "/17\n")
cat("ampiezza mediana banda: delta=1", round(median(b1[2, ] - b1[1, ]), 3),
    "| delta=2", round(median(b2[2, ] - b2[1, ]), 3), "\n")

vecchio_riass <- readRDS("output/tanh_recupero_rep1.rds")
cat("\nconfronto col riassunto orfano del 20/7 (mediana delta=1 h0-4):",
    paste(sprintf("%.3f", vecchio_riass$b1[1:5, 2]), collapse = " "), "\n")
cat("RIFIT TANH 23/7 COMPLETATO\n")
