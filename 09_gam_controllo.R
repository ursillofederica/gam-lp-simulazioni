# GAM frequentista di controllo per E3/GHKP: stessa replica del reperto
# bayesiano (rep_001), superficie tensoriale gemella (bs, k = c(6,8)),
# REML, controlli lineari per orizzonte. Domanda: anche la macchina
# frequentista spegne la non linearita'? Eseguire dalla cartella degli script.
library(mgcv)
source("01_genera.R")   # TT, H, H1

d <- readRDS("output/dati/nonlineare/rep_001.rds")
colnames(d$X_ctrl) <- paste0("ctrl", seq_len(ncol(d$X_ctrl)))
lungo <- do.call(rbind, lapply(0:H, function(hh) {
  data.frame(y = d$y_mat[, hh + 1], s = d$s, hh = hh, d$X_ctrl)
}))
lungo$hf <- factor(lungo$hh)
ctrl <- paste0("hf:", colnames(d$X_ctrl), collapse = " + ")
frm <- as.formula(paste("y ~ te(s, hh, bs = 'bs', k = c(6, 8)) +", ctrl))
fit <- gam(frm, data = lungo, method = "REML")

nd <- function(sv) data.frame(s = sv, hh = 0:H, hf = factor(0, levels = levels(lungo$hf)),
                              matrix(0, H1, ncol(d$X_ctrl),
                                     dimnames = list(NULL, colnames(d$X_ctrl))))
f  <- function(sv) predict(fit, nd(sv))
irf_pos <- as.numeric(f(2)  - f(0))
irf_neg <- as.numeric(f(-2) - f(0))

out <- list(fit_summary = summary(fit), sp = fit$sp,
            edf_te = sum(summary(fit)$edf),
            tab = data.frame(h = 0:H,
                             gam_pos = round(irf_pos, 3), vera_pos = round(d$irf_vera$pos, 3),
                             gam_neg = round(irf_neg, 3), vera_neg = round(d$irf_vera$neg, 3)))
saveRDS(out, "output/gam_controllo_E3.rds")
cat("EDF della superficie te(s,h):", round(out$edf_te, 2), "\n")
cat("parametri di lisciatura (REML):", format(fit$sp, digits = 3), "\n\n")
print(out$tab, row.names = FALSE)
cat("\nmax |IRF gam| pos:", round(max(abs(irf_pos)), 3),
    "| neg:", round(max(abs(irf_neg)), 3), "\n")
