# Reperti di E3 (fuori dalle celle MC: una replica per disegno).
# Tre pezzi che mancavano al rapporto e alla tesi:
#   1. superficie bayesiana su GHKP  (output/fit/nonlineare/superficie_strutturata)
#   2. GAM frequentista di controllo (output/gam_controllo_E3.rds)
#   3. recupero a segnale forte, DGP tanh (output/fit/tanh/superficie_strutturata)
# Produce: output/reperti_E3.rds (riassunti leggeri, cosi' il rapporto non
# rilegge i fit da 2 MB), le figure di tesi e tab_e3_reperti.tex.
# Eseguire dalla cartella degli script.
library(ggplot2)
FIG <- "../../../05_tesi/figure"    # destinazione delle figure di tesi
if (!dir.exists(FIG)) { FIG <- "figure"; dir.create(FIG, showWarnings = FALSE) }   # fuori dal progetto tesi
TAB <- "../../../05_tesi/tabelle"   # destinazione delle tabelle di tesi
if (!dir.exists(TAB)) { TAB <- "tabelle"; dir.create(TAB, showWarnings = FALSE) }   # fuori dal progetto tesi
h <- 0:16
col_str <- "#1a1a1a"; col_acc <- "#e67e22"
tema <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        axis.title = element_text(size = 10, color = "grey20"),
        axis.text = element_text(color = "grey20"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8.5, color = "grey20"),
        legend.position = "bottom",
        strip.text = element_text(size = 10.5, color = "grey20"))

# ---- riassunto di una serie di draw: mediana, banda 90%, mediane per catena --
riassumi <- function(dr) {
  q <- apply(dr, 2, quantile, probs = c(0.05, 0.5, 0.95))
  nc <- attr(dr, "nchains"); per <- nrow(dr) / nc
  per_catena <- sapply(1:nc, function(i)
    apply(dr[((i - 1) * per + 1):(i * per), , drop = FALSE], 2, median))
  list(lo = q[1, ], med = q[2, ], hi = q[3, ], per_catena = per_catena)
}

# =============================================================================
# 1. Superficie bayesiana su GHKP
# =============================================================================
sup <- readRDS("output/fit/nonlineare/superficie_strutturata/rep_001.rds")
S_pos <- riassumi(sup$irf_draws$pos); S_neg <- riassumi(sup$irf_draws$neg)

cat("== superficie su GHKP, rep_001\n")
cat("   rhat", round(sup$diagn$rhat_max, 3), "| div", sup$diagn$divergenze,
    "| ess", round(sup$diagn$ess_bulk_min), "\n")
cat("   max |mediana| pos:", signif(max(abs(S_pos$med)), 3),
    "| neg:", signif(max(abs(S_neg$med)), 3), "\n")
cat("   max |mediana per catena| (unanimita'): pos",
    signif(max(abs(S_pos$per_catena)), 3), "| neg",
    signif(max(abs(S_neg$per_catena)), 3), "\n")
cat("   ampiezza mediana banda 90%: pos", round(median(S_pos$hi - S_pos$lo), 3),
    "| neg", round(median(S_neg$hi - S_neg$lo), 3), "\n")
cat("   il vero entra nella banda? pos",
    sum(sup$irf_vera$pos >= S_pos$lo & sup$irf_vera$pos <= S_pos$hi), "/17 | neg",
    sum(sup$irf_vera$neg >= S_neg$lo & sup$irf_vera$neg <= S_neg$hi), "/17\n")

d_sup <- rbind(
  data.frame(h = h, med = S_pos$med, lo = S_pos$lo, hi = S_pos$hi,
             vera = sup$irf_vera$pos, contrasto = "delta == +2"),
  data.frame(h = h, med = S_neg$med, lo = S_neg$lo, hi = S_neg$hi,
             vera = sup$irf_vera$neg, contrasto = "delta == -2"))
d_sup$contrasto <- factor(d_sup$contrasto, levels = c("delta == +2", "delta == -2"))

p_sup <- ggplot(d_sup, aes(h)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "banda al 90% (superficie)"),
              alpha = 0.35) +
  geom_line(aes(y = med, color = "mediana a posteriori"), linewidth = 0.8) +
  geom_line(aes(y = vera, color = "vera"), linewidth = 0.8) +
  geom_point(aes(y = vera, color = "vera"), size = 1.1) +
  facet_wrap(~ contrasto, labeller = label_parsed) +
  scale_fill_manual(values = c("banda al 90% (superficie)" = "#bdc3c7")) +
  scale_color_manual(values = c("mediana a posteriori" = col_acc, "vera" = col_str)) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(x = "orizzonte h", y = "risposta d'impulso") + tema
ggsave(file.path(FIG, "fig_E3_superficie.pdf"), p_sup, width = 8, height = 3.6)
cat("   fig_E3_superficie scritta\n")

# =============================================================================
# 2. GAM frequentista di controllo (gia' stimato: 09_gam_controllo.R)
# =============================================================================
gam <- readRDS("output/gam_controllo_E3.rds")
cat("\n== GAM di controllo: EDF", round(gam$edf_te, 2),
    "| sp", paste(format(gam$sp, digits = 4), collapse = " / "),
    "| max |IRF| pos", max(abs(gam$tab$gam_pos)),
    "neg", max(abs(gam$tab$gam_neg)), "\n")

# =============================================================================
# 3. Recupero a segnale forte: DGP tanh
# =============================================================================
# deltas del tanh = c(2, 1) (taglia, non segno): irf_draws$pos e' delta = 2,
# $neg e' delta = 1. Il vero e' salvato come d1 (delta 1) e d2 (delta 2):
# l'accoppiamento va verificato, non assunto.
th <- readRDS("output/fit/tanh/superficie_strutturata/rep_001.rds")
T_d2 <- riassumi(th$irf_draws$pos)   # delta = 2
T_d1 <- riassumi(th$irf_draws$neg)   # delta = 1
# Controllo dell'accoppiamento: ogni serie deve stare piu' vicina al proprio
# vero che all'altro. Se il fit e' degenere (tutte le mediane a zero) le due
# distanze coincidono e il controllo non discrimina: va segnalato, non ignorato.
d_22 <- max(abs(T_d2$med - th$irf_vera$d2)); d_21 <- max(abs(T_d2$med - th$irf_vera$d1))
d_11 <- max(abs(T_d1$med - th$irf_vera$d1)); d_12 <- max(abs(T_d1$med - th$irf_vera$d2))
if (isTRUE(all.equal(d_22, d_21)) || isTRUE(all.equal(d_11, d_12)))
  stop("fit tanh degenere (mediane a zero): l'accoppiamento dei contrasti non e' ",
       "verificabile e il reperto non e' utilizzabile. Rifittare.")
stopifnot(d_22 < d_21, d_11 < d_12)

dentro <- function(R, v) sum(v >= R$lo & v <= R$hi)
cat("\n== tanh, rep_001 | rhat", round(th$diagn$rhat_max, 3),
    "| div", th$diagn$divergenze, "| ess", round(th$diagn$ess_bulk_min), "\n")
cat("   delta=1: dentro banda", dentro(T_d1, th$irf_vera$d1), "/17 | delta=2:",
    dentro(T_d2, th$irf_vera$d2), "/17\n")
cat("   max |errore mediana|: delta=1",
    round(max(abs(T_d1$med - th$irf_vera$d1)), 3), "| delta=2",
    round(max(abs(T_d2$med - th$irf_vera$d2)), 3), "\n")

d_th <- rbind(
  data.frame(h = h, med = T_d1$med, lo = T_d1$lo, hi = T_d1$hi,
             vera = th$irf_vera$d1, contrasto = "delta == 1"),
  data.frame(h = h, med = T_d2$med, lo = T_d2$lo, hi = T_d2$hi,
             vera = th$irf_vera$d2, contrasto = "delta == 2"))
d_th$contrasto <- factor(d_th$contrasto, levels = c("delta == 1", "delta == 2"))

p_th <- ggplot(d_th, aes(h)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = "banda al 90% (superficie)"),
              alpha = 0.35) +
  geom_line(aes(y = med, color = "mediana a posteriori"), linewidth = 0.8) +
  geom_line(aes(y = vera, color = "vera"), linewidth = 0.8) +
  geom_point(aes(y = vera, color = "vera"), size = 1.1) +
  facet_wrap(~ contrasto, labeller = label_parsed) +
  scale_fill_manual(values = c("banda al 90% (superficie)" = "#bdc3c7")) +
  scale_color_manual(values = c("mediana a posteriori" = col_acc, "vera" = col_str)) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(x = "orizzonte h", y = "risposta d'impulso") + tema
ggsave(file.path(FIG, "fig_E3_tanh.pdf"), p_th, width = 8, height = 3.6)
cat("   fig_E3_tanh scritta\n")

# =============================================================================
# Riassunti leggeri per il rapporto + tabella di tesi
# =============================================================================
rep_out <- list(
  superficie = list(diagn = sup$diagn, pos = S_pos, neg = S_neg,
                    vera = sup$irf_vera, log_g_catena = sapply(1:4, function(i)
                      mean(sup$log_g[((i - 1) * 2000 + 1):(i * 2000)]))),
  gam = gam[c("sp", "edf_te", "tab")],
  tanh = list(diagn = th$diagn, d1 = T_d1, d2 = T_d2, vera = th$irf_vera,
              dentro = c(d1 = dentro(T_d1, th$irf_vera$d1),
                         d2 = dentro(T_d2, th$irf_vera$d2))))
saveRDS(rep_out, "output/reperti_E3.rds")
cat("\nSalvato output/reperti_E3.rds\n")

it <- function(x, k = 2) gsub("\\.", "{,}", formatC(x, format = "f", digits = k))
writeLines(c(
  "\\begin{tabular}{llccc}", "\\toprule",
  "reperto & macchina di stima & $\\max|\\text{IRF}|$ stimata & $\\max|\\text{IRF}|$ vera & esito \\\\",
  "\\midrule",
  paste("superficie su GHKP & posterior, 4 catene &",
        paste0("$", it(max(abs(c(S_pos$med, S_neg$med))), 3), "$ &"),
        paste0("$", it(max(abs(unlist(sup$irf_vera)))), "$ &"),
        "superficie spenta \\\\"),
  paste("GAM di controllo & REML (\\texttt{mgcv}) &",
        paste0("$", it(max(abs(c(gam$tab$gam_pos, gam$tab$gam_neg)))), "$ &"),
        paste0("$", it(max(abs(c(gam$tab$vera_pos, gam$tab$vera_neg)))), "$ &"),
        "superficie spenta \\\\"),
  paste("superficie su tanh & posterior, 4 catene &",
        paste0("$", it(max(abs(c(T_d1$med, T_d2$med)))), "$ &"),
        paste0("$", it(max(abs(unlist(th$irf_vera)))), "$ &"),
        paste0("$", dentro(T_d1, th$irf_vera$d1) + dentro(T_d2, th$irf_vera$d2),
               "/34$ dentro banda \\\\")),
  "\\bottomrule", "\\end{tabular}",
  "\\\\[2pt]\\footnotesize Una replica per disegno (rep.\\ 001): dimostrazioni, non studi Monte Carlo. Le prime due righe sono la stessa replica stimata con due macchine diverse."),
  file.path(TAB, "tab_e3_reperti.tex"))
cat("tab_e3_reperti scritta\n")
