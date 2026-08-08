# Figure per § 5.3 (E2, alta persistenza) e § 5.4 (E3, non lineare).
# Stesso stile grafico di 05_figure_E1.R.
# Regola di lettura (protocollo § 5.1): a ogni c si usa la CELLA DIRETTA se
# esiste, altrimenti la ripesatura con Pareto-k <= 0,7.
# Output in FIG: fig_E2_*.pdf, fig_E3_*.pdf. Eseguire dalla cartella degli script.
library(ggplot2)
library(gridExtra)

res <- readRDS("output/risultati_completi_22lug.rds")
gam <- readRDS("output/gam_controllo_E3.rds")
h <- 0:16
FIG <- "../../../05_tesi/figure"    # destinazione delle figure di tesi
if (!dir.exists(FIG)) { FIG <- "figure"; dir.create(FIG, showWarnings = FALSE) }   # fuori dal progetto tesi

col_str <- "#1a1a1a"; col_lkj <- "#8c8c8c"; col_acc <- "#e67e22"
col_c <- c("1" = "#1a1a1a", "0,8" = "#7f8c8d", "0,7" = "#e67e22")

tema <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 10, color = "grey20"),
        axis.text  = element_text(color = "grey20"),
        legend.title = element_blank(),
        legend.text = element_text(size = 8.5, color = "grey20"),
        legend.key.height = unit(4, "mm"),
        strip.text = element_text(size = 10.5, color = "grey20"))

# ---- Fonte valida per ogni (esperimento, modello, c) -------------------------
# cella = nome della cella da cui leggere; c_letto = riga da leggere in quella
# cella. Cella diretta dove esiste (riga base), altrimenti ripesatura k <= 0,7.
FONTE <- list(
  E2 = list(
    strutturata = list("1"   = c("persistente_strutturata", "1"),
                       "0.9" = c("persistente_strutturata", "0.9"),
                       "0.8" = c("persistente_strutturata_c80", "0.8"),
                       "0.7" = c("persistente_strutturata_c70", "0.7")),
    lkj         = list("1"   = c("persistente_lkj", "1"),
                       "0.9" = c("persistente_lkj_c90", "0.9"),
                       "0.8" = c("persistente_lkj_c80", "0.8"),
                       "0.7" = c("persistente_lkj_c70", "0.7"))),
  E3 = list(
    strutturata = list("1"   = c("nonlineare_strutturata", "1"),
                       "0.9" = c("nonlineare_strutturata", "0.9"),
                       "0.8" = c("nonlineare_strutturata_c70", "0.8"),
                       "0.7" = c("nonlineare_strutturata_c70", "0.7")),
    lkj         = list("1"   = c("nonlineare_lkj", "1"),
                       "0.9" = c("nonlineare_lkj_c80", "0.9"),
                       "0.8" = c("nonlineare_lkj_c80", "0.8"),
                       "0.7" = c("nonlineare_lkj_c70", "0.7"))))

leggi <- function(esp, mod, cc, serie = NULL) {
  f <- FONTE[[esp]][[mod]][[cc]]
  x <- res[[f[1]]][[f[2]]]
  if (!is.null(serie)) x[[serie]] else x
}
etich <- c(strutturata = "strutturata", lkj = "LKJ")

# =============================================================================
# E2 — § 5.3
# =============================================================================

## Figura 1: profilo di copertura per orizzonte, c=1 contro il c ottimo di
## ciascun modello (strutturata 0,8, LKJ 0,7). I due c* hanno COLORI DISTINTI
## in legenda (correzione richiesta: prima erano dello stesso colore).
dprof <- rbind(
  data.frame(h = h, copertura = leggi("E2", "strutturata", "1")$copertura_h,
             serie = "c = 1", modello = "strutturata"),
  data.frame(h = h, copertura = leggi("E2", "strutturata", "0.8")$copertura_h,
             serie = "c* = 0,8", modello = "strutturata"),
  data.frame(h = h, copertura = leggi("E2", "lkj", "1")$copertura_h,
             serie = "c = 1", modello = "LKJ"),
  data.frame(h = h, copertura = leggi("E2", "lkj", "0.7")$copertura_h,
             serie = "c* = 0,7", modello = "LKJ"))
dprof$modello <- factor(dprof$modello, levels = c("strutturata", "LKJ"))
dprof$serie <- factor(dprof$serie, levels = c("c = 1", "c* = 0,8", "c* = 0,7"))

p_e2_h <- ggplot(dprof, aes(h, copertura, color = serie)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_wrap(~ modello) +
  scale_color_manual(values = c("c = 1" = "#1a1a1a", "c* = 0,8" = "#2980b9",
                                "c* = 0,7" = "#e67e22"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  scale_y_continuous(limits = c(0.65, 1), breaks = seq(0.7, 1, 0.1)) +
  labs(x = "orizzonte h", y = "copertura empirica") +
  tema + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_E2_orizzonte.pdf"), p_e2_h, width = 8, height = 3.6)
cat("fig_E2_orizzonte scritta\n")

## Figura ampiezza / distorsione / RMSE per orizzonte a c=1 (gemella E1)
diag_amp <- readRDS("output/diagnostiche_E2E3.rds")$E2
dAR <- do.call(rbind, lapply(c("strutturata", "lkj"), function(m) {
  x <- diag_amp[[m]]$calibrazione$`1`
  rbind(data.frame(h = h, val = x$amp_h,  serie = "ampiezza",    modello = etich[m]),
        data.frame(h = h, val = x$bias_h, serie = "distorsione", modello = etich[m]),
        data.frame(h = h, val = x$rmse_h, serie = "RMSE",        modello = etich[m]))
}))
dAR$modello <- factor(dAR$modello, levels = c("strutturata", "LKJ"))
dAR$serie <- factor(dAR$serie, levels = c("ampiezza", "distorsione", "RMSE"))
p_ar2 <- ggplot(dAR, aes(h, val, color = modello)) +
  geom_hline(data = data.frame(serie = factor("distorsione", levels = levels(dAR$serie))),
             aes(yintercept = 0), color = "grey60", linewidth = 0.3, inherit.aes = FALSE) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  facet_wrap(~ serie, scales = "free_y") +
  scale_color_manual(values = c(strutturata = col_str, LKJ = col_lkj)) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(x = "orizzonte h", y = NULL) +
  tema + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_E2_amp_rmse.pdf"), p_ar2, width = 8.6, height = 3.4)
cat("fig_E2_amp_rmse scritta\n")

## Figura 2: curva di copertura media in c (il punto c*)
dcur <- do.call(rbind, lapply(c("strutturata", "lkj"), function(m)
  do.call(rbind, lapply(c("1", "0.9", "0.8", "0.7"), function(cc) {
    x <- leggi("E2", m, cc)
    data.frame(c = as.numeric(cc), copertura = mean(x$copertura_h),
               diretta = x$base_c == as.numeric(cc), modello = etich[m])
  }))))
dcur$modello <- factor(dcur$modello, levels = c("strutturata", "LKJ"))

p_e2_c <- ggplot(dcur, aes(c, copertura, color = modello)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.71, y = 0.893, label = "nominale 90%",
           color = "grey40", size = 2.7, hjust = 1, vjust = 1) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(shape = diretta), size = 2.2, fill = "white", stroke = 0.7) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_color_manual(values = c(strutturata = col_str, LKJ = col_lkj)) +
  scale_x_reverse(breaks = c(1, 0.9, 0.8, 0.7),
                  labels = c("1", "0,9", "0,8", "0,7")) +
  scale_y_continuous(limits = c(0.75, 0.94), breaks = seq(0.76, 0.94, 0.04),
                     labels = function(x) sub("\\.", ",", sprintf("%.2f", x))) +
  labs(x = "learning rate c", y = "copertura media sugli orizzonti") +
  tema + theme(legend.position = c(0.18, 0.86),
               legend.background = element_rect(fill = "white", color = NA))
ggsave(file.path(FIG, "fig_E2_curva_c.pdf"), p_e2_c, width = 5.4, height = 3.6)
cat("fig_E2_curva_c scritta\n")

## Figura divergenze per replica, sette celle (gemella di fig_E1_divergenze)
diag_e2e3 <- readRDS("output/diagnostiche_E2E3.rds")
etich_cell <- c(strutturata = "strutturata, c=1", strutturata_c80 = "strutturata, c=0,8",
                strutturata_c70 = "strutturata, c=0,7", lkj = "LKJ, c=1",
                lkj_c90 = "LKJ, c=0,9", lkj_c80 = "LKJ, c=0,8", lkj_c70 = "LKJ, c=0,7")
ordine <- c("strutturata", "strutturata_c80", "strutturata_c70",
            "lkj", "lkj_c90", "lkj_c80", "lkj_c70")
dv2 <- do.call(rbind, lapply(ordine, function(m)
  data.frame(cella = etich_cell[m], div = diag_e2e3$E2[[m]]$repliche$div)))
dv2$cella <- factor(dv2$cella, levels = etich_cell[ordine])
p_div2 <- ggplot(dv2, aes(div)) +
  geom_bar(fill = "#8c8c8c", width = 0.7) +
  facet_wrap(~ cella, ncol = 4) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "transizioni divergenti", y = "numero di repliche") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        strip.text = element_text(size = 10, color = "grey20"))
ggsave(file.path(FIG, "fig_E2_divergenze.pdf"), p_div2, width = 8.6, height = 4.2)
cat("fig_E2_divergenze scritta\n")

# =============================================================================
# E3 — § 5.4
# =============================================================================

## Figura 3: copertura per orizzonte nei due contrasti, c = 1
dE3 <- do.call(rbind, lapply(c("strutturata", "lkj"), function(m)
  do.call(rbind, lapply(c("pos", "neg"), function(s)
    data.frame(h = h, copertura = leggi("E3", m, "1", s)$copertura_h,
               contrasto = ifelse(s == "pos", "delta = +2", "delta = -2"),
               modello = etich[m])))))
dE3$modello <- factor(dE3$modello, levels = c("strutturata", "LKJ"))
dE3$contrasto <- factor(dE3$contrasto, levels = c("delta = +2", "delta = -2"))

p_e3_h <- ggplot(dE3, aes(h, copertura, color = contrasto)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_wrap(~ modello) +
  scale_color_manual(values = c("delta = +2" = col_str, "delta = -2" = col_acc),
                     labels = c(expression(delta == +2), expression(delta == -2))) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "orizzonte h", y = "copertura empirica") +
  tema + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_E3_contrasti.pdf"), p_e3_h, width = 8, height = 3.6)
cat("fig_E3_contrasti scritta\n")

## Figura 4: temperare non muove il contrasto negativo (strutturata)
dtemp <- do.call(rbind, lapply(c("pos", "neg"), function(s)
  do.call(rbind, lapply(c("1", "0.7"), function(cc)
    data.frame(h = h, copertura = leggi("E3", "strutturata", cc, s)$copertura_h,
               c = sub("\\.", ",", cc),
               contrasto = ifelse(s == "pos", "delta == +2", "delta == -2"))))))
dtemp$contrasto <- factor(dtemp$contrasto, levels = c("delta == +2", "delta == -2"))
dtemp$c <- factor(dtemp$c, levels = c("1", "0,7"))

p_e3_temp <- ggplot(dtemp, aes(h, copertura, color = c)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_wrap(~ contrasto, labeller = label_parsed) +
  scale_color_manual(values = c("1" = "#1a1a1a", "0,7" = col_acc), name = "c") +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "orizzonte h", y = "copertura empirica") +
  tema + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_E3_temperare.pdf"), p_e3_temp, width = 8, height = 3.6)
cat("fig_E3_temperare scritta\n")

## Figura 5: GAM di controllo — la superficie spenta anche in frequentista
dgam <- rbind(
  data.frame(h = gam$tab$h, irf = gam$tab$vera_pos, serie = "vera",
             contrasto = "delta == +2"),
  data.frame(h = gam$tab$h, irf = gam$tab$gam_pos, serie = "GAM (REML)",
             contrasto = "delta == +2"),
  data.frame(h = gam$tab$h, irf = gam$tab$vera_neg, serie = "vera",
             contrasto = "delta == -2"),
  data.frame(h = gam$tab$h, irf = gam$tab$gam_neg, serie = "GAM (REML)",
             contrasto = "delta == -2"))
dgam$contrasto <- factor(dgam$contrasto, levels = c("delta == +2", "delta == -2"))
dgam$serie <- factor(dgam$serie, levels = c("vera", "GAM (REML)"))

p_gam <- ggplot(dgam, aes(h, irf, color = serie, linetype = serie)) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.1) +
  facet_wrap(~ contrasto, labeller = label_parsed) +
  scale_color_manual(values = c(vera = col_str, `GAM (REML)` = col_acc)) +
  scale_linetype_manual(values = c(vera = "solid", `GAM (REML)` = "22")) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(x = "orizzonte h", y = "risposta d'impulso") +
  tema + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_E3_gam_controllo.pdf"), p_gam, width = 8, height = 3.4)
cat("fig_E3_gam_controllo scritta\n")

## Figura divergenze per replica, celle E3 (gemella di fig_E2_divergenze)
etich_e3 <- c(strutturata = "strutturata, c=1", strutturata_c70 = "strutturata, c=0,7",
              lkj = "LKJ, c=1", lkj_c80 = "LKJ, c=0,8", lkj_c70 = "LKJ, c=0,7")
ordine3 <- c("strutturata", "strutturata_c70", "lkj", "lkj_c80", "lkj_c70")
dv3 <- do.call(rbind, lapply(ordine3, function(m)
  data.frame(cella = etich_e3[m], div = diag_e2e3$E3[[m]]$repliche$div)))
dv3$cella <- factor(dv3$cella, levels = etich_e3[ordine3])
p_div3 <- ggplot(dv3, aes(div)) +
  geom_bar(fill = "#8c8c8c", width = 0.7) +
  facet_wrap(~ cella, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 10, 2)) +
  labs(x = "transizioni divergenti", y = "numero di repliche") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        strip.text = element_text(size = 10, color = "grey20"))
ggsave(file.path(FIG, "fig_E3_divergenze.pdf"), p_div3, width = 8.6, height = 4.6)
cat("fig_E3_divergenze scritta\n")

## Figura copertura media vs learning rate, per contrasto e modello (E3)
dE3c <- do.call(rbind, lapply(c("strutturata", "lkj"), function(m)
  do.call(rbind, lapply(c("1", "0.9", "0.8", "0.7"), function(cc) rbind(
    data.frame(c = as.numeric(cc), copertura = mean(leggi("E3", m, cc, "pos")$copertura_h),
               contrasto = "delta = +2", modello = etich[m]),
    data.frame(c = as.numeric(cc), copertura = mean(leggi("E3", m, cc, "neg")$copertura_h),
               contrasto = "delta = -2", modello = etich[m]))))))
dE3c$modello <- factor(dE3c$modello, levels = c("strutturata", "LKJ"))
dE3c$contrasto <- factor(dE3c$contrasto, levels = c("delta = +2", "delta = -2"))
p_e3_curva <- ggplot(dE3c, aes(c, copertura, color = contrasto)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_wrap(~ modello) +
  scale_color_manual(values = c("delta = +2" = col_str, "delta = -2" = col_acc),
                     breaks = c("delta = +2", "delta = -2"),
                     labels = c(expression(delta == +2), expression(delta == -2))) +
  scale_x_reverse(breaks = c(1, 0.9, 0.8, 0.7), labels = c("1","0,9","0,8","0,7")) +
  labs(x = "learning rate c", y = "copertura media sugli orizzonti") +
  tema + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig_E3_curva_c.pdf"), p_e3_curva, width = 8, height = 3.8)
cat("fig_E3_curva_c scritta\n")
