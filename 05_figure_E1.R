# Figure per § 5.2 (E1, disegno di riferimento), in ggplot2.
# Input:  output/risultati_completi_22lug.rds, output/diagnostiche_E1.rds.
# Output in FIG: fig_E1_divergenze, fig_E1_pareto_k, fig_E1_copertura,
# fig_E1_amp_rmse. Eseguire dalla cartella degli script.
library(ggplot2)
FIG <- "../../../05_tesi/figure"    # destinazione delle figure di tesi
if (!dir.exists(FIG)) { FIG <- "figure"; dir.create(FIG, showWarnings = FALSE) }   # fuori dal progetto tesi

res  <- readRDS("output/risultati_completi_22lug.rds")
h <- 0:16

col_str <- "#1a1a1a"; col_lkj <- "#8c8c8c"

tema <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        plot.title = element_text(size = 11, hjust = 0.5),
        axis.title = element_text(size = 10, color = "grey20"),
        axis.text  = element_text(color = "grey20"),
        legend.position = c(0.72, 0.14),
        legend.title = element_blank(),
        legend.text = element_text(size = 8.5, color = "grey20"),
        legend.key.height = unit(4, "mm"),
        legend.background = element_rect(fill = "white", color = NA))

# --- Pannello A: copertura per orizzonte, c = 1 -------------------------------
dA <- rbind(
  data.frame(h = h, copertura = res$riferimento_strutturata$`1`$copertura_h,
             modello = "strutturata"),
  data.frame(h = h, copertura = res$riferimento_lkj$`1`$copertura_h,
             modello = "LKJ"))
dA$modello <- factor(dA$modello, levels = c("strutturata", "LKJ"))

pA <- ggplot(dA, aes(h, copertura, color = modello)) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  annotate("text", x = 0, y = 0.928, label = "nominale 90%",
           color = "grey40", size = 2.7, hjust = 0) +
  geom_line(linewidth = 0.8) +
  geom_point(data = subset(dA, modello == "LKJ"), size = 1.1) +
  scale_color_manual(values = c(strutturata = col_str, LKJ = col_lkj)) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  scale_y_continuous(limits = c(0.4, 1), breaks = seq(0.4, 1, 0.1)) +
  labs(title = "Copertura per orizzonte (c = 1)",
       x = "orizzonte h", y = "copertura empirica") +
  tema

dir.create(FIG, showWarnings = FALSE)

# --- Figura divergenze per replica -------------------------------------------
dg <- readRDS("output/diagnostiche_E1.rds")
dv <- rbind(data.frame(modello = "LKJ", div = dg$lkj$repliche$div),
            data.frame(modello = "strutturata", div = dg$strutturata$repliche$div))
dv$modello <- factor(dv$modello, levels = c("strutturata", "LKJ"))
p_div <- ggplot(dv, aes(div)) +
  geom_bar(fill = "#8c8c8c", width = 0.7) +
  facet_wrap(~ modello) +
  scale_x_continuous(breaks = 0:10) +
  labs(x = "transizioni divergenti", y = "numero di repliche") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        strip.text = element_text(size = 10.5, color = "grey20"))
ggsave(file.path(FIG, "fig_E1_divergenze.pdf"),
       p_div, width = 8, height = 3)
cat("figura divergenze scritta\n")

# --- Figura Pareto-k per replica ----------------------------------------------
kk <- rbind(
  data.frame(modello = "strutturata",
             c = rep(c("0,9","0,8","0,7"), each = 100),
             k = as.vector(dg$strutturata$k_repliche)),
  data.frame(modello = "LKJ",
             c = rep(c("0,9","0,8","0,7"), each = 100),
             k = as.vector(dg$lkj$k_repliche)))
kk$modello <- factor(kk$modello, levels = c("strutturata","LKJ"))
kk$c <- factor(kk$c, levels = c("0,9","0,8","0,7"))
p_k <- ggplot(kk, aes(c, k)) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "#b03a2e", linewidth = 0.4) +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "grey45", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.45, color = "#4a4a4a") +
  annotate("text", x = 0.55, y = 0.74, label = "0,7", color = "#b03a2e", size = 2.8, hjust = 0) +
  annotate("text", x = 0.55, y = 0.54, label = "0,5", color = "grey45", size = 2.8, hjust = 0) +
  facet_wrap(~ modello) +
  labs(x = "learning rate c bersaglio", y = expression(hat(k))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
        strip.text = element_text(size = 10.5, color = "grey20"))
ggsave(file.path(FIG, "fig_E1_pareto_k.pdf"),
       p_k, width = 8, height = 3.2)
cat("figura pareto-k scritta\n")

# --- Figura copertura per orizzonte, standalone --------------------------------
ggsave(file.path(FIG, "fig_E1_copertura.pdf"),
       pA + labs(title = NULL), width = 8, height = 3.6)
cat("figura copertura standalone scritta\n")

# --- Figura ampiezza + RMSE per orizzonte --------------------------------------
dC <- rbind(
  data.frame(h = h, val = dg$strutturata$calibrazione$`1`$amp_h,  serie = "ampiezza", modello = "strutturata"),
  data.frame(h = h, val = dg$lkj$calibrazione$`1`$amp_h,          serie = "ampiezza", modello = "LKJ"),
  data.frame(h = h, val = dg$strutturata$calibrazione$`1`$rmse_h, serie = "RMSE", modello = "strutturata"),
  data.frame(h = h, val = dg$lkj$calibrazione$`1`$rmse_h,         serie = "RMSE", modello = "LKJ"))
dC$modello <- factor(dC$modello, levels = c("strutturata","LKJ"))
p_ar <- ggplot(dC, aes(h, val, color = modello)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.3) +
  facet_wrap(~ serie, scales = "free_y") +
  scale_color_manual(values = c(strutturata = "#1a1a1a", LKJ = "#8c8c8c")) +
  scale_x_continuous(breaks = seq(0, 16, 4)) +
  labs(x = "orizzonte h", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(size = 10.5, color = "grey20"))
ggsave(file.path(FIG, "fig_E1_amp_rmse.pdf"),
       p_ar, width = 8, height = 3.4)
cat("figura ampiezza/RMSE scritta\n")
