# Grafico 1 del recupero tanh: i due contrasti di taglia (delta = 1 e delta = 2)
#  con mediana a posteriori, banda al 90% e valore vero.
# Fonte: output/reperti_E3.rds (nessun rifit). Output: fig_E3_tanh_unico.pdf
# Eseguire dalla cartella degli script.
library(ggplot2)
FIG <- "../../../05_tesi/figure"    
if (!dir.exists(FIG)) { FIG <- "figure"; dir.create(FIG, showWarnings = FALSE) }   
h <- 0:16
col1 <- "#1a1a1a"   # delta = 1
col2 <- "#e67e22"   # delta = 2

TH <- readRDS("output/reperti_E3.rds")$tanh
dd <- rbind(
  data.frame(h = h, med = TH$d1$med, lo = TH$d1$lo, hi = TH$d1$hi,
             vera = TH$vera$d1, contrasto = "delta == 1"),
  data.frame(h = h, med = TH$d2$med, lo = TH$d2$lo, hi = TH$d2$hi,
             vera = TH$vera$d2, contrasto = "delta == 2"))
dd$contrasto <- factor(dd$contrasto, levels = c("delta == 1", "delta == 2"))

p <- ggplot(dd, aes(h, group = contrasto)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = contrasto), alpha = 0.22) +
  geom_line(aes(y = med, color = contrasto), linewidth = 0.9) +
  geom_line(aes(y = vera, color = contrasto), linewidth = 0.5, linetype = "22") +
  geom_point(aes(y = vera, color = contrasto), size = 1.5, shape = 16) +
  scale_color_manual(values = c("delta == 1" = col1, "delta == 2" = col2),
                     labels = c(expression(delta == 1), expression(delta == 2))) +
  scale_fill_manual(values = c("delta == 1" = col1, "delta == 2" = col2),
                    labels = c(expression(delta == 1), expression(delta == 2))) +
  scale_x_continuous(breaks = seq(0, 16, 2)) +
  annotate("text", x = 11.5, y = 0.42, hjust = 0, size = 2.9, color = "grey30",
           label = "linea piena: mediana a posteriori") +
  annotate("text", x = 11.5, y = 0.24, hjust = 0, size = 2.9, color = "grey30",
           label = "punti + tratteggio: risposta vera") +
  labs(x = "orizzonte h", y = "risposta d'impulso",
       color = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        axis.title = element_text(size = 10, color = "grey20"),
        axis.text = element_text(color = "grey20"),
        legend.position = c(0.12, 0.86),
        legend.background = element_rect(fill = "white", color = NA))
ggsave(file.path(FIG, "fig_E3_tanh_unico.pdf"), p, width = 6.4, height = 4)

