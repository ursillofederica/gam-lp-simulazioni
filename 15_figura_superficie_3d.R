# Grafico 2 del recupero tanh: la superficie f(s,h) in 3D, vera contro stimata.
# Fonte: output/superficie_tanh_griglia.rds (14_superficie_tanh_3d.R).
# Output: fig_E3_superficie_3d.pdf (persp base R, due pannelli).
# Eseguire dalla cartella degli script.
FIG <- "../../../05_tesi/figure"    # destinazione delle figure di tesi
if (!dir.exists(FIG)) { FIG <- "figure"; dir.create(FIG, showWarnings = FALSE) }   # fuori dal progetto tesi
g <- readRDS("output/superficie_tanh_griglia.rds")

# persp vuole z con righe = x (s) e colonne = y (h): traspongo (dati sono h x s)
z_true <- t(g$vera); z_est <- t(g$est)
zlim <- range(z_true, z_est)

# colore per faccetta secondo l'altezza media (grigio -> arancio)
pal <- colorRampPalette(c("#f0f0f0", "#f6c28b", "#e67e22"))(64)
faccia <- function(z) {
  nr <- nrow(z); nc <- ncol(z)
  zc <- (z[-1, -1] + z[-1, -nc] + z[-nr, -1] + z[-nr, -nc]) / 4
  pal[cut(zc, breaks = seq(zlim[1], zlim[2], length.out = 65), include.lowest = TRUE)]
}

pdf(file.path(FIG, "fig_E3_superficie_3d.pdf"), width = 9, height = 4.4)
op <- par(mfrow = c(1, 2), mar = c(1.2, 1.2, 2.2, 1.2))
for (cfg in list(list(z = z_true, t = "Superficie vera  f(s,h)"),
                 list(z = z_est,  t = "Superficie stimata (media a posteriori)"))) {
  persp(x = g$s, y = g$h, z = cfg$z, zlim = zlim,
        theta = -38, phi = 22, expand = 0.62, col = faccia(cfg$z),
        border = "grey55", lwd = 0.3, shade = 0.35, ticktype = "detailed",
        nticks = 4, xlab = "shock s", ylab = "orizzonte h", zlab = "risposta",
        cex.axis = 0.6, cex.lab = 0.8, main = cfg$t, cex.main = 0.95)
}
par(op); dev.off()
cat("fig_E3_superficie_3d scritta | scarto max stimata-vera:",
    round(max(abs(g$est - g$vera)), 3), "\n")
