# ==============================================================================
# 00_calibrazione_riferimento.R — calibrazione UNA TANTUM dei phi_p (Processo 1)
#
# Serie: GDPC1 (PIL reale USA, trimestrale, FRED), trasformata in crescita:
#   g_t = 100 * diff(log(GDPC1))
# Campione: 1960Q1-2019Q4 (esclusi i trimestri COVID: outlier che
#   distorcerebbero l'AR — scelta dichiarata).
# Stima: AR(5) via OLS con intercetta; si conservano i soli phi_1..phi_5
#   (il DGP e' a media zero, l'intercetta si scarta).
# Output: output/phi_calibrati.rds (coefficienti + provenienza, CONGELATI)
#         output/GDPC1_raw.csv (dati grezzi congelati per riproducibilita')
# ==============================================================================

url_fred <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=GDPC1"
raw_path <- file.path("output", "GDPC1_raw.csv")
dir.create("output", showWarnings = FALSE)

if (!file.exists(raw_path)) {
    download.file(url_fred, raw_path, quiet = TRUE)   # congelato al primo run
}
dati <- read.csv(raw_path)
names(dati) <- c("data", "gdpc1")
dati$data <- as.Date(dati$data)
dati <- dati[dati$data >= as.Date("1959-10-01") &
             dati$data <= as.Date("2019-12-31"), ]   # un trimestre extra per il diff

g <- 100 * diff(log(dati$gdpc1))                      # crescita trimestrale, %

# AR(5) via OLS
n  <- length(g)
X  <- sapply(1:5, function(p) g[(6 - p):(n - p)])
gy <- g[6:n]
fit <- lm(gy ~ X)
phi <- unname(coef(fit)[-1])

# Controllo di stazionarieta': radici del polinomio caratteristico
radici <- polyroot(c(1, -phi))
stopifnot(all(Mod(radici) > 1))

cal <- list(
    phi            = phi,
    serie          = "GDPC1 (FRED)",
    trasformazione = "100 * diff(log)",
    campione       = paste(format(min(dati$data)), "-", format(max(dati$data))),
    n_oss          = length(gy),
    data_download  = format(Sys.Date()),
    radici_mod_min = min(Mod(radici))
)
saveRDS(cal, file.path("output", "phi_calibrati.rds"))

cat("phi calibrati:", round(phi, 4), "\n")
cat("modulo minimo radici:", round(min(Mod(radici)), 3), "(>1 = stazionario)\n")
cat("n oss:", length(gy), "| campione:", cal$campione, "\n")
