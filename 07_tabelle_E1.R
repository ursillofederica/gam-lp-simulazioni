# Tabelle LaTeX per il § 5.2 (E1), generate dai dati (riproducibili).
# Output in TAB: tab_e1_{protocollo,diagnostiche,ripesatura,cella,profilo}.tex
# Eseguire dalla cartella cap5.
TAB <- "../../../05_tesi/tabelle"    # destinazione delle tabelle di tesi
if (!dir.exists(TAB)) { TAB <- "tabelle"; dir.create(TAB, showWarnings = FALSE) }   # fuori dal progetto tesi
diag <- readRDS("output/diagnostiche_E1.rds")
dir.create(TAB, showWarnings = FALSE)
it <- function(x, d = 3) gsub("\\.", "{,}", formatC(x, format = "f", digits = d))

## Tabella 1: protocollo di ristima e certificazione
riga <- function(nome, m) {
  ct <- diag[[m]]$cert
  paste(nome, ct$R, ct$rifit_init, ct$rifit_div, ct$div_max_pre,
        paste0("$", it(ct$rhat_max), "$"), ct$div_max_post, ct$ess_min,
        paste0("$", it(ct$quota_div0, 2), "$"), sep = " & ")
}
writeLines(c(
  "\\begin{tabular}{lcccccccc}", "\\toprule",
  "modello & $R$ & \\multicolumn{2}{c}{ristime} & div.\\ max & $\\widehat{R}$ max & div.\\ max & ESS & quota \\\\",
  " & & init & div. & (pre) & (post) & (post) & min & div.\\ $=0$ \\\\",
  "\\midrule",
  paste0(riga("LKJ", "lkj"), " \\\\"),
  paste0(riga("strutturata", "strutturata"), " \\\\"),
  "\\bottomrule", "\\end{tabular}"),
  file.path(TAB, "tab_e1_protocollo.tex"))

## Tabella 2: distribuzioni delle tre diagnostiche
qr <- function(x, d) sapply(quantile(x, c(0, .25, .5, .75, 1)),
                            function(v) paste0("$", it(v, d), "$"))
qi <- function(x) sapply(quantile(x, c(0, .25, .5, .75, 1)), function(v) round(v))
blocco <- function(m, nome) {
  r <- diag[[m]]$repliche
  c(paste(nome, "& $\\widehat{R}$ &", paste(qr(r$rhat, 3), collapse = " & "), "\\\\"),
    paste(" & divergenze &", paste(qi(r$div), collapse = " & "), "\\\\"),
    paste(" & ESS bulk &", paste(qi(r$ess), collapse = " & "), "\\\\"))
}
writeLines(c(
  "\\begin{tabular}{llccccc}", "\\toprule",
  "modello & diagnostica & min & $q_{25}$ & mediana & $q_{75}$ & max \\\\",
  "\\midrule",
  blocco("lkj", "LKJ"), "\\midrule",
  blocco("strutturata", "strutturata"),
  "\\bottomrule", "\\end{tabular}"),
  file.path(TAB, "tab_e1_diagnostiche.tex"))
cat("tabelle E1 scritte\n")

## Tabella 3: affidabilita' della ripesatura
rr <- function(m, nome) {
  rp <- diag[[m]]$ripesatura
  sapply(1:nrow(rp), function(i) paste(
    if (i == 1) nome else "", paste0("$", it(rp$c[i], 1), "$"),
    paste0("$", it(rp$k_med[i], 2), "$"), paste0("$", it(rp$k_max[i], 2), "$"),
    paste0("$", it(rp$quota_k_alto[i], 2), "$"),
    paste0("$", it(rp$ess_rel_med[i], 3), "$"), sep = " & "))
}
writeLines(c(
  "\\begin{tabular}{llcccc}", "\\toprule",
  "modello & $c$ & $\\hat{k}$ mediano & $\\hat{k}$ max & quota $\\hat{k} > 0{,}7$ & ESS rel.\\ mediano \\\\",
  "\\midrule",
  paste0(rr("lkj", "LKJ"), " \\\\"), "\\midrule",
  paste0(rr("strutturata", "strutturata"), " \\\\"),
  "\\bottomrule", "\\end{tabular}"),
  file.path(TAB, "tab_e1_ripesatura.tex"))
cat("tabella ripesatura scritta\n")

## Tabella 4: vista di cella (calibrazione)
cl <- function(cc) {
  s <- diag$strutturata$calibrazione[[cc]]; l <- diag$lkj$calibrazione[[cc]]
  paste(paste0("$", gsub("\\.", "{,}", cc), "$"),
        paste0("$", it(mean(s$cop_h)), "$"), paste0("$", it(mean(l$cop_h)), "$"),
        paste0("$", it(mean(s$rmse_h)), "$"), paste0("$", it(mean(l$rmse_h)), "$"),
        paste0("$", it(mean(s$amp_h)), "$"), paste0("$", it(mean(l$amp_h)), "$"), sep = " & ")
}
writeLines(c(
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{copertura} & \\multicolumn{2}{c}{RMSE} & \\multicolumn{2}{c}{ampiezza} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "$c$ & strutt. & LKJ & strutt. & LKJ & strutt. & LKJ \\\\",
  "\\midrule",
  paste0(cl("1"), " \\\\"), paste0(cl("0.9"), " \\\\"),
  "\\bottomrule", "\\end{tabular}"),
  file.path(TAB, "tab_e1_cella.tex"))

## Tabella 5: profilo per orizzonte, doppio blocco
blocco_prof <- function(m) {
  x <- diag[[m]]$calibrazione$`1`
  sapply(1:17, function(i) paste(i - 1,
    paste0("$", it(x$cop_h[i]), "$"), paste0("$", it(x$mcse_h[i]), "$"),
    paste0("$", it(x$bias_h[i]), "$"), paste0("$", it(x$rmse_h[i]), "$"),
    paste0("$", it(x$amp_h[i]), "$"), sep = " & "))
}
writeLines(c(
  "\\begin{tabular}{cccccc}", "\\toprule",
  "$h$ & copertura & MC-SE & distorsione & RMSE & ampiezza \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\emph{strutturata, $c = 1$}} \\\\",
  paste0(blocco_prof("strutturata"), " \\\\"), "\\midrule",
  "\\multicolumn{6}{l}{\\emph{LKJ, $c = 1$}} \\\\",
  paste0(blocco_prof("lkj"), " \\\\"),
  "\\bottomrule", "\\end{tabular}"),
  file.path(TAB, "tab_e1_profilo.tex"))
cat("tabelle calibrazione scritte\n")
