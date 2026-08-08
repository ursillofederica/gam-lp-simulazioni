# Tabelle LaTeX per 5.3 (E2) e 5.4 (E3), generate dai dati.
# Fonte: output/diagnostiche_E2E3.rds (08_diagnostiche_E2E3.R).
d <- readRDS("output/diagnostiche_E2E3.rds")
TAB <- "../../../05_tesi/tabelle"  
if (!dir.exists(TAB)) { TAB <- "tabelle"; dir.create(TAB, showWarnings = FALSE) }   
it <- function(x, k = 3) gsub("\\.", "{,}", formatC(x, format = "f", digits = k))
m_ <- function(x) paste0("$", x, "$")

FONTE <- list(
  E2 = list(strutturata = list("1" = c("strutturata", "1"),
                               "0.9" = c("strutturata", "0.9"),
                               "0.8" = c("strutturata_c80", "0.8"),
                               "0.7" = c("strutturata_c70", "0.7")),
            lkj         = list("1" = c("lkj", "1"),
                               "0.9" = c("lkj_c90", "0.9"),
                               "0.8" = c("lkj_c80", "0.8"),
                               "0.7" = c("lkj_c70", "0.7"))),
  E3 = list(strutturata = list("1" = c("strutturata", "1"),
                               "0.9" = c("strutturata", "0.9"),
                               "0.8" = c("strutturata_c70", "0.8"),
                               "0.7" = c("strutturata_c70", "0.7")),
            lkj         = list("1" = c("lkj", "1"),
                               "0.9" = c("lkj_c80", "0.9"),
                               "0.8" = c("lkj_c80", "0.8"),
                               "0.7" = c("lkj_c70", "0.7"))))
leggi <- function(esp, mod, cc) {
  f <- FONTE[[esp]][[mod]][[cc]]
  d[[esp]][[f[1]]]$calibrazione[[f[2]]]
}
diretta <- function(esp, mod, cc) {
  f <- FONTE[[esp]][[mod]][[cc]]
  d[[esp]][[f[1]]]$cert$base_c == as.numeric(cc)
}
GC <- c("1", "0.9", "0.8", "0.7")

#  Ristime 
reps_backup <- function(dirs, dgp, cella) {
  f <- unlist(lapply(dirs, function(dd)
    list.files(file.path("output/fit", dd),
               pattern = sprintf("^%s_%s_rep_[0-9]+", dgp, cella))))
  length(unique(sub(".*rep_([0-9]+).*", "\\1", f)))
}
DGP <- c(E2 = "persistente", E3 = "nonlineare")

tab_protocollo <- function(esp, celle, etich, file) {
  riga <- function(m, nome) {
    ct <- d[[esp]][[m]]$cert
    n_div  <- reps_backup(c("_backup_div_21lug", "_backup_div_22lug"), DGP[esp], m)
    n_rhat <- ct$rifit_div - n_div
    paste(nome, m_(it(ct$base_c, 1)), ct$R, ct$rifit_init, n_div, n_rhat,
          ifelse(is.na(ct$div_max_pre), "--", ct$div_max_pre),
          m_(it(ct$rhat_max, 4)), ct$div_max_post, ct$ess_min,
          m_(it(ct$quota_div0, 2)), sep = " & ")
  }
  writeLines(c(
    "\\begin{tabular}{lcccccccccc}", "\\toprule",
    "cella & $c$ base & $R$ & \\multicolumn{3}{c}{ristime} & div.\\ max & $\\widehat{R}$ max & div.\\ max & ESS & quota \\\\",
    "\\cmidrule(lr){4-6}",
    " & & & init & div. & $\\widehat{R}$ & (pre) & (post) & (post) & min & div.\\ $=0$ \\\\",
    "\\midrule",
    paste0(mapply(riga, celle, etich), " \\\\"),
    "\\bottomrule", "\\end{tabular}"),
    file.path(TAB, file))
}

tab_protocollo("E2",
  c("strutturata", "strutturata_c80", "strutturata_c70", "lkj", "lkj_c90", "lkj_c80", "lkj_c70"),
  c("strutturata", "strutturata", "strutturata", "LKJ", "LKJ", "LKJ", "LKJ"),
  "tab_e2_protocollo.tex")
tab_protocollo("E3",
  c("strutturata", "strutturata_c70", "lkj", "lkj_c80", "lkj_c70"),
  c("strutturata", "strutturata", "LKJ", "LKJ", "LKJ"),
  "tab_e3_protocollo.tex")

# Tabella: distribuzione delle diagnostiche 
tab_diagnostiche <- function(esp, celle, etich, file) {
  qr <- function(x, dd) sapply(quantile(x, c(0, .25, .5, .75, 1)),
                               function(v) m_(it(v, dd)))
  qi <- function(x) sapply(quantile(x, c(0, .25, .5, .75, 1)), function(v) round(v))
  blocco <- function(m, nome) {
    r <- d[[esp]][[m]]$repliche
    c(paste(nome, "& $\\widehat{R}$ &", paste(qr(r$rhat, 3), collapse = " & "), "\\\\"),
      paste(" & divergenze &", paste(qi(r$div), collapse = " & "), "\\\\"),
      paste(" & ESS bulk &", paste(qi(r$ess), collapse = " & "), "\\\\"))
  }
  writeLines(c(
    "\\begin{tabular}{llccccc}", "\\toprule",
    "cella & diagnostica & min & $q_{25}$ & mediana & $q_{75}$ & max \\\\",
    "\\midrule",
    unlist(lapply(seq_along(celle), function(j)
      c(blocco(celle[j], etich[j]), if (j < length(celle)) "\\midrule"))),
    "\\bottomrule", "\\end{tabular}"), file.path(TAB, file))
}
tab_diagnostiche("E2",
  c("strutturata", "strutturata_c80", "strutturata_c70", "lkj", "lkj_c90", "lkj_c80", "lkj_c70"),
  c("strutt.\\ ($c=1$)", "strutt.\\ ($c=0{,}8$)", "strutt.\\ ($c=0{,}7$)",
    "LKJ ($c=1$)", "LKJ ($c=0{,}9$)", "LKJ ($c=0{,}8$)", "LKJ ($c=0{,}7$)"),
  "tab_e2_diagnostiche.tex")
tab_diagnostiche("E3",
  c("strutturata", "strutturata_c70", "lkj", "lkj_c80", "lkj_c70"),
  c("strutt.\\ ($c=1$)", "strutt.\\ ($c=0{,}7$)",
    "LKJ ($c=1$)", "LKJ ($c=0{,}8$)", "LKJ ($c=0{,}7$)"),
  "tab_e3_diagnostiche.tex")

# Tabella: affidabilità della ripesatura
tab_ripesatura <- function(esp, celle, etich, file) {
  blocco <- function(m, nome) {
    rp <- d[[esp]][[m]]$ripesatura
    sapply(seq_len(nrow(rp)), function(i) paste(
      if (i == 1) nome else "", m_(it(rp$c[i], 1)), m_(it(rp$k_med[i], 2)),
      m_(it(rp$k_max[i], 2)), m_(it(rp$quota_k_alto[i], 2)),
      m_(it(rp$ess_rel_med[i], 3)), sep = " & "))
  }
  writeLines(c(
    "\\begin{tabular}{llcccc}", "\\toprule",
    "cella (base) & $c$ & $\\hat{k}$ mediano & $\\hat{k}$ max & quota $\\hat{k} > 0{,}7$ & ESS rel.\\ mediano \\\\",
    "\\midrule",
    unlist(lapply(seq_along(celle), function(j)
      c(paste0(blocco(celle[j], etich[j]), " \\\\"),
        if (j < length(celle)) "\\midrule"))),
    "\\bottomrule", "\\end{tabular}"), file.path(TAB, file))
}

tab_ripesatura("E2",
  c("lkj", "lkj_c90", "lkj_c80", "lkj_c70",
    "strutturata", "strutturata_c80", "strutturata_c70"),
  c("LKJ ($c=1$)", "LKJ ($c=0{,}9$)", "LKJ ($c=0{,}8$)", "LKJ ($c=0{,}7$)",
    "strutt.\\ ($c=1$)", "strutt.\\ ($c=0{,}8$)", "strutt.\\ ($c=0{,}7$)"),
  "tab_e2_ripesatura.tex")
cat("tabella ripesatura E2 scritta\n")
tab_ripesatura("E3",
  c("lkj","lkj_c80","lkj_c70","strutturata","strutturata_c70"),
  c("LKJ ($c=1$)","LKJ ($c=0{,}8$)","LKJ ($c=0{,}7$)","strutt.\\ ($c=1$)","strutt.\\ ($c=0{,}7$)"),
  "tab_e3_ripesatura.tex")

# Tabella: "convalida incrociata" (solo E2: piu' strade allo stesso c) ---
# Per ogni (modello, c bersaglio) elenca le vie valide: fit diretto oppure
# ripesatura con k_max <= 0.7 e la copertura che ciascuna restituisce.
tab_crossval <- function(file) {
  celle_mod <- list(LKJ = c("lkj", "lkj_c90", "lkj_c80", "lkj_c70"),
                    strutturata = c("strutturata", "strutturata_c80", "strutturata_c70"))
  righe <- character(); primo_mod <- TRUE
  for (mod in names(celle_mod)) {
    if (!primo_mod) righe <- c(righe, "\\midrule")
    primo_mod <- FALSE; prima_riga_mod <- TRUE; primo_gruppo <- TRUE
    for (cc in c(1, 0.9, 0.8, 0.7)) {
      vie <- list()
      for (m in celle_mod[[mod]]) {
        bc <- d$E2[[m]]$cert$base_c
        cop <- mean(d$E2[[m]]$calibrazione[[as.character(cc)]]$cop_h)
        if (bc == cc)
          vie[[length(vie) + 1]] <- list(base = bc, cop = cop, via = "fit diretto", ord = 2)
        else {
          rp <- d$E2[[m]]$ripesatura; km <- rp$k_max[rp$c == cc]
          if (length(km) && km <= 0.7)
            vie[[length(vie) + 1]] <- list(base = bc, cop = cop,
              via = paste0("ripesatura ($\\hat{k}=", it(km, 2), "$)"), ord = 1)
        }
      }
      vie <- vie[order(-sapply(vie, `[[`, "ord"), -sapply(vie, `[[`, "base"))]
      if (!primo_gruppo) righe <- c(righe, "\\addlinespace")
      primo_gruppo <- FALSE
      for (j in seq_along(vie)) {
        v <- vie[[j]]
        col_mod <- if (prima_riga_mod) mod else ""
        col_c <- if (j == 1) m_(it(cc, 1)) else ""
        righe <- c(righe, paste0(paste(col_mod, col_c, m_(it(v$base, 1)),
          m_(it(v$cop, 3)), v$via, sep = " & "), " \\\\"))
        prima_riga_mod <- FALSE
      }
    }
  }
  writeLines(c("\\begin{tabular}{llccl}", "\\toprule",
    "modello & $c$ bersaglio & base $c$ & copertura & fonte \\\\", "\\midrule",
    righe, "\\bottomrule", "\\end{tabular}"), file.path(TAB, file))
}
tab_crossval("tab_e2_crossval.tex")

# E2: curva in c (copertura, RMSE, ampiezza per i due modelli) 
riga_c <- function(cc) {
  s <- leggi("E2", "strutturata", cc); l <- leggi("E2", "lkj", cc)
  marca <- function(esp, mod) if (diretta(esp, mod, cc)) "$^{\\dagger}$" else ""
  paste(m_(gsub("\\.", "{,}", cc)),
        paste0(m_(it(mean(s$cop_h))), marca("E2", "strutturata")),
        paste0(m_(it(mean(l$cop_h))), marca("E2", "lkj")),
        m_(it(mean(s$rmse_h))), m_(it(mean(l$rmse_h))),
        m_(it(mean(s$amp_h))), m_(it(mean(l$amp_h))), sep = " & ")
}
writeLines(c(
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{copertura} & \\multicolumn{2}{c}{RMSE} & \\multicolumn{2}{c}{ampiezza} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "$c$ & strutt. & LKJ & strutt. & LKJ & strutt. & LKJ \\\\",
  "\\midrule",
  paste0(sapply(GC, riga_c), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\\\[2pt]\\footnotesize $^{\\dagger}$ cella stimata direttamente alla potenza $c$; le altre voci sono ripesature con $\\hat{k} \\le 0{,}7$."),
  file.path(TAB, "tab_e2_curva.tex"))

# E2: profilo per orizzonte, i due modelli a c = 1 e c = 0,7
blocco_prof <- function(esp, mod, cc) {
  x <- leggi(esp, mod, cc)
  sapply(seq_along(x$cop_h), function(i) paste(i - 1,
    m_(it(x$cop_h[i])), m_(it(x$mcse_h[i])), m_(it(x$bias_h[i])),
    m_(it(x$rmse_h[i])), m_(it(x$amp_h[i])), sep = " & "))
}
writeLines(c(
  "\\begin{tabular}{cccccc}", "\\toprule",
  "$h$ & copertura & MC-SE & distorsione & RMSE & ampiezza \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\emph{strutturata, $c = 1$}} \\\\",
  paste0(blocco_prof("E2", "strutturata", "1"), " \\\\"), "\\midrule",
  "\\multicolumn{6}{l}{\\emph{strutturata, $c = 0{,}7$}} \\\\",
  paste0(blocco_prof("E2", "strutturata", "0.7"), " \\\\"), "\\midrule",
  "\\multicolumn{6}{l}{\\emph{LKJ, $c = 1$}} \\\\",
  paste0(blocco_prof("E2", "lkj", "1"), " \\\\"), "\\midrule",
  "\\multicolumn{6}{l}{\\emph{LKJ, $c = 0{,}7$}} \\\\",
  paste0(blocco_prof("E2", "lkj", "0.7"), " \\\\"),
  "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_e2_profilo.tex"))

# E3: cella, due contrasti 
riga_c3 <- function(cc) {
  s <- leggi("E3", "strutturata", cc); l <- leggi("E3", "lkj", cc)
  marca <- function(mod) if (diretta("E3", mod, cc)) "$^{\\dagger}$" else ""
  paste(m_(gsub("\\.", "{,}", cc)),
        paste0(m_(it(mean(s$pos$cop_h))), marca("strutturata")),
        m_(it(mean(s$neg$cop_h))),
        paste0(m_(it(mean(l$pos$cop_h))), marca("lkj")),
        m_(it(mean(l$neg$cop_h))),
        m_(it(mean(s$pos$amp_h))), m_(it(mean(l$pos$amp_h))), sep = " & ")
}
writeLines(c(
  "\\begin{tabular}{lcccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{strutturata} & \\multicolumn{2}{c}{LKJ} & \\multicolumn{2}{c}{ampiezza} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
  "$c$ & $\\delta = +2$ & $\\delta = -2$ & $\\delta = +2$ & $\\delta = -2$ & strutt. & LKJ \\\\",
  "\\midrule",
  paste0(sapply(GC, riga_c3), " \\\\"),
  "\\bottomrule", "\\end{tabular}",
  "\\\\[2pt]\\footnotesize L'ampiezza delle bande è identica nei due contrasti per costruzione. $^{\\dagger}$ cella stimata direttamente alla potenza $c$."),
  file.path(TAB, "tab_e3_cella.tex"))

# E3: profilo per orizzonte nei due contrasti, c = 1 
blocco_prof3 <- function(mod, cc) {
  x <- leggi("E3", mod, cc)
  sapply(seq_along(x$pos$cop_h), function(i) paste(i - 1,
    m_(it(x$pos$cop_h[i])), m_(it(x$pos$bias_h[i])),
    m_(it(x$neg$cop_h[i])), m_(it(x$neg$bias_h[i])),
    m_(it(x$pos$amp_h[i])), sep = " & "))
}
writeLines(c(
  "\\begin{tabular}{cccccc}", "\\toprule",
  " & \\multicolumn{2}{c}{$\\delta = +2$} & \\multicolumn{2}{c}{$\\delta = -2$} & \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "$h$ & copertura & distorsione & copertura & distorsione & ampiezza \\\\",
  "\\midrule",
  "\\multicolumn{6}{l}{\\emph{strutturata, $c = 1$}} \\\\",
  paste0(blocco_prof3("strutturata", "1"), " \\\\"), "\\midrule",
  "\\multicolumn{6}{l}{\\emph{strutturata, $c = 0{,}7$}} \\\\",
  paste0(blocco_prof3("strutturata", "0.7"), " \\\\"), "\\midrule",
  "\\multicolumn{6}{l}{\\emph{LKJ, $c = 1$}} \\\\",
  paste0(blocco_prof3("lkj", "1"), " \\\\"),
  "\\bottomrule", "\\end{tabular}"), file.path(TAB, "tab_e3_profilo.tex"))

