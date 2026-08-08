# Estrazione della superficie stimata f(s,h) su una griglia di s per il grafico
# 3D del recupero tanh. Il fit salvato tiene solo i due contrasti: qui si
# riesegue UNA volta il fit tanh buono (stesso seed del reperto, seed=1500)
# solo per estrarre i coefficienti gamma. NON sovrascrive il fit su disco:
# salva la superficie in output/superficie_tanh_griglia.rds.
# Eseguire dalla cartella degli script.
source("02_fitta.R")   # basi (B_H, nodi_s, J_S, K_H), MODELLI, costanti

d <- readRDS(file.path(DIR_OUT, "tanh", "rep_001.rds"))
sd_stan <- costruisci_stan_data_superficie(d, deltas = c(2, 1))
fit <- MODELLI$superficie_strutturata$sample(
  data = sd_stan, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000,
  adapt_delta = 0.95, seed = 1500, init = 0.1,
  refresh = 0, show_messages = FALSE)

cat("estratto fit: rhat", round(max(fit$summary("irf_pos")$rhat, na.rm = TRUE), 3),
    "| div", sum(fit$diagnostic_summary()$num_divergent), "\n")

# gamma: media a posteriori, riordinata in matrice K x J (col-major = (j-1)*K+k)
gamma_mean <- colMeans(fit$draws("gamma", format = "matrix"))
G <- matrix(gamma_mean, nrow = K_H, ncol = J_S)   # G[k,j] = gamma[(j-1)*K+k]

# base su s ricostruita dai dati e valutata su una griglia fitta
B_S <- splines::bs(d$s, knots = nodi_s, degree = 3, intercept = TRUE,
                   Boundary.knots = c(-4, 4))
s_grid <- seq(-4, 4, length.out = 60)
BSg <- predict(B_S, s_grid)                       # 60 x J_S

# superficie stimata e vera su (h, s): [H1 x 60]
surf_est  <- B_H %*% G %*% t(BSg)
ff <- function(s, hh) 1.5 * hh * exp(-0.2 * hh) * tanh(0.5 * s)
surf_true <- outer(0:H, s_grid, function(hh, s) ff(s, hh))

saveRDS(list(h = 0:H, s = s_grid, est = surf_est, vera = surf_true),
        "output/superficie_tanh_griglia.rds")
cat("range superficie stimata:", round(range(surf_est), 2),
    "| vera:", round(range(surf_true), 2), "\n")
cat("scarto massimo stimata-vera:", round(max(abs(surf_est - surf_true)), 3), "\n")
cat("Salvato output/superficie_tanh_griglia.rds\n")
