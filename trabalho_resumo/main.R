# =============================================================================
# main.R — Análise principal
# Classificação bayesiana de proficiência em Biologia (níveis tipo PAEBES)
# JPGTEIXEIRA — EST613 Estatística Bayesiana / UFV
#
# Antes: rodar simulacao.R uma vez (gera a pasta data/).
# Com dados reais, basta trocar data/alunos.csv e data/notas.csv
# (mesmo formato); as partes que usam data/verdade.csv deixam de existir.
#
# Saídas (pasta results/): figuras .png, tabelas .csv e os ajustes do Stan (.rds)
# =============================================================================

library(ggplot2)
library(rstan)
library(coda)       # diagnósticos de convergência (como no Módulo V/VI)
library(loo)        # WAIC e LOO
library(bayesplot)  # gráficos da verificação preditiva
source("funcoes.R")
dir.create("results", showWarnings = FALSE)
options(mc.cores = 4)        # uma cadeia por núcleo
rstan_options(auto_write = TRUE)

# -----------------------------------------------------------------------------
# 1. Leitura dos dados
# -----------------------------------------------------------------------------

alunos  <- read.csv("data/alunos.csv")
notas   <- read.csv("data/notas.csv")
verdade <- read.csv("data/verdade.csv")

alunos <- alunos[order(alunos$id), ]
notas  <- notas[order(notas$id, notas$aval), ]

# -----------------------------------------------------------------------------
# 2. Parâmetros fixos do M1 (escolhidos antes de ver os resultados)
# -----------------------------------------------------------------------------

sd_prior <- 10   # incerteza da nota do histórico como prior
sd_prova <- 6    # erro de medida de uma avaliação
sd_ruido <- 3    # quanto o aluno pode mudar entre duas avaliações
limiar   <- 0.7  # probabilidade mínima para trocar o aluno de nível

# -----------------------------------------------------------------------------
# 3. M0 e M1
# -----------------------------------------------------------------------------

post_M1 <- m1_todos(alunos, notas, sd_prior, sd_prova, sd_ruido)

# -----------------------------------------------------------------------------
# 4. M2 em Stan — o caminho até o modelo final
#    Cada tentativa muda UMA coisa em relação à anterior.
#    Os ajustes ficam salvos em results/ (.rds); para rodar de novo, apague o .rds.
#    Verdade da simulação, para comparar: PR = -10, PM = 0, PF = +8; sigma_hist = 5.
#    Histórico: só é aprovado quem tem >= 60; quem chega com 60 foi "aprovado no
#    limite" (histórico verdadeiro <= 60). Todas as tentativas já tratam isso como
#    dado censurado (censura = 1); a seção 4b mostra o que acontece sem isso.
# -----------------------------------------------------------------------------

modelo_M2 <- stan_model(model_code = MODELO_STAN)
d_stan    <- dados_stan(alunos, notas)     # dados no formato do Stan (usado várias vezes)

# Checagem preditiva A PRIORI (Tópico 9): antes de ajustar, as priors geram
# notas plausíveis? Simula 1.000 "anos de notas" só a partir das priors.
set.seed(2026)
y_prior <- prior_preditiva(d_stan, n = 1000)
prior_check <- data.frame(
  estatistica = c("% de notas < 0", "% de notas > 100", "percentil 1%", "mediana", "percentil 99%"),
  valor = round(c(100 * mean(y_prior < 0), 100 * mean(y_prior > 100),
                  quantile(y_prior, c(0.01, 0.5, 0.99))), 1))
cat("\n===== Checagem preditiva a priori =====\n")
print(prior_check, row.names = FALSE)
write.csv(prior_check, "results/prior_preditiva.csv", row.names = FALSE)

ajustar_M2 <- function(arquivo, prior_tau = 1, ancora = 1, prior_sigma_hist = 2,
                       prof_nc = 1, com_efeitos = 1, censura = 1,
                       guardar = c("theta", "gamma", "delta_o", "delta_a", "tau_esc",
                                   "tau_prof", "omega", "sigma", "sigma_hist")) {
  if (file.exists(arquivo)) return(readRDS(arquivo))
  fit <- sampling(modelo_M2,
                  data = dados_stan(alunos, notas, prior_tau, ancora, prior_sigma_hist,
                                    prof_nc, com_efeitos, censura),
                  chains = 4, iter = 2000, warmup = 1000, seed = 2026,
                  control = list(adapt_delta = 0.99),
                  pars = guardar)          # guarda só o necessário (arquivo menor)
  saveRDS(fit, arquivo)
  fit
}
so_efeitos <- c("gamma", "delta_a", "tau_esc", "tau_prof", "omega", "sigma", "sigma_hist")

# Tentativa 1: modelo "ingênuo" (sem âncora; sigma_hist ~ half-normal)
fit_t1 <- ajustar_M2("results/fit_M2_t1.rds", ancora = 0, prior_sigma_hist = 1, guardar = so_efeitos)

# Tentativa 2: + âncora (professores atuais somam zero; delta centrado)
fit_t2 <- ajustar_M2("results/fit_M2_t2.rds", ancora = 1, prior_sigma_hist = 1, prof_nc = 0,
                     guardar = so_efeitos)

# Tentativa 3: + prior Gama(2; 0,2) para sigma_hist
fit_t3 <- ajustar_M2("results/fit_M2_t3.rds", prof_nc = 0, guardar = so_efeitos)
# ... e a mesma tentativa 3 com a half-Cauchy (sensibilidade): uma cadeia trava
fit_t3hc <- ajustar_M2("results/fit_M2_t3_halfcauchy.rds", prior_tau = 2, prof_nc = 0,
                       guardar = so_efeitos)

# Tentativa 4 (FINAL): delta dos professores atuais não centrado
fit_hn <- ajustar_M2("results/fit_M2_halfnormal.rds")                  # principal
fit_hc <- ajustar_M2("results/fit_M2_halfcauchy.rds", prior_tau = 2)   # sensibilidade

tentativas <- rbind(
  diagnostico_curto(fit_t1,   "1: ingênuo (half-normal)", d_stan$nomes_pa),
  diagnostico_curto(fit_t2,   "2: + âncora (half-normal)", d_stan$nomes_pa),
  diagnostico_curto(fit_t3,   "3: + Gama em sigma_hist (half-normal)", d_stan$nomes_pa),
  diagnostico_curto(fit_t3hc, "3: idem, com half-Cauchy", d_stan$nomes_pa),
  diagnostico_curto(fit_hn,   "4: + delta não centrado (half-normal) FINAL", d_stan$nomes_pa),
  diagnostico_curto(fit_hc,   "4: idem, com half-Cauchy", d_stan$nomes_pa)
)
cat("\n===== Caminho até o modelo final (M2) =====\n")
print(tentativas, row.names = FALSE)
write.csv(tentativas, "results/tentativas_M2.csv", row.names = FALSE)

# Comparação de modelos: o mesmo modelo final SEM escola e professor (seção 9)
fit_se <- ajustar_M2("results/fit_M2_semefeitos.rds", com_efeitos = 0)

# 4b. O modelo final IGNORANDO a regra dos 60 ("ingênuo": o 60 vale como nota exata)
fit_ing <- ajustar_M2("results/fit_M2_ingenuo.rds", censura = 0)
cat("\nAlunos com histórico = 60 (aprovados no limite):", sum(alunos$nota_hist <= 60), "de", nrow(alunos), "\n")

post_M2     <- m2_posterior(fit_hn, alunos)
post_M2_hc  <- m2_posterior(fit_hc, alunos)
post_M2_ing <- m2_posterior(fit_ing, alunos)

# -----------------------------------------------------------------------------
# 5. Classificação: uma linha por aluno x avaliação
# -----------------------------------------------------------------------------

res <- merge(notas, verdade[, c("id", "aval", "entrada", "trajetoria", "theta_real")],
             by = c("id", "aval"))
res <- merge(res, alunos[, c("id", "escola_origem", "prof_atual", "nota_hist")], by = "id")
res <- res[order(res$id, res$aval), ]

res$nivel_real <- as.character(nota_para_nivel(res$theta_real))
res$m0_nota    <- ave(res$nota, res$id, FUN = m0_media_simples)
res$M0         <- as.character(nota_para_nivel(res$m0_nota))

c1 <- classificar(post_M1, alunos, limiar)
c2 <- classificar(post_M2, alunos, limiar)
c3 <- classificar(post_M2_hc, alunos, limiar)
c4 <- classificar(post_M2_ing, alunos, limiar)
res$M1_maior_prob <- c1$maior_prob;  res$M1_limiar <- c1$limiar
res$M2_maior_prob <- c2$maior_prob;  res$M2_limiar <- c2$limiar
res$M2hc_limiar   <- c3$limiar       # sensibilidade (half-Cauchy)
res$M2ing_limiar  <- c4$limiar       # ignorando a censura do histórico

# Médias da posterior, para as figuras
res$M1_media <- post_M1$media[order(post_M1$id, post_M1$aval)]
res$M1_s     <- post_M1$s[order(post_M1$id, post_M1$aval)]
res$M2_media <- post_M2$media
res$M2_s     <- post_M2$s

modelos <- c("M0", "M1_maior_prob", "M1_limiar", "M2_maior_prob", "M2_limiar", "M2hc_limiar",
             "M2ing_limiar")

# -----------------------------------------------------------------------------
# 6. Métricas de classificação
# -----------------------------------------------------------------------------

acerto <- function(mod, grupo = NULL) {
  certo <- res[[mod]] == res$nivel_real
  if (is.null(grupo)) mean(certo) else tapply(certo, res[[grupo]], mean)
}

acerto_aval <- do.call(rbind, lapply(modelos, function(mod) {
  data.frame(modelo = mod, aval = 1:9, acerto = acerto(mod, "aval"))
}))
acerto_trajetoria <- sapply(modelos, acerto, grupo = "trajetoria")
acerto_professor  <- sapply(modelos, acerto, grupo = "prof_atual")
res$historico <- ifelse(res$nota_hist <= 60, "hist = 60 (censurado)", "hist > 60")
acerto_historico <- sapply(modelos, acerto, grupo = "historico")

# Trocas de nível ("piscadas") e atraso na detecção, por aluno
por_aluno <- do.call(rbind, lapply(split(res, res$id), function(d) {
  linha <- data.frame(id = d$id[1], trocas_real = n_trocas(d$nivel_real))
  for (mod in modelos) {
    linha[[paste0("trocas_", mod)]] <- n_trocas(d[[mod]])
    linha[[paste0("atraso_", mod)]] <- atraso_deteccao(d$nivel_real, d[[mod]])
  }
  linha
}))

tabela <- do.call(rbind, lapply(modelos, function(mod) {
  atr <- por_aluno[[paste0("atraso_", mod)]]
  atr <- atr[!is.na(atr)]                          # só alunos que mudaram de nível
  data.frame(modelo = mod,
             acerto_perc   = round(100 * acerto(mod), 1),
             trocas_medias = round(mean(por_aluno[[paste0("trocas_", mod)]]), 2),
             atraso_medio  = round(mean(atr[is.finite(atr)]), 2),
             perc_nao_detectado = round(100 * mean(is.infinite(atr)), 1))
}))

cat("\nTrocas reais de nível por aluno (média):", round(mean(por_aluno$trocas_real), 2), "\n")
cat("Alunos que mudam de nível de verdade:", sum(!is.na(por_aluno$atraso_M0)), "de", nrow(por_aluno), "\n\n")
print(tabela)
cat("\nAcerto por trajetória (%):\n");     print(round(100 * acerto_trajetoria, 1))
cat("\nAcerto por professor atual (%):\n"); print(round(100 * acerto_professor, 1))
cat("\nAcerto por histórico (%):\n"); print(round(100 * acerto_historico, 1))

write.csv(tabela, "results/tabela_modelos.csv", row.names = FALSE)
write.csv(round(100 * acerto_trajetoria, 1), "results/acerto_trajetoria.csv")
write.csv(round(100 * acerto_professor, 1), "results/acerto_professor.csv")
write.csv(round(100 * acerto_historico, 1), "results/acerto_historico.csv")

# -----------------------------------------------------------------------------
# 7. Convergência do M2 com os diagnósticos do coda (como no Módulo V/VI)
#    Gelman-Rubin (~1), Geweke (|z| < 1,96), Heidelberger-Welch (passou),
#    ESS (> 400), MCSE e intervalo HPD de 95%.
# -----------------------------------------------------------------------------

pars <- c("tau_esc", "tau_prof", "omega", "sigma", "sigma_hist", "gamma", "delta_a")

coda_hn <- diagnosticos_coda(fit_hn, pars)
coda_hc <- diagnosticos_coda(fit_hc, pars)
cat("\n===== Convergência (coda) — half-normal =====\n")
print(coda_hn$tabela, digits = 3, row.names = FALSE)
cat("\n===== Convergência (coda) — half-Cauchy =====\n")
print(coda_hc$tabela, digits = 3, row.names = FALSE)
cat("\nRaftery-Lewis:\n"); print(coda_hn$raftery)
write.csv(rbind(cbind(prior = "half-normal(0, 10)", coda_hn$tabela),
                cbind(prior = "half-Cauchy(0, 5)",  coda_hc$tabela)),
          "results/diagnostico_M2.csv", row.names = FALSE)

# Autocorrelação das cadeias (cadeia 1), como o autocorr.plot do Módulo VI
# (mesmos 4 parâmetros do traceplot da fig5, para os dois gráficos ficarem lado a lado no slide)
png("results/fig8_autocorrelacao.png", width = 2400, height = 1500, res = 250)
par(mfrow = c(2, 2), mar = c(4, 4.5, 2.5, 1), cex = 1.35)  # fonte grande para o slide
autocorr.plot(coda_hn$mc[[1]][, c("tau_esc", "tau_prof", "omega", "sigma")],
              auto.layout = FALSE, lag.max = 20, lwd = 3)
dev.off()

# -----------------------------------------------------------------------------
# 8. Recuperação dos efeitos (média e HPD 95%) — só existe na simulação
# -----------------------------------------------------------------------------

verd_esc <- tapply(verdade$gamma_origem, substr(verdade$prof_origem_real, 1, 2), function(x) x[1])
verd_po  <- tapply(verdade$delta_origem, verdade$prof_origem_real, function(x) x[1])
verd_pa  <- tapply(verdade$delta_atual, alunos$prof_atual[match(verdade$id, alunos$id)], function(x) x[1])

efeito <- function(par, nomes, verdadeiro, grupo) {
  x   <- as.matrix(fit_hn, pars = par)
  hpd <- HPDinterval(as.mcmc(x), prob = 0.95)
  data.frame(grupo = grupo, nome = nomes, real = as.numeric(verdadeiro[nomes]),
             media = colMeans(x), inf = hpd[, "lower"], sup = hpd[, "upper"])
}
recuperacao <- rbind(
  efeito("gamma",   d_stan$nomes_esc, verd_esc, "Escola de origem (γ)"),
  efeito("delta_a", d_stan$nomes_pa,  verd_pa,  "Professor atual (δ)"),
  efeito("delta_o", d_stan$nomes_po,  verd_po,  "Professor de origem (δ)")
)
write.csv(recuperacao, "results/recuperacao_efeitos.csv", row.names = FALSE)

# Censura x ingênuo: o que cada modelo conclui sobre as escolas e o histórico
resumo_par <- function(fit, par) {
  x <- as.matrix(fit, pars = par); hpd <- HPDinterval(as.mcmc(x), prob = 0.95)
  data.frame(media = colMeans(x), inf = hpd[, 1], sup = hpd[, 2])
}
censura_tab <- rbind(
  cbind(modelo = "com censura", nome = c(d_stan$nomes_esc, "sigma_hist"),
        real = c(as.numeric(verd_esc[d_stan$nomes_esc]), 5),
        rbind(resumo_par(fit_hn, "gamma"), resumo_par(fit_hn, "sigma_hist"))),
  cbind(modelo = "ingênuo", nome = c(d_stan$nomes_esc, "sigma_hist"),
        real = c(as.numeric(verd_esc[d_stan$nomes_esc]), 5),
        rbind(resumo_par(fit_ing, "gamma"), resumo_par(fit_ing, "sigma_hist"))))
cat("\n===== Escolas e erro do histórico: com censura x ingênuo =====\n")
print(censura_tab, digits = 3, row.names = FALSE)
write.csv(censura_tab, "results/censura_vs_ingenuo.csv", row.names = FALSE)

# Partial pooling / shrinkage (Tópico 9)
# Compara, para cada professor de origem:
#   "sem pooling": média, só entre os alunos dele, do que sobra do histórico
#                  depois de tirar theta_1 e gamma (cada professor isolado);
#   "partial pooling": estimativa do M2 (delta_o ~ N(0, tau_prof^2)).
# Teoria: partial = fator x sem pooling, com fator = tau^2 / (tau^2 + sigma_h^2 / n).
# Professor com poucos alunos tem pouca informação própria -> fator menor ->
# estimativa mais "puxada" para 0 (a média do grupo).
am_hn_sh <- amostras_m2(fit_hn, d_stan)
th1   <- colMeans(am_hn_sh$theta[, , 1])
gam   <- colMeans(am_hn_sh$gamma)[d_stan$escola]
sobra <- d_stan$nota_hist - th1 - gam
conh  <- !is.na(alunos$prof_origem) & alunos$nota_hist > 60   # só histórico observado
sem_pooling <- tapply(sobra[conh], alunos$prof_origem[conh], mean)
n_po  <- table(alunos$prof_origem[conh])
shrinkage <- data.frame(
  nome        = d_stan$nomes_po,
  n_alunos    = as.integer(n_po[d_stan$nomes_po]),
  real        = as.numeric(verd_po[d_stan$nomes_po]),
  sem_pooling = as.numeric(sem_pooling[d_stan$nomes_po]),
  partial     = colMeans(am_hn_sh$delta_o))
shrinkage$fator <- shrinkage$partial / shrinkage$sem_pooling
shrinkage$fator_teorico <- mean(am_hn_sh$tau_prof)^2 /
  (mean(am_hn_sh$tau_prof)^2 + mean(am_hn_sh$sigma_hist)^2 / shrinkage$n_alunos)
cat("\n===== Partial pooling dos professores de origem =====\n")
print(shrinkage[order(shrinkage$n_alunos), ], digits = 2, row.names = FALSE)
write.csv(shrinkage, "results/shrinkage_professores.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 9. Ajuste e capacidade preditiva (Tópico 7) — NÃO usa a verdade da simulação
#    M2 completo x M2 sem escola e professor (mesmos dados, mesmas priors)
# -----------------------------------------------------------------------------

am_hn <- amostras_m2(fit_hn, d_stan)
d_se  <- dados_stan(alunos, notas, com_efeitos = 0)
am_se <- amostras_m2(fit_se, d_se)

# 9.1 Verossimilhança: amostras x observações (1.080 notas do ano + 120 históricos)
ll_hn <- m2_loglik(am_hn, d_stan)
ll_se <- m2_loglik(am_se, d_se)

# 9.2 DIC (como no Módulo VI)
dic_hn <- dic(ll_hn, m2_loglik(media_amostras(am_hn), d_stan))
dic_se <- dic(ll_se, m2_loglik(media_amostras(am_se), d_se))

# 9.3 WAIC e LOO POR NOTA (pacote loo); r_eff corrige pela autocorrelação das cadeias
#     Primeira tentativa: cada nota tem seu próprio theta, então esconder uma nota
#     muda muito a posterior e a aproximação do LOO falha (Pareto-k > 0,7).
cadeia <- rep(1:4, each = 1000)
waic_hn <- waic(ll_hn);  waic_se <- waic(ll_se)
loo_hn  <- loo(ll_hn, r_eff = relative_eff(exp(ll_hn), cadeia))
loo_se  <- loo(ll_se, r_eff = relative_eff(exp(ll_se), cadeia))

# 9.4 LOO POR ALUNO, com os theta integrados pelo filtro de Kalman (funcoes.R, 8b)
#     Pergunta: "chegou um aluno novo; o modelo prevê bem as notas dele?"
lla_hn <- m2_loglik_aluno(am_hn, d_stan)       # amostras x 120 alunos

# Conferência: para um aluno censurado e uma amostra, a conta do Kalman +
# suavizador tem de bater com a integral numérica em theta_1
i_c <- which(d_stan$nota_hist <= 60)[1]; s1 <- 1
am1 <- lapply(am_hn, function(x) if (is.matrix(x)) x[s1, , drop = FALSE]
              else if (length(dim(x)) == 3) x[s1, , , drop = FALSE] else x[s1])
d1  <- within(d_stan, { y <- y[i_c, , drop = FALSE]; nota_hist <- nota_hist[i_c]
  escola <- escola[i_c]; prof_orig <- prof_orig[i_c]; prof_atual <- prof_atual[i_c]; N <- 1 })
desl <- am1$gamma[1, d1$escola] + if (d1$prof_orig > 0) am1$delta_o[1, d1$prof_orig] else 0
s2h  <- am1$sigma_hist^2 + if (d1$prof_orig > 0) 0 else am1$tau_prof^2
integrando <- function(th1) sapply(th1, function(t1) {
  m <- t1; v <- 0; l <- 0
  for (t in 1:d1$T) { if (t > 1) v <- v + am1$omega^2
    p <- passo_kalman(m, v, d1$y[1, t], am1$delta_a[1, d1$prof_atual], am1$sigma^2)
    l <- l + p$ll; m <- p$m; v <- p$v }
  exp(l) * dnorm(t1, 50, 25) * pnorm(60, t1 + desl, sqrt(s2h)) })
cat("Conferência do LOO censurado (Kalman x integral):",
    round(m2_loglik_aluno(am1, d1)[1, 1], 6),
    round(log(integrate(integrando, -50, 150, subdivisions = 5000, rel.tol = 1e-12)$value), 6), "\n")
lla_se <- m2_loglik_aluno(am_se, d_se)
looa_hn <- loo(lla_hn, r_eff = relative_eff(exp(lla_hn), cadeia))
looa_se <- loo(lla_se, r_eff = relative_eff(exp(lla_se), cadeia))

comparacao <- data.frame(
  modelo = c("M2 completo", "M2 sem escola e professor"),
  DIC    = round(c(dic_hn["DIC"], dic_se["DIC"]), 1),
  pD     = round(c(dic_hn["pD"],  dic_se["pD"]), 1),
  WAIC_nota  = round(c(waic_hn$estimates["waic", "Estimate"], waic_se$estimates["waic", "Estimate"]), 1),
  LOOIC_nota = round(c(loo_hn$estimates["looic", "Estimate"], loo_se$estimates["looic", "Estimate"]), 1),
  k_ruins_nota  = c(sum(loo_hn$diagnostics$pareto_k > 0.7), sum(loo_se$diagnostics$pareto_k > 0.7)),
  LOOIC_aluno   = round(c(looa_hn$estimates["looic", "Estimate"], looa_se$estimates["looic", "Estimate"]), 1),
  k_ruins_aluno = c(sum(looa_hn$diagnostics$pareto_k > 0.7), sum(looa_se$diagnostics$pareto_k > 0.7)),
  k_max_aluno   = round(c(max(looa_hn$diagnostics$pareto_k), max(looa_se$diagnostics$pareto_k)), 2)
)
cat("\n===== Comparação de modelos (menor = melhor) =====\n")
print(comparacao, row.names = FALSE)
cat("\nLOO por nota (1.200 observações):\n")
print(loo_compare(list(completo = loo_hn, sem_efeitos = loo_se)))
cat("\nLOO por aluno (120 alunos, theta integrados):\n")
print(loo_compare(list(completo = looa_hn, sem_efeitos = looa_se)))

# Robustez: os poucos alunos com Pareto-k > 0,7 são os de professor de origem
# que só aparece neles (o delta dele depende só daquele aluno). Mudam a conclusão?
dif_aluno <- looa_hn$pointwise[, "elpd_loo"] - looa_se$pointwise[, "elpd_loo"]
k_ruim    <- looa_hn$diagnostics$pareto_k > 0.7
cat("Alunos com Pareto-k > 0,7:", alunos$id[k_ruim], "\n")
cat("Diferença de ELPD sem esses alunos:", round(sum(dif_aluno[!k_ruim]), 1), "\n")
cat("Alunos em que o modelo completo prevê melhor:", sum(dif_aluno > 0), "de", length(dif_aluno), "\n")
write.csv(comparacao, "results/comparacao_modelos.csv", row.names = FALSE)

# 9.5 Verificação preditiva a posteriori, em três níveis (Tópico 9):
#     global (todas as notas), grupo (professor e escola) e individual (cada aluno)
set.seed(2026)
y_obs  <- as.vector(t(d_stan$y))
h_obs  <- d_stan$nota_hist
rep_hn <- simular_rep(am_hn, d_stan)
yrep_hn <- rep_hn$y
grupo_prof   <- rep(alunos$prof_atual, each = d_stan$T)
grupo_escola <- alunos$escola_origem

# Nível individual: o quanto as 9 notas de cada aluno variam (DP por aluno)
dp_aluno <- function(y) apply(matrix(y, ncol = d_stan$T, byrow = TRUE), 1, sd)
dp_obs   <- dp_aluno(y_obs)
dp_rep   <- t(apply(yrep_hn, 1, dp_aluno))       # simulações x 120 alunos

# 9.6 p-valor bayesiano: P(estatística simulada >= observada); ~0,5 é bom, perto de 0 ou 1 é ruim
p_bayes <- function(sim, obs) mean(sim >= obs)
ppc_p <- data.frame(
  nivel = c(rep("global", 4), "individual", "grupo (professor)", "grupo (escola)"),
  estatistica = c("média", "desvio-padrão", "mínimo", "máximo",
                  "média dos DPs por aluno",
                  "maior diferença entre médias dos professores",
                  "maior diferença entre médias do histórico por escola"),
  p_bayesiano = round(c(
    p_bayes(apply(yrep_hn, 1, mean), mean(y_obs)),
    p_bayes(apply(yrep_hn, 1, sd),   sd(y_obs)),
    p_bayes(apply(yrep_hn, 1, min),  min(y_obs)),
    p_bayes(apply(yrep_hn, 1, max),  max(y_obs)),
    p_bayes(rowMeans(dp_rep), mean(dp_obs)),
    p_bayes(apply(yrep_hn, 1, function(y) diff(range(tapply(y, grupo_prof, mean)))),
            diff(range(tapply(y_obs, grupo_prof, mean)))),
    p_bayes(apply(rep_hn$h, 1, function(h) diff(range(tapply(h, grupo_escola, mean)))),
            diff(range(tapply(h_obs, grupo_escola, mean))))), 3))
cat("\n===== Verificação preditiva: p-valores bayesianos =====\n")
print(ppc_p, row.names = FALSE)
write.csv(ppc_p, "results/ppc_pvalores.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 10. Teste de hipóteses bayesiano (Tópico 8)
#     Hipóteses = regiões do parâmetro; a resposta é a probabilidade a posteriori.
#     ROPE de ±2 pontos: diferença considerada desprezível na prática.
# -----------------------------------------------------------------------------

hipoteses <- rbind(
  do.call(rbind, lapply(seq_along(d_stan$nomes_pa), function(j)
    teste_hipotese(am_hn$delta_a[, j], paste0("delta_", d_stan$nomes_pa[j])))),
  do.call(rbind, lapply(seq_along(d_stan$nomes_esc), function(j)
    teste_hipotese(am_hn$gamma[, j], paste0("gamma_", d_stan$nomes_esc[j]))))
)
cat("\n===== Teste de hipóteses bayesiano =====\n")
print(hipoteses, row.names = FALSE)
write.csv(hipoteses, "results/teste_hipoteses.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 11. Figuras
# -----------------------------------------------------------------------------

cores  <- c("Real" = "black", "M0" = "#E69F00", "M1" = "#56B4E9", "M2" = "#0072B2")
faixas <- data.frame(corte = CORTES[2:4])

# 11.1 Um aluno de entrada média de cada trajetória
exemplos <- sapply(c("evolui", "involui", "estagna", "gangorra"), function(tr) {
  res$id[res$entrada == "medio" & res$trajetoria == tr][1]
})
d1 <- res[res$id %in% exemplos, ]
d1$painel <- paste0(d1$trajetoria, " (", d1$id, ", prof. ", d1$prof_atual, ")")

fig1 <- ggplot(d1, aes(x = aval)) +
  geom_hline(data = faixas, aes(yintercept = corte), linetype = "dashed", colour = "grey60") +
  geom_ribbon(aes(ymin = pmax(M2_media - 1.96 * M2_s, 0), ymax = pmin(M2_media + 1.96 * M2_s, 100)),
              fill = "#0072B2", alpha = 0.15) +
  geom_point(aes(y = nota), colour = "grey40") +
  geom_line(aes(y = theta_real, colour = "Real"), linewidth = 1) +
  geom_line(aes(y = m0_nota, colour = "M0")) +
  geom_line(aes(y = M1_media, colour = "M1")) +
  geom_line(aes(y = M2_media, colour = "M2"), linewidth = 0.8) +
  facet_wrap(~ painel) +
  scale_colour_manual(values = cores, name = NULL) +
  scale_x_continuous(breaks = 1:9) +
  labs(x = "Avaliação", y = "Proficiência (0-100)",
       title = "Proficiência real e estimada por M0, M1 e M2 (faixa = IC 95% do M2)",
       subtitle = "Pontos = notas observadas; linhas tracejadas = cortes 40 / 60 / 80") +
  theme_bw()
ggsave("results/fig1_trajetorias.png", fig1, width = 9, height = 6, dpi = 300)

# 11.2 Probabilidade de cada nível ao longo do ano (aluno gangorra, M1 x M2)
id_g <- exemplos["gangorra"]
empilhar <- function(post, nome) {
  d <- post[post$id == id_g, c("aval", NIVEIS)]
  data.frame(modelo = nome, aval = rep(d$aval, 4),
             nivel = factor(rep(NIVEIS, each = nrow(d)), levels = rev(NIVEIS)),
             prob = unlist(d[, NIVEIS]))
}
d2 <- rbind(empilhar(post_M1, "M1"), empilhar(post_M2, "M2"))

fig2 <- ggplot(d2, aes(x = aval, y = prob, fill = nivel)) +
  geom_area() +
  geom_hline(yintercept = limiar, linetype = "dashed") +
  facet_wrap(~ modelo) +
  scale_x_continuous(breaks = 1:9) +
  scale_fill_manual(values = c("Avançado" = "#1A9850", "Proficiente" = "#91CF60",
                               "Básico" = "#FC8D59", "Abaixo do Básico" = "#D73027"),
                    name = "Nível") +
  labs(x = "Avaliação", y = "Probabilidade posterior",
       title = paste("Probabilidade de cada nível ao longo do ano — aluno", id_g, "(gangorra)"),
       subtitle = "Linha tracejada = limiar de 0,7 para trocar de nível") +
  theme_bw()
ggsave("results/fig2_probabilidades.png", fig2, width = 9, height = 5, dpi = 300)

# 11.3 Taxa de acerto por avaliação (regra do limiar + M0)
d3 <- acerto_aval[acerto_aval$modelo %in% c("M0", "M1_limiar", "M2_limiar"), ]
fig3 <- ggplot(d3, aes(x = aval, y = 100 * acerto, colour = modelo)) +
  geom_line(linewidth = 1) + geom_point() +
  scale_x_continuous(breaks = 1:9) +
  scale_colour_manual(values = c("M0" = "#E69F00", "M1_limiar" = "#56B4E9",
                                 "M2_limiar" = "#0072B2"), name = "Modelo") +
  labs(x = "Avaliação", y = "Alunos classificados no nível certo (%)",
       title = "Taxa de acerto da classificação ao longo do ano") +
  theme_bw()
ggsave("results/fig3_acerto.png", fig3, width = 8, height = 5, dpi = 300)

# 11.4 Efeitos verdadeiros x estimados (M2)
d4 <- recuperacao[recuperacao$grupo != "Professor de origem (δ)", ]
fig4 <- ggplot(d4, aes(y = nome)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = inf, xmax = sup), height = 0.2, colour = "#0072B2") +
  geom_point(aes(x = media, shape = "Estimado (M2)"), colour = "#0072B2", size = 2.5) +
  geom_point(aes(x = real, shape = "Verdadeiro"), colour = "black", size = 2.5) +
  facet_wrap(~ grupo, scales = "free_y") +
  scale_shape_manual(values = c("Estimado (M2)" = 16, "Verdadeiro" = 4), name = NULL) +
  labs(x = "Efeito na nota (pontos)", y = NULL,
       title = "Recuperação dos efeitos de escola e professor (média e HPD 95%)") +
  theme_bw()
ggsave("results/fig4_efeitos.png", fig4, width = 9, height = 4, dpi = 300)

# 11.5 Diagnóstico visual: cadeias dos desvios-padrão (half-normal)
fig5 <- stan_trace(fit_hn, pars = c("tau_esc", "tau_prof", "omega", "sigma")) +
  theme(text = element_text(size = 17), legend.position = "bottom")  # fonte grande para o slide
ggsave("results/fig5_cadeias.png", fig5, width = 9, height = 5.5, dpi = 300)

# 11.6 Sensibilidade: posterior de tau com half-normal x half-Cauchy
d6 <- rbind(
  data.frame(prior = "half-normal(0, 10)", par = rep(c("tau_esc", "tau_prof"), each = 4000),
             valor = c(extract(fit_hn)$tau_esc, extract(fit_hn)$tau_prof)),
  data.frame(prior = "half-Cauchy(0, 5)", par = rep(c("tau_esc", "tau_prof"), each = 4000),
             valor = c(extract(fit_hc)$tau_esc, extract(fit_hc)$tau_prof)))
fig6 <- ggplot(d6, aes(x = valor, fill = prior)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~ par, scales = "free") +
  coord_cartesian(xlim = c(0, 40)) +
  labs(x = "Desvio-padrão do efeito (pontos)", y = "Densidade posterior", fill = "Prior",
       title = "Sensibilidade à prior dos desvios-padrão") +
  theme_bw()
ggsave("results/fig6_sensibilidade.png", fig6, width = 9, height = 4, dpi = 300)

# 11.7 Verificação preditiva a posteriori nos três níveis
# Azul mais escuro que o padrão do bayesplot (as réplicas y_rep ficavam quase invisíveis no projetor)
color_scheme_set(c("#8fb3d9", "#6a9bd1", "#3f7fbf", "#2c6aa8", "#1b4f8a", "#0b2545"))
fig7a <- ppc_dens_overlay(y_obs, yrep_hn[1:50, ]) +
  labs(title = "Global: notas reais (y) x simuladas (y_rep)", x = "Nota")
fig7b <- ppc_stat_grouped(y_obs, yrep_hn, group = grupo_prof, stat = "mean", binwidth = 0.2) +
  labs(title = "Grupo: média das notas por professor", x = "Média")
fig7c <- ppc_stat_grouped(h_obs, rep_hn$h, group = grupo_escola, stat = "mean", binwidth = 0.5) +
  labs(title = "Grupo: média do histórico por escola de origem", x = "Média")
fig7d <- ppc_dens_overlay(dp_obs, dp_rep[1:50, ]) +
  labs(title = "Individual: DP das 9 notas de cada aluno", x = "DP")
ggsave("results/fig7a_ppc_densidade.png", fig7a, width = 8, height = 4.5, dpi = 300, bg = "white")
ggsave("results/fig7b_ppc_professor.png", fig7b, width = 9, height = 3.5, dpi = 300, bg = "white")
ggsave("results/fig7c_ppc_escola.png", fig7c, width = 11, height = 3.5, dpi = 300, bg = "white")
ggsave("results/fig7d_ppc_individual.png", fig7d, width = 8, height = 4.5, dpi = 300, bg = "white")
ggsave("results/fig7_ppc_niveis.png",
       bayesplot_grid(fig7a, fig7d, fig7b, fig7c, grid_args = list(ncol = 2)),
       width = 14, height = 8, dpi = 250, bg = "white")

# 11.8 Checagem preditiva a priori: notas geradas só pelas priors x faixa 0-100
d9 <- data.frame(nota = as.vector(y_prior[1:200, ]))
fig9 <- ggplot(d9, aes(x = nota)) +
  geom_histogram(binwidth = 5, fill = "grey70", colour = "white") +
  geom_vline(xintercept = c(0, 100), linetype = "dashed", colour = "#D73027") +
  labs(x = "Nota simulada só a partir das priors", y = "Frequência",
       title = "Checagem preditiva a priori",
       subtitle = paste0("Linhas tracejadas = limites 0 e 100; fora da faixa: ",
                         format(round(100 * mean(y_prior < 0 | y_prior > 100), 1), decimal.mark = ","),
                         "% das notas")) +
  theme_bw()
ggsave("results/fig9_prior_preditiva.png", fig9, width = 8, height = 4.5, dpi = 300)

# 11.9 Partial pooling: estimativa sem pooling x M2, tamanho = nº de alunos
fig10 <- ggplot(shrinkage[!is.na(shrinkage$sem_pooling), ], aes(x = sem_pooling, y = partial)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
  geom_segment(aes(xend = sem_pooling, yend = sem_pooling), colour = "#0072B2", alpha = 0.5) +
  geom_point(aes(size = n_alunos), colour = "#0072B2") +
  scale_size_continuous(range = c(2, 7), breaks = c(1, 3, 5, 7), name = "Alunos") +
  labs(x = "Sem pooling (cada professor só com os seus alunos)",
       y = "Partial pooling (M2)",
       title = "Partial pooling: professores com poucos alunos são puxados para 0",
       subtitle = "Tracejado: sem encolhimento. Fator teórico = 0,68 com 1 aluno; 0,91 com 5 alunos") +
  theme_bw()
ggsave("results/fig10_shrinkage.png", fig10, width = 8, height = 5.5, dpi = 300)

# 11.10 Desenho da simulação (figura para a apresentação)
#   escolas de origem (E5 em vermelho: infla +12), 12 perfis de alunos
#   (3 entradas x 4 trajetórias) e professores atuais como carinhas
circulo <- function(cx, cy, r, id) {
  a <- seq(0, 2 * pi, length.out = 60)
  data.frame(x = cx + r * cos(a), y = cy + r * sin(a), id = id)
}
boca <- function(cx, cy, curva, id) {       # curva > 0 sorriso; < 0 triste; 0 reta
  x <- seq(-0.45, 0.45, length.out = 20)
  data.frame(x = cx + x, y = cy - 0.35 + curva * (x^2 - 0.1), id = id)
}
casa <- function(cx, cy, id) {               # quadrado + telhado
  data.frame(x = cx + c(-0.5, 0.5, 0.5, 0.65, 0, -0.65, -0.5, -0.5),
             y = cy + c(-0.5, -0.5, 0.3, 0.3, 0.9, 0.3, 0.3, -0.5), id = id)
}
# escolas
esc <- data.frame(nome = paste0("E", 1:5), x = 1, y = c(9, 7.4, 5.8, 4.2, 2.6),
                  efeito = c("0", "−4", "+3", "0", "+12 (infla)"))
casas <- do.call(rbind, lapply(1:5, function(j) casa(esc$x[j], esc$y[j], esc$nome[j])))
casas$cor <- ifelse(casas$id == "E5", "infla", "normal")
# perfis: 3 entradas (linhas) x 4 trajetórias (colunas)
traj_forma <- list(evolui = function(t) t, involui = function(t) -t,
                   estagna = function(t) 0 * t, gangorra = function(t) sin(4 * pi * t))
ent <- c(forte = 8, medio = 5.5, fraco = 3)
mini <- do.call(rbind, lapply(names(traj_forma), function(tr) {
  do.call(rbind, lapply(names(ent), function(e) {
    t <- seq(0, 1, length.out = 40)
    cx <- 6.3 + 1.6 * (match(tr, names(traj_forma)) - 1)
    data.frame(x = cx - 0.6 + 1.2 * t, y = ent[e] + 0.6 * traj_forma[[tr]](t),
               id = paste(tr, e))
  }))
}))
# professores
prof <- data.frame(nome = c("PF", "PM", "PR"), x = 14.4, y = c(8, 5.5, 3),
                   curva = c(2.2, 0, -2.2), efeito = c("facilita (+8)", "mediano (0)", "rigoroso (−10)"))
rostos <- do.call(rbind, lapply(1:3, function(j) circulo(prof$x[j], prof$y[j], 0.8, prof$nome[j])))
bocas  <- do.call(rbind, lapply(1:3, function(j) boca(prof$x[j], prof$y[j], prof$curva[j], prof$nome[j])))
olhos  <- data.frame(x = rep(prof$x, each = 2) + c(-0.28, 0.28), y = rep(prof$y, each = 2) + 0.22)

fig0 <- ggplot() +
  geom_polygon(data = casas, aes(x, y, group = id, fill = cor), colour = NA) +
  scale_fill_manual(values = c(normal = "grey15", infla = "#D73027"), guide = "none") +
  geom_text(data = esc, aes(x + 0.95, y, label = paste0(nome, "  ", efeito)), hjust = 0, size = 3.6,
            colour = ifelse(esc$nome == "E5", "#D73027", "grey15")) +
  geom_path(data = mini, aes(x, y, group = id), colour = "#0072B2", linewidth = 1.1) +
  annotate("text", x = 6.3 + 1.6 * (0:3), y = 10, label = names(traj_forma), size = 3.8, fontface = "bold") +
  annotate("text", x = 5.3, y = ent, label = c("forte\n~85", "médio\n~55", "fraco\n~35"), hjust = 1, size = 3.6) +
  geom_polygon(data = rostos, aes(x, y, group = id), fill = "#FFE08A", colour = "grey20") +
  geom_point(data = olhos, aes(x, y), size = 1.8) +
  geom_path(data = bocas, aes(x, y, group = id), linewidth = 1) +
  geom_text(data = prof, aes(x + 1.05, y, label = paste0(nome, "\n", efeito)), hjust = 0, size = 3.6) +
  annotate("text", x = c(1.8, 8.7, 15.2), y = 11.3, size = 4.6, fontface = "bold",
           label = c("Escola de origem", "120 alunos: 12 perfis × 10", "Professor atual")) +
  annotate("text", x = 8.7, y = 1.3, size = 3.6, colour = "grey30",
           label = "nota do histórico (mínimo 60: aprovação)  →  9 avaliações no ano (3 trimestres: 2 provas + 1 trabalho)") +
  coord_equal(xlim = c(0.2, 17.6), ylim = c(0.8, 11.6)) +
  theme_void()
ggsave("results/fig0_desenho.png", fig0, width = 12, height = 7, dpi = 300, bg = "white")

# 11.11 Probabilidade de cada nível por perfil (M2), média dos 10 alunos do perfil
#   (versão anterior da fig12; não está na apresentação)
pp <- merge(post_M2[, c("id", "aval", NIVEIS)],
            unique(res[, c("id", "entrada", "trajetoria")]), by = "id")
pp <- aggregate(pp[, NIVEIS], by = list(entrada = pp$entrada, trajetoria = pp$trajetoria,
                                        aval = pp$aval), FUN = mean)
d11 <- data.frame(entrada = rep(pp$entrada, 4), trajetoria = rep(pp$trajetoria, 4),
                  aval = rep(pp$aval, 4),
                  nivel = factor(rep(NIVEIS, each = nrow(pp)), levels = rev(NIVEIS)),
                  prob = unlist(pp[, NIVEIS]))
d11$entrada    <- factor(d11$entrada, levels = c("forte", "medio", "fraco"),
                         labels = c("Entrada forte (~85)", "Entrada média (~55)", "Entrada fraca (~35)"))
d11$trajetoria <- factor(d11$trajetoria, levels = c("evolui", "involui", "estagna", "gangorra"))
fig11 <- ggplot(d11, aes(x = aval, y = prob, fill = nivel)) +
  geom_area(colour = "white", linewidth = 0.1) +
  facet_grid(entrada ~ trajetoria) +
  scale_x_continuous(breaks = c(2, 5, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), labels = c("0,25", "0,5", "0,75"), expand = c(0, 0)) +
  scale_fill_manual(values = c("Avançado" = "#1A9850", "Proficiente" = "#91CF60",
                               "Básico" = "#FC8D59", "Abaixo do Básico" = "#D73027"), name = "Nível") +
  labs(x = "Avaliação", y = "Probabilidade posterior (média dos 10 alunos)") +
  theme_bw(base_size = 13) + theme(panel.spacing = unit(0.5, "lines"), strip.text = element_text(face = "bold"))
ggsave("results/fig11_prob_perfis.png", fig11, width = 11, height = 6.5, dpi = 300)

# 11.12 Probabilidade de cada nível para UM aluno de cada perfil (M2)
#   faixa de cima = nível REAL do aluno (da simulação), para comparar
alunos_ex <- sapply(split(alunos$id, paste(res$entrada[match(alunos$id, res$id)],
                                           res$trajetoria[match(alunos$id, res$id)])), `[`, 1)
pa <- post_M2[post_M2$id %in% alunos_ex, c("id", "aval", NIVEIS)]
d12 <- data.frame(id = rep(pa$id, 4), aval = rep(pa$aval, 4),
                  nivel = factor(rep(NIVEIS, each = nrow(pa)), levels = rev(NIVEIS)),
                  prob = unlist(pa[, NIVEIS]))
info <- unique(res[res$id %in% alunos_ex, c("id", "entrada", "trajetoria", "prof_atual")])
d12  <- merge(d12, info, by = "id")
real <- res[res$id %in% alunos_ex, c("id", "aval", "nivel_real")]
real <- merge(real, info, by = "id")
real$nivel <- factor(real$nivel_real, levels = rev(NIVEIS))
for (d in c("d12", "real")) {
  x <- get(d)
  x$entrada    <- factor(x$entrada, levels = c("forte", "medio", "fraco"),
                         labels = c("Entrada forte (~85)", "Entrada média (~55)", "Entrada fraca (~35)"))
  x$trajetoria <- factor(x$trajetoria, levels = c("evolui", "involui", "estagna", "gangorra"))
  assign(d, x)
}
rotulos <- unique(d12[, c("entrada", "trajetoria", "id", "prof_atual")])
fig12 <- ggplot(d12, aes(x = aval, y = prob, fill = nivel)) +
  geom_area(colour = "white", linewidth = 0.1) +
  geom_tile(data = real, aes(x = aval, y = 1.09, fill = nivel), height = 0.12, width = 1,
            inherit.aes = FALSE) +
  geom_text(data = rotulos, aes(x = 9, y = -0.06, label = paste0(id, " · prof. ", prof_atual)),
            inherit.aes = FALSE, hjust = 1, size = 3, colour = "grey30") +
  facet_grid(entrada ~ trajetoria) +
  scale_x_continuous(breaks = c(2, 5, 8), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75), labels = c("0,25", "0,5", "0,75"),
                     expand = c(0, 0), limits = c(-0.12, 1.16)) +
  scale_fill_manual(values = c("Avançado" = "#1A9850", "Proficiente" = "#91CF60",
                               "Básico" = "#FC8D59", "Abaixo do Básico" = "#D73027"), name = "Nível") +
  labs(x = "Avaliação", y = "Probabilidade posterior (M2)",
       subtitle = "Faixa de cima = nível REAL do aluno; área = o que o modelo diz") +
  theme_bw(base_size = 13) + theme(panel.spacing = unit(0.5, "lines"),
                                   strip.text = element_text(face = "bold"))
ggsave("results/fig12_prob_alunos.png", fig12, width = 11, height = 7, dpi = 300)

# 11.13 Ciclo bayesiano do modelo (figura para a apresentação)
caixa <- function(x, y, txt, cor = "grey20", fundo = "white", tam = 5.6)
  annotate("label", x = x, y = y, label = txt, colour = cor, fill = fundo, size = tam,
           label.padding = unit(0.6, "lines"), lineheight = 0.95, fontface = "bold")
seta <- function(x, y, xend, yend, curv = 0)
  annotate("curve", x = x, y = y, xend = xend, yend = yend, curvature = curv,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"), colour = "grey35", linewidth = 0.8)
fig13 <- ggplot() + xlim(-0.3, 13.6) + ylim(-0.2, 6.3) +
  caixa(2, 5.4, "HISTÓRICO\nnota que o aluno traz") +
  seta(2, 4.85, 2, 3.85) +
  annotate("text", x = 2.25, y = 4.35, hjust = 0, size = 5, colour = "#D73027", fontface = "bold",
           label = "M2: desconta a\nescola de origem") +
  caixa(2, 3.2, "PRIOR\nonde acho que\no aluno está", fundo = "#DCEBF7") +
  annotate("text", x = 4.3, y = 3.2, label = "×", size = 12) +
  caixa(6.5, 3.2, "PROVA t\nverossimilhança", fundo = "#FFF2CC") +
  annotate("text", x = 6.5, y = 4.45, size = 5, colour = "#D73027", fontface = "bold",
           label = "M2: desconta o professor") +
  annotate("text", x = 8.7, y = 3.2, label = "=", size = 12) +
  caixa(11, 3.2, "POSTERIORI\nonde o aluno está agora\n→ P(cada nível)", fundo = "#D9EAD3") +
  seta(11, 2.45, 8.3, 0.75, curv = -0.3) +
  caixa(6.5, 0.6, "+ incerteza: o aluno pode ter mudado\n→ vira a PRIOR da próxima prova", tam = 5.2) +
  seta(4.7, 0.75, 2, 2.4, curv = -0.3) +
  theme_void()
ggsave("results/fig13_ciclo.png", fig13, width = 11.5, height = 5.4, dpi = 300, bg = "white")

# 11.14 Censura x ingênuo: efeito de cada escola de origem
d14 <- censura_tab[censura_tab$nome != "sigma_hist", ]
d14$modelo <- factor(d14$modelo, levels = c("ingênuo", "com censura"),
                     labels = c("Ingênuo (60 = nota exata)", "Com censura (60 = \"60 ou menos\")"))
fig14 <- ggplot(d14, aes(y = nome)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = inf, xmax = sup, colour = modelo), height = 0.25,
                 position = position_dodge(width = 0.5)) +
  geom_point(aes(x = media, colour = modelo), size = 2.8, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = real), shape = 4, size = 4, stroke = 1.3) +
  scale_colour_manual(values = c("#D73027", "#0072B2"), name = NULL) +
  labs(x = "Efeito da escola de origem no histórico (pontos)", y = NULL,
       subtitle = "X = valor verdadeiro; ponto e barra = média e HPD 95%") +
  theme_bw(base_size = 13) + theme(legend.position = "bottom")
ggsave("results/fig14_censura.png", fig14, width = 9, height = 5, dpi = 300)

cat("\nFiguras e tabelas salvas em results/\n")
