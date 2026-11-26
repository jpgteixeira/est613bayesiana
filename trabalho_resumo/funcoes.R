# =============================================================================
# funcoes.R — Funções do estudo (carregado pela main.R com source())
# Classificação bayesiana de proficiência em Biologia (níveis tipo PAEBES)
# JPGTEIXEIRA — EST613 Estatística Bayesiana / UFV
#
# Conteúdo:
#   1. Níveis e probabilidades (normal truncada em [0, 100])
#   2. M0 — média simples
#   3. M1 — Normal-Normal dinâmico (filtro de Kalman simples)
#   4. Regras de classificação (maior probabilidade x limiar 0,7)
#   5. Métricas de comparação
#   6. M2 — modelo hierárquico em Stan
#   7-8. Verossimilhança do M2, DIC, simulação preditiva e checagem a priori
#   9. Diagnósticos do coda
#   10. Teste de hipóteses bayesiano
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Níveis e probabilidades
# -----------------------------------------------------------------------------

# Pontos de corte na escala 0-100 e nomes dos níveis
CORTES  <- c(0, 40, 60, 80, 100)
NIVEIS  <- c("Abaixo do Básico", "Básico", "Proficiente", "Avançado")

# Converte uma nota (ou proficiência) no nível correspondente
nota_para_nivel <- function(x) {
  cut(x, breaks = CORTES, labels = NIVEIS, include.lowest = TRUE)
}

# Probabilidade de o aluno estar em cada nível.
# A posterior é uma normal(m, s) TRUNCADA em [0, 100]:
# P(a < theta < b) = [Phi(b) - Phi(a)] / [Phi(100) - Phi(0)]
prob_niveis <- function(m, s) {
  acumulada <- pnorm(CORTES, mean = m, sd = s)
  p <- diff(acumulada) / (acumulada[5] - acumulada[1])
  names(p) <- NIVEIS
  p
}

# Média da normal truncada em [0, 100] (é a "nota estimada" do aluno)
media_truncada <- function(m, s) {
  a <- (0 - m) / s
  b <- (100 - m) / s
  m + s * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
}

# -----------------------------------------------------------------------------
# 2. M0 — Média simples (o que a escola faz hoje)
# -----------------------------------------------------------------------------

# Na avaliação t, a nota do aluno é a média das avaliações 1..t
m0_media_simples <- function(y) {
  cumsum(y) / seq_along(y)
}

# -----------------------------------------------------------------------------
# 3. M1 — Normal-Normal dinâmico (filtro de Kalman simples)
# -----------------------------------------------------------------------------
# Para UM aluno:
#   prior:     theta_0 ~ N(nota_hist, sd_prior^2)
#   dinâmica:  theta_t = theta_{t-1} + ruído,  ruído ~ N(0, sd_ruido^2)
#   dados:     y_t     = theta_t + erro,       erro  ~ N(0, sd_prova^2)
#
# Em cada avaliação são dois passos:
#   (a) PREVISÃO:    a incerteza cresce (v = v + sd_ruido^2), pois o aluno
#                    pode ter mudado desde a última prova;
#   (b) ATUALIZAÇÃO: fórmula fechada Normal-Normal (média ponderada
#                    pelas precisões, precisão = 1 / variância).
# A truncagem em [0, 100] só entra na hora de calcular as probabilidades
# (prob_niveis) e a média (media_truncada).

m1_kalman <- function(nota_hist, y, sd_prior, sd_prova, sd_ruido) {
  n_t <- length(y)
  m <- numeric(n_t)   # média da posterior em cada avaliação
  v <- numeric(n_t)   # variância da posterior em cada avaliação

  m_atual <- nota_hist
  v_atual <- sd_prior^2

  for (t in 1:n_t) {
    # (a) previsão: soma o ruído da dinâmica
    v_atual <- v_atual + sd_ruido^2

    # (b) atualização Normal-Normal
    precisao_prior <- 1 / v_atual
    precisao_dado  <- 1 / sd_prova^2
    v_atual <- 1 / (precisao_prior + precisao_dado)
    m_atual <- v_atual * (precisao_prior * m_atual + precisao_dado * y[t])

    m[t] <- m_atual
    v[t] <- v_atual
  }
  data.frame(aval = 1:n_t, m = m, s = sqrt(v))
}

# Roda o M1 para todos os alunos e devolve uma tabela longa com
# média, desvio e probabilidade de cada nível por aluno x avaliação
m1_todos <- function(alunos, notas, sd_prior, sd_prova, sd_ruido) {
  saida <- list()
  for (i in seq_len(nrow(alunos))) {
    id_i <- alunos$id[i]
    y_i  <- notas$nota[notas$id == id_i]          # já ordenado por aval
    post <- m1_kalman(alunos$nota_hist[i], y_i, sd_prior, sd_prova, sd_ruido)

    probs <- t(mapply(prob_niveis, post$m, post$s))
    post$media <- mapply(media_truncada, post$m, post$s)
    saida[[i]] <- data.frame(id = id_i, post, probs, check.names = FALSE)
  }
  do.call(rbind, saida)
}

# -----------------------------------------------------------------------------
# 4. Regras de classificação
# -----------------------------------------------------------------------------

# Regra A: nível com maior probabilidade em cada avaliação
regra_maior_prob <- function(probs) {
  NIVEIS[apply(probs, 1, which.max)]
}

# Regra B: só troca de nível quando a probabilidade do novo nível passa
# do limiar (0,7). Abaixo disso, o aluno continua no nível anterior.
# Evita que o aluno "pisque" entre níveis por causa de uma prova ruim.
regra_limiar <- function(probs, nivel_inicial, limiar = 0.7) {
  atual <- nivel_inicial
  saida <- character(nrow(probs))
  for (t in 1:nrow(probs)) {
    melhor <- which.max(probs[t, ])
    if (probs[t, melhor] > limiar) atual <- NIVEIS[melhor]
    saida[t] <- atual
  }
  saida
}

# Aplica as duas regras a uma posterior (M1 ou M2) e devolve as duas colunas.
# O nível de partida da regra do limiar é o nível da nota do histórico.
classificar <- function(post, alunos, limiar = 0.7) {
  post <- post[order(post$id, post$aval), ]
  probs <- as.matrix(post[, NIVEIS])
  maior <- regra_maior_prob(probs)
  lim <- character(nrow(post))
  for (i in seq_len(nrow(alunos))) {
    linhas  <- which(post$id == alunos$id[i])
    inicial <- as.character(nota_para_nivel(alunos$nota_hist[i]))
    lim[linhas] <- regra_limiar(probs[linhas, , drop = FALSE], inicial, limiar)
  }
  data.frame(id = post$id, aval = post$aval, maior_prob = maior, limiar = lim)
}

# -----------------------------------------------------------------------------
# 5. Métricas de comparação
# -----------------------------------------------------------------------------

# Número de trocas de nível numa sequência (ex.: B, B, P, B -> 2 trocas)
n_trocas <- function(niveis) {
  sum(niveis[-1] != niveis[-length(niveis)])
}

# Atraso na detecção da mudança de nível (em número de avaliações).
# t_mud = primeira avaliação em que o nível REAL difere do nível real da aval. 1.
# Atraso = quantas avaliações depois de t_mud o modelo passa a mostrar o
# novo nível real. NA = o aluno não mudou de nível; Inf = o modelo nunca detectou.
atraso_deteccao <- function(nivel_real, nivel_modelo) {
  t_mud <- which(nivel_real != nivel_real[1])[1]
  if (is.na(t_mud)) return(NA)
  acertos <- which(nivel_modelo == nivel_real[t_mud] & seq_along(nivel_modelo) >= t_mud)
  if (length(acertos) == 0) return(Inf)
  acertos[1] - t_mud
}

# -----------------------------------------------------------------------------
# 6. M2 — Modelo hierárquico dinâmico em Stan
# -----------------------------------------------------------------------------
# Mesma ideia do M1, mas agora com:
#   - gamma (escola de origem) e delta (professores) ESTIMADOS a partir dos dados;
#   - variâncias (sigma, omega, ...) também estimadas;
#   - um único tau_prof para todos os professores (origem e atuais), porque
#     com só 3 professores atuais seria quase impossível estimar o tau deles.
# O modelo tem duas "chaves" (opções) para refazer o caminho até a versão final:
#   ancora = 0/1           -> sem/com a restrição de soma zero nos professores atuais
#   prior_sigma_hist = 1/2 -> half-normal(0, 10) / Gama(2; 0,2) para sigma_hist
#   prof_nc = 0/1          -> delta dos professores atuais centrado / não centrado
#                             (tentativa 4: a half-Cauchy travava uma cadeia com
#                             o centrado; não centrado = delta = tau_prof * z)
# E mais uma chave, para a comparação de modelos (seção 9 da main.R):
#   com_efeitos = 0/1      -> sem/com os efeitos de escola e professor
# Ver a seção 4 da main.R (tentativas 1, 2 e 3).
# Parametrização não centrada (x = tau * x_raw) para evitar divergências.
# A truncagem em [0, 100] entra só na classificação (como no M1).

MODELO_STAN <- "
data {
  int<lower=1> N;                                  // alunos
  int<lower=1> T;                                  // avaliações
  matrix[N, T] y;                                  // notas do ano
  vector[N] nota_hist;                             // nota do histórico
  int<lower=1> n_esc;
  array[N] int<lower=1, upper=n_esc> escola;       // escola de origem
  int<lower=1> n_po;
  array[N] int<lower=0, upper=n_po> prof_orig;     // 0 = desconhecido
  int<lower=2> n_pa;
  array[N] int<lower=1, upper=n_pa> prof_atual;
  int<lower=1, upper=2> prior_tau;                 // 1 = half-normal; 2 = half-Cauchy
  int<lower=0, upper=1> ancora;                    // 1 = professores atuais somam zero
  int<lower=1, upper=2> prior_sigma_hist;          // 1 = half-normal; 2 = Gama(2; 0,2)
  int<lower=0, upper=1> prof_nc;                   // 1 = delta_a não centrado
  int<lower=0, upper=1> com_efeitos;               // 0 = sem escola e professor
  int<lower=0, upper=1> censura;                   // 1 = histórico 60 vale '60 ou menos'
}
parameters {
  vector[N] theta1;                  // proficiência na 1ª avaliação (chegada)
  matrix[N, T - 1] z;                // ruído padronizado da dinâmica
  vector[n_esc] gamma_raw;
  vector[n_po]  delta_o_raw;
  vector[n_pa - ancora] delta_a_livre; // com âncora, o último é definido pela soma zero
  real<lower=0> tau_esc;             // DP entre escolas
  real<lower=0> tau_prof;            // DP entre professores
  real<lower=0> omega;               // DP do ruído da dinâmica
  real<lower=0> sigma;               // DP do erro das avaliações
  real<lower=0> sigma_hist;          // DP do erro da nota do histórico
}
transformed parameters {
  // com_efeitos = 0 zera todos os efeitos (modelo de comparação)
  vector[n_esc] gamma   = com_efeitos * tau_esc  * gamma_raw;
  vector[n_po]  delta_o = com_efeitos * tau_prof * delta_o_raw;
  vector[n_pa]  delta_a_base;       // efeitos antes de multiplicar por tau
  vector[n_pa]  delta_a;
  if (ancora == 1) {
    // Âncora da escala: os professores atuais somam zero, ou seja,
    // theta é medido na régua do professor médio da escola atual
    delta_a_base = append_row(delta_a_livre, -sum(delta_a_livre));
  } else {
    delta_a_base = delta_a_livre;    // sem âncora: cada professor é livre
  }
  if (ancora == 1 && prof_nc == 0) {
    delta_a = delta_a_base;              // centrado (tentativas 2 e 3)
  } else {
    delta_a = tau_prof * delta_a_base;   // não centrado (tentativas 1 e 4)
  }
  delta_a = com_efeitos * delta_a;
  matrix[N, T] theta;                // theta_t = theta_{t-1} + omega * z_t
  theta[, 1] = theta1;
  for (t in 2:T) theta[, t] = theta[, t - 1] + omega * z[, t - 1];
}
model {
  // Priors
  theta1 ~ normal(50, 25);
  to_vector(z) ~ std_normal();
  gamma_raw   ~ std_normal();
  delta_o_raw ~ std_normal();
  if (ancora == 1 && prof_nc == 0) delta_a_base ~ normal(0, tau_prof);  // centrado
  else                             delta_a_base ~ std_normal();          // não centrado
  if (prior_tau == 1) {
    tau_esc  ~ normal(0, 10);
    tau_prof ~ normal(0, 10);
  } else {
    tau_esc  ~ cauchy(0, 5);
    tau_prof ~ cauchy(0, 5);
  }
  omega      ~ normal(0, 5);
  sigma      ~ normal(0, 10);
  if (prior_sigma_hist == 1) sigma_hist ~ normal(0, 10);
  else                       sigma_hist ~ gamma(2, 0.2);  // evita o zero: o histórico sempre tem erro

  // Nota do histórico: se o professor de origem é desconhecido,
  // a variância dele é somada (prior mais larga).
  // Censura: só é aprovado quem tem 60 ou mais; quem chega com 60 foi
  // 'aprovado no limite', então o histórico verdadeiro era 60 OU MENOS.
  for (i in 1:N) {
    real mu_h;
    real sd_h;
    if (prof_orig[i] > 0) {
      mu_h = theta1[i] + gamma[escola[i]] + delta_o[prof_orig[i]];
      sd_h = sigma_hist;
    } else {
      mu_h = theta1[i] + gamma[escola[i]];
      sd_h = sqrt(square(sigma_hist) + com_efeitos * square(tau_prof));
    }
    if (censura == 1 && nota_hist[i] <= 60)
      target += normal_lcdf(60 | mu_h, sd_h);   // P(histórico verdadeiro <= 60)
    else
      nota_hist[i] ~ normal(mu_h, sd_h);        // nota observada
  }

  // Avaliações do ano
  for (t in 1:T)
    y[, t] ~ normal(theta[, t] + delta_a[prof_atual], sigma);
}
"

# Organiza os dados no formato que o Stan pede (códigos viram números 1, 2, 3...)
dados_stan <- function(alunos, notas, prior_tau = 1, ancora = 1, prior_sigma_hist = 2,
                       prof_nc = 1, com_efeitos = 1, censura = 1) {
  esc <- factor(alunos$escola_origem)
  po  <- factor(alunos$prof_origem)                # NA continua NA
  pa  <- factor(alunos$prof_atual)
  y   <- matrix(notas$nota[order(notas$id, notas$aval)],
                nrow = nrow(alunos), byrow = TRUE)  # alunos em ordem de id
  list(
    N = nrow(alunos), T = ncol(y), y = y, nota_hist = alunos$nota_hist,
    n_esc = nlevels(esc), escola = as.integer(esc),
    n_po = nlevels(po), prof_orig = ifelse(is.na(po), 0L, as.integer(po)),
    n_pa = nlevels(pa), prof_atual = as.integer(pa),
    prior_tau = prior_tau, ancora = ancora, prior_sigma_hist = prior_sigma_hist,
    prof_nc = prof_nc, com_efeitos = com_efeitos, censura = censura,
    # guardados só para nomear os resultados depois
    nomes_esc = levels(esc), nomes_po = levels(po), nomes_pa = levels(pa)
  )
}

# Probabilidade de cada nível a partir das amostras do MCMC:
# é só a proporção de amostras de theta que caem em cada faixa
# (amostras abaixo de 0 contam como 0 e acima de 100 como 100 = truncagem)
m2_posterior <- function(fit, alunos) {
  th <- rstan::extract(fit, pars = "theta")$theta   # amostras x N x T
  th <- pmin(pmax(th, 0), 100)
  saida <- list()
  for (i in seq_len(nrow(alunos))) {
    for (t in seq_len(dim(th)[3])) {
      amostras <- th[, i, t]
      p <- table(nota_para_nivel(amostras)) / length(amostras)
      saida[[length(saida) + 1]] <- data.frame(
        id = alunos$id[i], aval = t, media = mean(amostras), s = sd(amostras),
        t(as.numeric(p)), check.names = FALSE)
    }
  }
  saida <- do.call(rbind, saida)
  names(saida)[5:8] <- NIVEIS
  saida
}

# Diagnóstico de um ajuste do Stan, em uma linha:
#   R-hat máximo  -> as cadeias concordam? (bom: < 1,01; aceitável: < 1,05)
#   ESS mínimo    -> quantas amostras "independentes" equivalentes (bom: > 400)
#   divergências  -> o amostrador "derrapou" em alguma região (bom: 0)
#   treedepth     -> quantas vezes o amostrador bateu no limite de passos (bom: 0)
# Também mostra o sigma_hist e os efeitos dos professores atuais, que é onde
# os problemas apareceram.
diagnostico_curto <- function(fit, nome, nomes_pa) {
  s <- summary(fit)$summary
  amostrador <- get_sampler_params(fit, inc_warmup = FALSE)
  treedepth  <- sum(sapply(amostrador, function(x) sum(x[, "treedepth__"] >= 10)))
  pa <- round(s[grep("^delta_a", rownames(s)), "mean"], 1)
  names(pa) <- paste0("delta_", nomes_pa)
  data.frame(
    versao        = nome,
    rhat_max      = round(max(s[, "Rhat"], na.rm = TRUE), 3),
    ess_min       = round(min(s[, "n_eff"], na.rm = TRUE)),
    divergencias  = get_num_divergent(fit),
    treedepth_max = treedepth,
    sigma_hist    = round(s["sigma_hist", "mean"], 2),
    ess_sigma_hist = round(s["sigma_hist", "n_eff"]),
    t(pa)
  )
}

# -----------------------------------------------------------------------------
# 7. Amostras do M2 em formato simples (uma linha por amostra, cadeias em ordem)
# -----------------------------------------------------------------------------
amostras_m2 <- function(fit, d) {
  pega <- function(p) as.matrix(fit, pars = p)
  sigma <- pega("sigma")[, 1]
  list(
    theta      = array(pega("theta"), dim = c(length(sigma), d$N, d$T)),
    gamma      = pega("gamma"),
    delta_o    = pega("delta_o"),
    delta_a    = pega("delta_a"),
    sigma      = sigma,
    omega      = pega("omega")[, 1],
    sigma_hist = pega("sigma_hist")[, 1],
    tau_prof   = pega("tau_prof")[, 1]
  )
}

# Média a posteriori de cada parâmetro, no mesmo formato (uma única "amostra")
media_amostras <- function(am) {
  list(
    theta      = array(apply(am$theta, c(2, 3), mean), dim = c(1, dim(am$theta)[2:3])),
    gamma      = matrix(colMeans(am$gamma), 1),
    delta_o    = matrix(colMeans(am$delta_o), 1),
    delta_a    = matrix(colMeans(am$delta_a), 1),
    sigma      = mean(am$sigma),
    omega      = mean(am$omega),
    sigma_hist = mean(am$sigma_hist),
    tau_prof   = mean(am$tau_prof)
  )
}

# -----------------------------------------------------------------------------
# 8. Verossimilhança do M2 (para DIC, WAIC e LOO)
# -----------------------------------------------------------------------------
# log p(dado | parâmetros) para cada amostra (linhas) e cada observação (colunas).
# São N x T notas do ano + N notas do histórico:
#   nota_it   ~ N(theta_it + delta_a, sigma^2)
#   nota_hist ~ N(theta_i1 + gamma + delta_o, sigma_hist^2)          (prof. conhecido)
#   nota_hist ~ N(theta_i1 + gamma, sigma_hist^2 + tau_prof^2)        (prof. desconhecido)
#   histórico = 60 (censurado): P(histórico verdadeiro <= 60) no lugar da densidade
# A verossimilhança completa é o produto de todas essas normais
# (no log, a soma de cada linha da matriz).
m2_loglik <- function(am, d) {
  S <- length(am$sigma)
  ll_ano  <- matrix(NA, S, d$N * d$T)
  ll_hist <- matrix(NA, S, d$N)
  k <- 0
  for (i in 1:d$N) {
    for (t in 1:d$T) {
      k <- k + 1
      media <- am$theta[, i, t] + am$delta_a[, d$prof_atual[i]]
      ll_ano[, k] <- dnorm(d$y[i, t], media, am$sigma, log = TRUE)
    }
    if (d$prof_orig[i] > 0) {
      media <- am$theta[, i, 1] + am$gamma[, d$escola[i]] + am$delta_o[, d$prof_orig[i]]
      dp    <- am$sigma_hist
    } else {
      media <- am$theta[, i, 1] + am$gamma[, d$escola[i]]
      dp    <- sqrt(am$sigma_hist^2 + d$com_efeitos * am$tau_prof^2)
    }
    ll_hist[, i] <- if (d$censura == 1 && d$nota_hist[i] <= 60)
      pnorm(60, media, dp, log.p = TRUE)                  # censurado
    else dnorm(d$nota_hist[i], media, dp, log = TRUE)     # observado
  }
  cbind(ll_ano, ll_hist)
}

# DIC como no Módulo VI (Spiegelhalter et al., 2002):
#   D(theta) = -2 log p(y | theta)
#   Dbar = média de D nas amostras;  Dhat = D na média a posteriori
#   pD = Dbar - Dhat (número efetivo de parâmetros);  DIC = Dbar + pD
dic <- function(ll, ll_media) {
  D_bar <- mean(-2 * rowSums(ll))
  D_hat <- -2 * sum(ll_media)
  c(Dbar = D_bar, Dhat = D_hat, pD = D_bar - D_hat, DIC = 2 * D_bar - D_hat)
}

# -----------------------------------------------------------------------------
# 8b. Verossimilhança POR ALUNO, com os theta integrados (para o LOO por aluno)
# -----------------------------------------------------------------------------
# Pergunta: "chegou um aluno novo; o modelo prevê bem as notas dele?"
# Para isso, os theta do aluno não podem ficar fixos: integramos sobre todas as
# trajetórias possíveis. Como tudo é normal e linear, a integral sai pelo
# mesmo filtro de Kalman do M1, observação por observação:
#   p(h, y1, ..., y9) = p(h) * p(y1 | h) * p(y2 | h, y1) * ...
# Cada fator é a "previsão" do Kalman antes de ver aquela nota.

# Um passo do Kalman (vale para vetores: uma posição por amostra do MCMC)
#   m, v   = média e variância atuais de theta
#   obs    = nota observada;  desloc = efeito que soma à nota (gamma, delta)
#   s2     = variância do erro da nota
passo_kalman <- function(m, v, obs, desloc, s2) {
  ll <- dnorm(obs, m + desloc, sqrt(v + s2), log = TRUE)  # previsão da nota
  k  <- v / (v + s2)                                       # ganho de Kalman
  list(ll = ll, m = m + k * (obs - desloc - m), v = (1 - k) * v)
}

# log p(histórico e 9 notas do aluno i | parâmetros): amostras x alunos
# m0, v0 = prior de theta_1 no Stan: N(50, 25^2)
#
# Histórico CENSURADO (= 60): ele não é uma normal, então a ordem muda:
#   p(y1..y9, h* <= 60) = p(y1..y9) * P(h* <= 60 | y1..y9)
# O 1º fator sai do Kalman só com as notas do ano. O 2º precisa de theta_1
# dado TODAS as notas: é o "suavizador" de Kalman (Rauch-Tung-Striebel), que
# volta da avaliação 9 até a 1 corrigindo a estimativa de theta_1.
# Conferido contra integração numérica (ver main.R, seção 9.4).
m2_loglik_aluno <- function(am, d, m0 = 50, v0 = 25^2) {
  S  <- length(am$sigma)
  ll <- matrix(0, S, d$N)
  for (i in 1:d$N) {
    # nota do histórico: h = theta_1 + gamma + delta_o + erro
    if (d$prof_orig[i] > 0) {
      desloc <- am$gamma[, d$escola[i]] + am$delta_o[, d$prof_orig[i]]
      s2     <- am$sigma_hist^2
    } else {
      desloc <- am$gamma[, d$escola[i]]
      s2     <- am$sigma_hist^2 + d$com_efeitos * am$tau_prof^2
    }
    censurado <- d$censura == 1 && d$nota_hist[i] <= 60
    m <- rep(m0, S); v <- rep(v0, S)
    if (!censurado) {                    # histórico observado: entra primeiro
      p <- passo_kalman(m, v, d$nota_hist[i], desloc, s2)
      ll[, i] <- p$ll; m <- p$m; v <- p$v
    }
    # notas do ano: y_t = theta_t + delta_a + erro;  theta_t = theta_{t-1} + ruído
    mf <- vf <- mp <- vp <- matrix(NA, S, d$T)   # filtrado e previsto, por avaliação
    for (t in 1:d$T) {
      if (t > 1) v <- v + am$omega^2              # o aluno pode ter mudado
      mp[, t] <- m; vp[, t] <- v
      p <- passo_kalman(m, v, d$y[i, t], am$delta_a[, d$prof_atual[i]], am$sigma^2)
      ll[, i] <- ll[, i] + p$ll; m <- p$m; v <- p$v
      mf[, t] <- m; vf[, t] <- v
    }
    if (censurado) {                     # P(h* <= 60 | notas do ano)
      ms <- mf[, d$T]; vs <- vf[, d$T]
      for (t in (d$T - 1):1) {           # suavizador: volta até a avaliação 1
        J  <- vf[, t] / vp[, t + 1]
        ms <- mf[, t] + J * (ms - mp[, t + 1])
        vs <- vf[, t] + J^2 * (vs - vp[, t + 1])
      }
      ll[, i] <- ll[, i] + pnorm(60, ms + desloc, sqrt(vs + s2), log.p = TRUE)
    }
  }
  ll
}

# Verificação preditiva a posteriori: simula dados novos a partir de n amostras
# da posterior. Devolve uma lista com:
#   y: notas do ano simuladas (n x N*T), na ordem de as.vector(t(d$y))
#      (aluno 1 aval 1..9, aluno 2 aval 1..9, ...)
#   h: notas do histórico simuladas (n x N)
simular_rep <- function(am, d, n = 200) {
  qual <- round(seq(1, length(am$sigma), length.out = n))
  y_rep <- matrix(NA, n, d$N * d$T)
  h_rep <- matrix(NA, n, d$N)
  conhecido <- d$prof_orig > 0
  for (r in seq_along(qual)) {
    s <- qual[r]
    media <- am$theta[s, , ] + am$delta_a[s, d$prof_atual]   # N x T
    y_rep[r, ] <- as.vector(t(media + matrix(rnorm(d$N * d$T, 0, am$sigma[s]), d$N)))
    # histórico: professor conhecido usa o delta dele; desconhecido soma a variância
    media_h <- am$theta[s, , 1] + am$gamma[s, d$escola]
    media_h[conhecido] <- media_h[conhecido] + am$delta_o[s, d$prof_orig[conhecido]]
    dp_h <- ifelse(conhecido, am$sigma_hist[s],
                   sqrt(am$sigma_hist[s]^2 + d$com_efeitos * am$tau_prof[s]^2))
    h_rep[r, ] <- rnorm(d$N, media_h, dp_h)
    if (d$censura == 1) h_rep[r, ] <- pmax(h_rep[r, ], 60)   # regra de aprovação
  }
  list(y = y_rep, h = h_rep)
}

# -----------------------------------------------------------------------------
# 8c. Checagem preditiva A PRIORI (Tópico 9)
# -----------------------------------------------------------------------------
# Antes de olhar os dados: sorteia os parâmetros SÓ das priors do M2 final e
# simula notas com a mesma estrutura do estudo (mesmos alunos, escolas e
# professores). Pergunta: as notas geradas são plausíveis (entre 0 e 100)?
prior_preditiva <- function(d, n = 1000) {
  y_prior <- matrix(NA, n, d$N * d$T)
  conhecido <- d$prof_orig > 0
  for (r in 1:n) {
    tau_esc  <- abs(rnorm(1, 0, 10));  tau_prof <- abs(rnorm(1, 0, 10))
    omega    <- abs(rnorm(1, 0, 5));   sigma    <- abs(rnorm(1, 0, 10))
    gamma    <- rnorm(d$n_esc, 0, tau_esc)
    z        <- rnorm(d$n_pa);  delta_a <- tau_prof * (z - mean(z))  # soma zero
    theta1   <- rnorm(d$N, 50, 25)
    theta    <- theta1 + t(apply(matrix(rnorm(d$N * d$T, 0, omega), d$N), 1,
                                 function(x) cumsum(c(0, x[-1]))))     # passeio
    y <- theta + delta_a[d$prof_atual] + matrix(rnorm(d$N * d$T, 0, sigma), d$N)
    y_prior[r, ] <- as.vector(t(y))
  }
  y_prior
}

# -----------------------------------------------------------------------------
# 9. Diagnósticos de convergência do pacote coda (os mesmos do Módulo V/VI)
# -----------------------------------------------------------------------------
diagnosticos_coda <- function(fit, pars) {
  mc <- As.mcmc.list(fit, pars = pars)             # Stan -> formato do coda
  gel <- gelman.diag(mc, multivariate = FALSE)$psrf
  gew <- sapply(geweke.diag(mc), function(g) g$z)  # um z por cadeia
  hei <- sapply(heidel.diag(mc), function(h) h[, "stest"])  # 1 = passou
  hpd <- HPDinterval(as.mcmc(as.matrix(mc)), prob = 0.95)
  # MCSE: erro do próprio MCMC na média (= DP da posterior / raiz do ESS)
  mcse <- summary(mc)$statistics[, "Time-series SE"]
  tab <- data.frame(
    parametro    = colnames(as.matrix(mc)),
    media        = colMeans(as.matrix(mc)),
    HPD_inf      = hpd[, "lower"],
    HPD_sup      = hpd[, "upper"],
    ESS          = effectiveSize(mc),
    MCSE         = mcse,
    gelman_rubin = gel[, "Point est."],
    geweke_max_z = apply(abs(gew), 1, max),         # |z| > 1,96 indicaria problema
    heidel_passou = apply(hei, 1, function(x) paste0(sum(x == 1), "/", length(x)))
  )
  rownames(tab) <- NULL
  # Raftery-Lewis: exige pelo menos 3746 amostras por cadeia (q = 0,025; r = 0,005)
  raftery <- tryCatch(raftery.diag(mc[[1]]), error = function(e) conditionMessage(e))
  list(tabela = tab, raftery = raftery, mc = mc)
}

# -----------------------------------------------------------------------------
# 10. Teste de hipóteses bayesiano (Tópico 8)
# -----------------------------------------------------------------------------
# Para um efeito x (amostras da posterior):
#   P(x > 0 | dados), P(x < 0 | dados) e
#   P(|x| < rope | dados): região de equivalência prática (ROPE), "efeito desprezível"
teste_hipotese <- function(x, nome, rope = 2) {
  data.frame(parametro = nome, media = round(mean(x), 2),
             P_maior_0 = round(mean(x > 0), 3), P_menor_0 = round(mean(x < 0), 3),
             P_ROPE = round(mean(abs(x) < rope), 3))
}
