# =============================================================================
# simulacao.R — Gera os dados simulados do estudo
# Classificação bayesiana de proficiência em Biologia (níveis tipo PAEBES)
# JPGTEIXEIRA — EST613 Estatística Bayesiana / UFV
#
# Saídas (pasta data/):
#   alunos.csv   -> o que a escola sabe na chegada do aluno
#   notas.csv    -> as 9 avaliações do ano (formato longo)
#   verdade.csv  -> valores verdadeiros (só existe na simulação)
# =============================================================================

set.seed(2026)
dir.create("data", showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 1. Parâmetros da simulação (tudo que é "verdade" fica declarado aqui)
# -----------------------------------------------------------------------------

n_por_perfil <- 10            # 12 perfis x 10 = 120 alunos
n_aval       <- 9             # 3 trimestres x (2 provas + 1 trabalho)
sd_entrada   <- 5             # variação individual em torno do nível de entrada
sd_passeio   <- 2             # pequeno ruído na trajetória (aprendizagem não é linear)
sd_hist      <- 5             # erro da nota do histórico (série anterior)
sd_prova     <- 6             # erro de medida de cada avaliação do ano
p_prof_conhecido <- 0.4       # chance de saber quem foi o professor de origem

# Nível de entrada (proficiência real no início do ano, escala 0-100)
entrada <- c(fraco = 35, medio = 55, forte = 85)

# Variação total ao longo do ano por trajetória (da aval. 1 à 9)
# (o forte evolui menos e cai mais, para não estourar o teto de 100)
variacao <- list(
  evolui   = c(fraco = 20,  medio = 20,  forte = 10),
  involui  = c(fraco = -15, medio = -15, forte = -20),
  estagna  = c(fraco = 0,   medio = 0,   forte = 0)
)
amp_gangorra <- 10            # amplitude da oscilação (± pontos)

# Efeito das escolas de origem (gamma): E5 infla nota
gamma_escola <- c(E1 = 0, E2 = -4, E3 = 3, E4 = 0, E5 = 12)

# Efeito dos professores de origem (delta_origem): 3 por escola
prof_origem_efeito <- c(
  E1a =  0, E1b = -8, E1c =  4,
  E2a =  2, E2b =  0, E2c = -5,
  E3a =  6, E3b = -3, E3c =  0,
  E4a = -10, E4b = 3, E4c =  5,
  E5a =  4, E5b =  0, E5c =  8
)

# Professores da escola de destino (delta_atual), nomeados pelo rigor:
#   PR = rigoroso (tira nota), PM = mediano (neutro), PF = facilita (dá nota)
prof_atual_efeito <- c(PR = -10, PM = 0, PF = 8)

# Calendário das avaliações
calendario <- data.frame(
  aval      = 1:n_aval,
  trimestre = rep(1:3, each = 3),
  tipo      = rep(c("prova", "prova", "trabalho"), times = 3)
)

# -----------------------------------------------------------------------------
# 2. Alunos: perfil, origem e professores
# -----------------------------------------------------------------------------

perfis <- expand.grid(
  trajetoria = c("evolui", "involui", "estagna", "gangorra"),
  entrada    = names(entrada),
  stringsAsFactors = FALSE
)
perfis <- perfis[rep(seq_len(nrow(perfis)), each = n_por_perfil), c("entrada", "trajetoria")]

n <- nrow(perfis)
alunos <- data.frame(
  id         = sprintf("A%03d", 1:n),
  entrada    = perfis$entrada,
  trajetoria = perfis$trajetoria,
  perfil     = paste(perfis$entrada, perfis$trajetoria, sep = "_"),
  stringsAsFactors = FALSE
)

# Origem sorteada ao acaso (independente do perfil)
alunos$escola_origem <- sample(names(gamma_escola), n, replace = TRUE)
alunos$prof_origem_real <- paste0(alunos$escola_origem,
                                  sample(c("a", "b", "c"), n, replace = TRUE))
alunos$prof_atual <- sample(names(prof_atual_efeito), n, replace = TRUE)

# Professor de origem só é conhecido para parte dos alunos
conhecido <- runif(n) < p_prof_conhecido
alunos$prof_origem <- ifelse(conhecido, alunos$prof_origem_real, NA)

# -----------------------------------------------------------------------------
# 3. Trajetória verdadeira da proficiência (theta) ao longo das 9 avaliações
# -----------------------------------------------------------------------------

limitar <- function(x) pmin(pmax(x, 0), 100)

theta_inicial <- entrada[alunos$entrada] + rnorm(n, 0, sd_entrada)
t_rel <- (0:(n_aval - 1)) / (n_aval - 1)        # 0 na aval. 1, 1 na aval. 9

theta <- matrix(NA, n, n_aval)
for (i in 1:n) {
  if (alunos$trajetoria[i] == "gangorra") {
    # sobe e desce duas vezes no ano
    tendencia <- amp_gangorra * sin(2 * pi * 2 * t_rel)
  } else {
    tendencia <- variacao[[alunos$trajetoria[i]]][alunos$entrada[i]] * t_rel
  }
  passeio <- cumsum(c(0, rnorm(n_aval - 1, 0, sd_passeio)))
  theta[i, ] <- limitar(theta_inicial[i] + tendencia + passeio)
}

# -----------------------------------------------------------------------------
# 4. Nota do histórico (série anterior): theta inicial + viés da origem
# -----------------------------------------------------------------------------

alunos$gamma_origem <- gamma_escola[alunos$escola_origem]
alunos$delta_origem <- prof_origem_efeito[alunos$prof_origem_real]
alunos$nota_hist <- round(limitar(
  theta_inicial + alunos$gamma_origem + alunos$delta_origem + rnorm(n, 0, sd_hist)
), 1)
# Regra de aprovação: só passa de série quem tem 60 ou mais. Quem ficaria abaixo
# chega "aprovado no limite", com 60 (recuperação / conselho de classe).
alunos$nota_hist <- pmax(alunos$nota_hist, 60)

# -----------------------------------------------------------------------------
# 5. Notas do ano: theta do momento + rigor do professor atual + erro
#    (prova e trabalho simulados da mesma forma; a coluna tipo fica só como info)
# -----------------------------------------------------------------------------

alunos$delta_atual <- prof_atual_efeito[alunos$prof_atual]

notas <- merge(data.frame(id = rep(alunos$id, each = n_aval),
                          aval = rep(1:n_aval, times = n)),
               calendario, by = "aval")
notas <- notas[order(notas$id, notas$aval), ]
idx <- match(notas$id, alunos$id)
theta_long <- theta[cbind(idx, notas$aval)]
notas$nota <- round(limitar(
  theta_long + alunos$delta_atual[idx] + rnorm(nrow(notas), 0, sd_prova)
), 1)

# -----------------------------------------------------------------------------
# 6. Gravação
# -----------------------------------------------------------------------------

# O que a escola sabe (sem nenhum valor verdadeiro)
write.csv(alunos[, c("id", "escola_origem", "prof_origem", "nota_hist", "prof_atual")],
          "data/alunos.csv", row.names = FALSE, na = "NA")

write.csv(notas[, c("id", "aval", "trimestre", "tipo", "nota")],
          "data/notas.csv", row.names = FALSE)

# Verdade: só para avaliar o modelo no final
verdade <- data.frame(
  id           = notas$id,
  perfil       = alunos$perfil[idx],
  entrada      = alunos$entrada[idx],
  trajetoria   = alunos$trajetoria[idx],
  prof_origem_real = alunos$prof_origem_real[idx],
  gamma_origem = alunos$gamma_origem[idx],
  delta_origem = alunos$delta_origem[idx],
  delta_atual  = alunos$delta_atual[idx],
  aval         = notas$aval,
  theta_real   = round(theta_long, 2)
)
write.csv(verdade, "data/verdade.csv", row.names = FALSE)

cat("Simulação concluída:", n, "alunos x", n_aval, "avaliações\n")
