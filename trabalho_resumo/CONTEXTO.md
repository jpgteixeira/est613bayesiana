# Contexto do projeto — Classificação bayesiana de proficiência (níveis PAEBES)

**Autor:** João Paulo Gusmão Teixeira (JPGTEIXEIRA)
**Disciplina:** EST613 – Estatística Bayesiana (UFV), Profa. Camila Ferreira Azevedo
**Atualizado em:** 21/11/2026

> Este arquivo resume tudo o que já foi decidido. Use-o para retomar o trabalho numa conversa nova.

---

## 1. Prazos da disciplina

| Entrega | Data |
|---|---|
| Resumo em tópicos com o passo a passo de uma análise bayesiana | **26/11/2026 (qui)** |
| Resumo expandido (máx. 5 páginas) | **30/11/2026 (seg)** |
| Apresentação oral (10–15 min) — João Paulo é o **1º** | **30/11/2026 (seg)** |

Seções do resumo expandido: Título e identificação · Introdução · Objetivo · Metodologia · Resultados e discussão · Conclusão · Referências.

Na apresentação: contexto breve, objetivo, escolhas metodológicas; **a maior parte do tempo vai para os resultados**; no fim, retomar a pergunta e concluir.

---

## 2. Ideia central

Um modelo bayesiano que classifica os estudantes de Biologia nos níveis do PAEBES (**Abaixo do Básico, Básico, Proficiente, Avançado**) e **se atualiza a cada avaliação do ano**, permitindo que o aluno mude de nível.

**Pergunta do trabalho:** o modelo bayesiano dinâmico, que considera a escola de origem e o rigor do professor, classifica os alunos melhor que a média simples e detecta mais cedo quando o aluno muda de nível?

Princípio prático: o modelo se adapta à rotina da escola (notas de 0–100 lançadas pelos professores), e não o contrário.

---

## 3. Decisões tomadas

### Escala
- Fica em **0–100** neste trabalho. A ligação com a escala do PAEBES (média 250, DP 50, TRI) fica como **trabalho futuro / publicação**, com dados reais.
- Motivos: com dados simulados, converter para a escala do PAEBES é só uma transformação linear e não muda nenhuma conclusão; o PAEBES de Ciências da Natureza é pontual (ex.: 2019, 3ª série EM), então não há uma prior de Biologia para a 1ª série; os pontos de corte que o Gemini sugeriu cobrem só Língua Portuguesa e Matemática e não foram conferidos.
- Notas com ponto extra: usar a proporção da nota máxima, limitada a 100.

### Modelo
- **Proficiência do aluno (θ):** normal **truncada em [0, 100]**. Com prior normal truncada e verossimilhança normal, a posterior é a mesma normal de antes, só que truncada no mesmo intervalo, então a fórmula fechada continua valendo.
- **Prior:** nota do histórico (série anterior). É legítimo usá-la como prior, porque é informação que existe antes dos dados novos (prior informativa empírica).
- **Na chegada (prior):**
  `nota_hist = θ + γ_escola_origem + δ_prof_origem + erro`
  - A escola de origem é sempre conhecida e corrige a escola que infla nota.
  - O professor de origem só às vezes é conhecido (em cidade pequena, sim; em cidade grande, raramente). Quando não é conhecido, **a variância do efeito do professor é somada à variância da prior** (a prior fica mais larga). Não precisa de nenhuma regra especial para o dado faltante.
- **Durante o ano (verossimilhança):**
  `nota_prova = θ_t + δ_prof_atual + erro`
  - O professor atual é sempre conhecido. A escola atual não entra no modelo (é uma escola só).
- **Dinâmica:** antes de cada avaliação, soma-se uma pequena variância à proficiência (θ_t = θ_{t−1} + ruído; é um filtro de Kalman simples). Isso faz as provas recentes pesarem mais e permite que o aluno mude de nível.
- **Rigor do professor:** é modelado como efeito δ. Isso é melhor que padronizar por turma (escore Z), porque o Z confunde professor rigoroso com turma fraca. Para conseguir separar os dois, é preciso ligação entre os dados: alunos que trocam de professor, professores com várias turmas, ou uma régua externa comum (o PAEBES, no futuro).
- **Classificação:** probabilidade posterior de o aluno estar em cada faixa.

### Três modelos comparados
| Modelo | Descrição |
|---|---|
| **M0** | Média simples das notas (o que a escola faz hoje) |
| **M1** | Fórmula fechada Normal–Normal dinâmica, **sem** efeitos de escola e professor |
| **M2** | Stan, com γ_escola e δ_professor **estimados** (hierárquico, professor aninhado na escola) |

- Priors das variâncias no M2: **half-normal(0, 10)** em τ_escola e τ_prof como modelo principal; **half-Cauchy** como análise de sensibilidade (Gelman 2006; Dias et al. 2022, TAG 135:1385–1399). Com 2 ou 3 professores, a cauda da Cauchy pode deixar a posterior instável.
- O M1 é rápido, então dá para repetir a simulação umas 100 vezes. O M2 roda em uma réplica representativa e mostra os diagnósticos (R-hat, ESS, MCSE, divergências), que é o conteúdo da EST613 (Módulo V usa rstan, rjags, coda, bayesplot, truncnorm).
- **M2 como implementado (19/11):** código Stan dentro da `funcoes.R` (objeto `MODELO_STAN`). Um único τ_prof para professores de origem e atuais. **Âncora:** δ dos professores atuais somam zero (sem isso, uma constante passa de θ para os efeitos e os δ saem ~+5 deslocados). θ_1 ~ N(50, 25²) ligado direto à nota do histórico; dinâmica θ_t = θ_{t−1} + ω·z_t (não centrada). Priors: τ ~ half-normal(0,10) [sensib.: half-Cauchy(0,5)]; σ ~ half-normal(0,10); ω ~ half-normal(0,5); **σ_hist ~ Gama(2; 0,2)** (evita o zero; com half-normal ficava preso perto de 0). 4 cadeias × 2000 (1000 warmup), adapt_delta 0,99, seed 2026. ~4 min para os dois ajustes; salvos em `results/fit_M2_*.rds` (90 MB cada, não versionar). Probabilidade dos níveis no M2 = proporção de amostras de θ em cada faixa.
- **Caminho do M2 fica no código** (main.R, seção 4; chaves no Stan: `ancora`, `prior_sigma_hist`, `prof_nc`, `com_efeitos`), uma mudança por vez, adapt_delta 0,99: T1 ingênuo (3 div., R-hat 1,04, δ PR/PM/PF = −4,9/4,7/11,9) → T2 + âncora (6 div., δ = −8,8/0,8/7,9) → T3 + Gama em σ_hist (0 div.; mas com half-Cauchy uma cadeia travou com τ_prof≈0, R-hat 5,0) → **T4 final: δ_atual não centrado** (0 div., R-hat ≤ 1,01 nas duas priors). Tabela em `results/tentativas_M2.csv`. João Paulo quer ver o processo, não só a solução final.
- **Avaliação conforme a disciplina (19/11, conferido no Drive da EST613):** verossimilhança explícita (função `m2_loglik`); diagnósticos do coda (Gelman-Rubin, Geweke, Heidelberger-Welch, ESS, autocorr, HPD; Raftery-Lewis não se aplica: exige ≥3746 amostras/cadeia); DIC (como Módulo VI), WAIC e LOO comparando M2 × M2 sem efeitos (`com_efeitos = 0`); verificação preditiva (y_rep, p-valores bayesianos); teste de hipóteses (P(δ<0|y), ROPE ±2). A professora usa JAGS com 100.000 iterações; com Stan/NUTS 2.000 bastam (critério = ESS/MCSE).
- O tipo de avaliação (prova × trabalho) **não entra no modelo** por enquanto; fica só como comentário na apresentação (trabalho costuma inflar a nota e diferenciar menos os alunos).

### Estrutura do código (simples, poucos arquivos)
```
simulacao.R   -> gera data/ (roda uma vez, set.seed(2026))
funcoes.R     -> funções (níveis, M0, M1/Kalman, regras de classificação, métricas)
main.R        -> source("funcoes.R"), lê data/, salva em results/
data/         -> alunos.csv, notas.csv, verdade.csv
results/      -> figuras e tabelas
docs/         -> PDFs de referência (PAEBES)
presentation/ -> apresentacao.qmd (Quarto revealjs; figuras lidas de ../results/)
resumo/       -> resumo_expandido.tex (LaTeX, pdflatex; figuras de ../results/)
```
A `main.R` **lê os CSV**, não chama a simulação. Com dados reais da Sedu, basta trocar `alunos.csv` e `notas.csv` (mesmo formato); o `verdade.csv` não existe com dados reais.

---

## 4. Dados simulados (já gerados — `simulacao.R`)

- **120 alunos** = 12 perfis × 10 alunos; **9 avaliações** (3 trimestres × 2 provas + 1 trabalho) = 1.080 notas.
- **Entrada:** fraco ~35, médio ~55, forte ~85 (DP individual 5).
- **Trajetórias (variação da aval. 1 à 9):** evolui (+20; forte +10), involui (−15; forte −20), estagna (0), gangorra (±10, sobe e desce duas vezes no ano). Há um pequeno passeio aleatório (DP 2) em todas.
- **Escolas de origem (γ):** E1 = 0, E2 = −4, E3 = +3, E4 = 0, **E5 = +12 (infla)**.
- **Professores de origem (δ):** 3 por escola (E1a…E5c), entre −10 e +8. O professor é conhecido para **~40%** dos alunos; nos demais aparece `NA`.
- **Professores atuais (nomes mnemônicos, 19/11):** **PR = −10 (rigoroso)**, **PM = 0 (mediano)**, **PF = +8 (facilita)**. (Antes PA/PB/PC; dados idênticos, só o rótulo mudou.)
- **Erros:** histórico DP 5; avaliações DP 6. Notas limitadas a [0, 100].
- A origem (escola e professor) é sorteada de forma independente do perfil.

### Formato dos arquivos
- `alunos.csv`: `id, escola_origem, prof_origem, nota_hist, prof_atual`
- `notas.csv`: `id, aval, trimestre, tipo, nota` (formato longo)
- `verdade.csv`: `id, perfil, entrada, trajetoria, prof_origem_real, gamma_origem, delta_origem, delta_atual, aval, theta_real`

---

## 5. Pontos ainda em aberto

1. ~~Pontos de corte~~ → **decidido (19/11): 40 / 60 / 80.**
2. ~~Critério de reclassificação~~ → **decidido (19/11): mostrar os dois.** Principal = troca só quando P > 0,7; comparação = maior probabilidade.
3. ~~Ruído da dinâmica~~ → fixo no M1 (sd_ruido = 3; sd_prior = 10; sd_prova = 6); **estimado no M2** (ω ≈ 3,9).
4. ~~Escopo~~ → o M2 entra inteiro (é o melhor modelo e traz os diagnósticos).
5. Métricas usadas: acerto, trocas de nível por aluno, % de mudanças reais não detectadas, recuperação de γ e δ. O "atraso médio" ficou ambíguo com uma réplica só (fica de fora do texto).
6. **Em aberto:** repetir a simulação ~100× (só M0/M1, ou M2 em poucas réplicas)? Vale a pena até 30/11?

## 5b. Resultados (19/11, uma réplica, modelo final T4)

| Modelo | Acerto (%) | Trocas/aluno | % mudanças não detectadas |
|---|---|---|---|
| M0 | 66,4 | 0,58 | 48,1 |
| M1 maior prob. | 69,6 | 0,98 | 32,1 |
| M1 limiar 0,7 | 70,3 | 0,55 | 39,5 |
| M2 maior prob. | 85,6 | 0,68 | 29,6 |
| M2 limiar 0,7 | 84,8 | 0,38 | 37,0 |
| M2 limiar (half-Cauchy) | 84,9 | 0,39 | 34,6 |

Trocas reais: 1,55/aluno; 81 de 120 mudam de nível. Acerto por professor (M0 / M1 lim / M2 lim): PR 55,3 / 59,6 / 82,2 · PM 79,2 / 80,9 / 86,0 · PF 61,7 / 67,9 / 86,1.
Recuperação (HPD 95%): γ_E5 = 13,8 (10,4–17,1; real 12); δ PR/PM/PF = −8,8 / 0,8 / 7,9 (real −10 / 0 / 8).
Convergência (coda): Gelman-Rubin ≤ 1,01; ESS mín. > 350; Heidelberger aprovado; Geweke 5/52 com |z|>1,96.
Comparação (menor = melhor): DIC 8167 × 8383. **LOO por aluno (opção B, 19/11):** a verossimilhança de cada aluno integra os θ pelo filtro de Kalman do M1 (função `m2_loglik_aluno`, conferida contra a normal multivariada: igual até 1e-14); ΔELPD = 82,3 (EP 9,1); completo prevê melhor 104/120 alunos; Pareto-k > 0,7 só em 3 alunos (A070, A092, A094: professor de origem que só aparece neles); sem eles ΔELPD = 85,3. O LOO/WAIC por nota (θ fixo) ficou no código como 1ª tentativa descartada: 68/1200 com k > 0,7.
**Tópico 9 (19/11, Drive):** checklist feito.
- Checagem preditiva **a priori** (`prior_preditiva`): mediana 50, mas 9,1% das notas fora de 0–100 (p1 = −20, p99 = 120). **Priors NÃO alteradas** (João Paulo decide); mesma questão do teto.
- Verificação preditiva **a posteriori em 3 níveis** (`simular_rep`, agora com h_rep): global média 0,40, DP 0,53, mín. 0,20, **máx. 1,00** (sem teto de 100); grupo professor 0,46, escola (histórico) 0,32; individual (DP das 9 notas do aluno) 0,74. Figura `fig7_ppc_niveis.png`.
- MCSE na tabela do coda (≤ 0,06).
- **Partial pooling** dos professores de origem: estimativa M2 = fator × "sem pooling"; fator 0,70 (1 aluno) → 0,94 (7 alunos), igual à teoria τ²/(τ²+σ_h²/n). Comparar com o valor *verdadeiro* NÃO mostra o shrinkage (mistura ruído do histórico e confusão escola×professor da E5). Figura `fig10_shrinkage.png`.
- Mensagem "convergência ≠ adequação": M2 converge e é o melhor no LOO, mas falha no máximo.
- Identificabilidade (aula 19/11): tentativa 1 = modelo não identificável, como y = μ + τ_i + e → restrição Σδ = 0; citado no slide e no resumo.
- Resumo agora com 5 páginas (só referências na 5ª).
- **Priors: decisão final (19/11) = manter como estão e mencionar o teto de 100 como limitação.**
- **Apresentação revisada (19/11):** sem o checklist; slides que estouravam corrigidos (`{.smaller}`, legendas menores, título menor via CSS no próprio .qmd), conferido com render + capturas. Figuras novas: `fig0_desenho.png` (escolas pretas e E5 vermelha, 12 trajetórias, professores como carinhas feliz/neutra/triste = PF/PM/PR) no slide de dados simulados; `fig11_prob_perfis.png` (12 painéis entrada × trajetória, probabilidade média de cada nível, "tipo STRUCTURE") no lugar do slide do aluno gangorra. João Paulo vai passar rápido por vários slides (convergência, comparação etc.).
- **Dúvidas conceituais esclarecidas (19/11):** (1) prior θ_i1 ~ N(50, 25²) pesa ~2% frente ao histórico (6% com professor de origem desconhecido); quem informa o aluno é o histórico corrigido (h − γ − δ), que faz o papel de prior; ele entra "como dado" porque γ e δ também são estimados. No M1, o histórico é a média da prior sem correção. (2) As notas são **censuradas** em 100 (28/1080 notas = 100; 2 = 0; 5 históricos = 100), não truncadas; θ é truncado só na classificação. Censura no Stan tiraria a fórmula fechada do LOO por aluno (Kalman) → mantida como limitação. (3) Hiperparâmetros = τ_esc, τ_prof (hiperpriori half-N(0,10)). Extensão possível: escola também na prior do aluno, θ_i1 ~ N(μ + u_escola, σ²) (escola que forma alunos melhores, ≠ γ, que é viés da nota; identificável); não incluída porque na simulação a escola é sorteada sem relação com o perfil.
- **Apresentação:** 23 slides, todos com notas de fala detalhadas (`::: {.notes}`, tecla S no revealjs); slide 6 virou tabela "parâmetro · o que é · prior", com τ marcado como hiperparâmetro.
Hipóteses: P(δ_PR<0|y)=1,00; P(δ_PF>0|y)=1,00; P(|δ_PM|<2|y)=0,92; P(γ_E5>0|y)=1,00.

## 6. Próximos passos

1. ~~Fechar os pontos 1 e 2~~ e ~~escrever `funcoes.R`/`main.R` (M0 e M1)~~ — feito em 19/11.
2. ~~M2 em Stan + diagnósticos~~, ~~figuras/tabelas~~, ~~1ª versão da apresentação e do resumo expandido (4 págs.)~~ — feito em 19/11.
3. ~~Verossimilhança explícita, diagnósticos coda + HPD, DIC/WAIC/LOO, verificação preditiva, teste de hipóteses, nomes PR/PM/PF~~ — feito em 19/11; resumo (4 págs.) e apresentação atualizados.
4. João Paulo instalando Stan (rstan) no VS Code para rodar a `main.R` localmente (~15 min na 1ª vez: 7 ajustes do Stan).
5. Revisar texto do resumo e da apresentação; decidir sobre as réplicas (ponto 6); possível melhoria: teto de 100 nas notas (censura).
6. Resumo das etapas bayesianas (26/11): usar o próprio trabalho como roteiro (prior → verossimilhança → posterior → MCMC → diagnóstico → inferência).

## 7. Trabalho futuro (publicação)
- Dados reais anonimizados da Sedu (notas do 9º ano e bimestrais da 1ª série), com termo de anuência e aprovação de comitê de ética.
- Ligação com a escala do PAEBES / avaliação diagnóstica como régua externa comum, para separar o rigor do professor do nível da turma.

## 21/11 — apresentação mais "entendível"
- Slide 6 virou "Como o modelo funciona": figura do ciclo (results/fig13_ciclo.png, main.R 11.13): HISTÓRICO → PRIOR × PROVA t = POSTERIORI → (+ incerteza ω) → prior da próxima prova; em vermelho, onde o M2 desconta escola (no histórico) e professor (em cada prova). Tabela de parâmetros foi para as notas do slide.
- Slide de perfis virou "Um aluno de cada perfil" (results/fig12_prob_alunos.png, main.R 11.12): 12 alunos individuais (A001, A011, ..., A111), faixa de cima = nível real, área = probabilidade do M2. fig11 (médias por perfil) continua sendo gerada, mas não está na apresentação.
- Escolas: um γ por escola, todas vindas da mesma N(0, τ_esc²); não precisa classificar escolas antes. Extensão com dados reais: covariável de escola na média (γ_e ~ N(β·tipo_e, τ_esc²)).

## 21/11 — revisão final (código, apresentação, resumo)
- Slides novos: "O modelo em 5 passos" (M1 × M2: parâmetro, prior, verossimilhança, dinâmica, posteriori, classificação) e "Distribuições e valores (M2)" (prior assumida × média a posteriori × valor real: σ 6,3; σ_h 3,5; ω 3,9; τ_esc 7,9; τ_prof 5,4). Saiu o slide "Modelos comparados e priors". 26 slides.
- Considerações finais (resumo + slide de conclusão): histórico 6º–9º ano (prior do 1º ano = posteriori do fim do EF) e PAEBES como régua externa (PAEBES_i ~ N(a + b·θ_i, σ_P²); calibrar efeitos e prever o padrão do aluno).
- Revisão: termo "censura em 100" adotado na limitação; ESS mín. corrigido para > 480; Heidelberger aprovado exceto ω em 1 de 4 cadeias (half-Normal); fronteira dos níveis (40, 60] consistente com o código; NUTS justificado (MH com gradiente; >1.000 parâmetros correlacionados, priors não conjugadas) no resumo e nas notas; comentário PA/PB/PC → PR/PM/PF na main.R.
- Resumo: 5 páginas, conclusão termina na p. 4, só referências na p. 5.
- Ideia conectada (fora do trabalho): Calibra (CAT/TRI 3PL) — este modelo daria a prior da diagnóstica e a dinâmica entre sessões; o Calibra daria a régua comum; δ do professor ≈ dificuldade b do item.
- Pendente: resumo em tópicos (26/11); réplicas da simulação descartadas (limitação "uma réplica" já citada).

## 23/11 — histórico censurado em 60 (aprovado pela professora)
- **Simulação:** só é aprovado quem tem histórico >= 60; quem ficaria abaixo chega com 60 ("aprovado no limite"), mantendo a proficiência real. 60 de 120 alunos. Notas do ano e proficiências idênticas à versão anterior.
- **M2:** histórico = 60 entra como censurado: `target += normal_lcdf(60 | mu_h, s_h)` (flag `censura` no Stan e em `dados_stan`). M1 usa o histórico cru. DIC, PPC (h_rep = max(h, 60)) e LOO por aluno atualizados; LOO por aluno censurado = p(y)·P(h* <= 60 | y), com suavizador de Kalman (RTS), conferido contra integração numérica. Partial pooling só com histórico observado.
- **M2 ingênuo** (60 como nota exata, `censura = 0`): classifica igual (82,4% × 82,8%), mas conclui que todas as escolas inflam (E1 +6,7; E2 +12,6; E3 +8,9; E4 +12,5; E5 +17,1) e σ_h = 9,5. Com censura: −2,0 / −0,9 / +1,8 / +2,2 / +11,5 (real 0/−4/3/0/12), σ_h = 4,0. Figura `fig14_censura.png`, tabela `censura_vs_ingenuo.csv`.
- **Resultados novos:** acerto M0 66,4; M1 lim 67,9; M2 lim 82,8; M2 max 84,7; hc 82,7. Por professor (M0/M1/M2): PR 55,3/59,1/79,2 · PM 79,2/79,0/85,3 · PF 61,7/63,0/83,3. Efeitos: E5 +11,5 (7,0–15,6); δ PR/PM/PF −9,1/0,9/8,3; σ 6,4; σ_h 4,0; ω 3,8; τ_esc 7,1; τ_prof 5,8. DIC 7835 × 7939; LOO aluno ΔELPD 41,1 (EP 7,6), 93/120, 4 alunos k > 0,7 (A013, A070, A092, A093), sem eles 45,3; LOO por nota 59/1200 k > 0,7. PPC: média 0,47, DP 0,56, mín 0,17, máx 1,00, individual 0,60, professor 0,45, escola 0,66. Coda: GR ≤ 1,01, ESS mín 545, MCSE ≤ 0,07, Heidelberger ok exceto γ_E2 em 1/4 cadeias, Geweke 5/52 (hn) e 10/52 (hc). P(|δ_PM| < 2) = 0,84. Shrinkage: fator teórico 0,68 (1 aluno) → 0,91 (5).
- **Tentativas (dados novos):** T1 9 div, R-hat 1,02, δ −5,5/4,3/12,0; T2 0 div (ESS σ_h 307); T3 0 div (ESS σ_h 854), half-Cauchy também converge; T4 final 1 div em 4.000. A cadeia travada com half-Cauchy (R-hat 5) só ocorreu na versão anterior dos dados; T4 mantida por robustez.
- **Apresentação:** 28 slides; slide novo "Histórico 60: censura" (fig14); verossimilhança com caso censurado; todos os números atualizados. **Resumo:** 5 páginas, método e resultados atualizados.
- Proposta de 1 página para a professora: `proposta_modelo.pdf` (enviada no chat).
