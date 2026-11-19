# ==============================================================
# Trabalho prático IV - Análise de Dados
# Regressão linear múltipla
# ==============================================================

# ---- 0. Pacotes ----------------------------------------------
#pacotes <- c("MASS", "car", "lmtest", "gamlss", "ggplot2")
#for (p in pacotes) {
 # if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
#}
library(MASS)     # stepAIC
library(car)      # vif
library(lmtest)   # bptest, dwtest
library(gamlss)   # stepGAIC
library(ggplot2)

# ---- 1. Importação e preparação dos dados --------------------
dados <- read.table("BasedeDados.txt", header = TRUE, 
                    sep = "\t", stringsAsFactors = FALSE)
str(dados)

pre  <- dados[dados$VISIT == " PRE",  ]
post <- dados[dados$VISIT == "POST", ]
pre  <- pre[order(pre$ID), ]
post <- post[order(post$ID), ]
stopifnot(all(pre$ID == post$ID))

# Variável dependente: variação absoluta (POST - PRE)
bd <- data.frame(
  ID         = pre$ID,
  VARIACAO   = post$WEIGHT - pre$WEIGHT,
  TREAT      = factor(pre$TREAT, levels = c("CONTROL", "PLACEBO", "TREATMENT")),
  AGE        = pre$AGE,
  LIVING     = factor(pre$LIVINGSTATUS),
  SMOKE      = factor(pre$SMOKESTATUS),
  GENDER     = factor(pre$GENDER)
)
bd$LIVING <- relevel(bd$LIVING, ref = "ALONE (SINGLE)")
bd$SMOKE  <- relevel(bd$SMOKE,  ref = "NON-SMOKER")
bd$GENDER <- relevel(bd$GENDER, ref = "FEMALE")

str(bd)
summary(bd)

# ---- 2. Análise exploratória ---------------------------------
summary(bd$VARIACAO)
sd(bd$VARIACAO)

hist(bd$VARIACAO, main = "Histograma da variação de peso",
     xlab = "Variação de peso (kg)", col = "lightblue", breaks = 15)

boxplot(VARIACAO ~ TREAT,  data = bd, main = "Variação por grupo",
        xlab = "Grupo", ylab = "Variação (kg)")
boxplot(VARIACAO ~ GENDER, data = bd, main = "Variação por sexo",
        xlab = "Sexo", ylab = "Variação (kg)")
boxplot(VARIACAO ~ SMOKE,  data = bd, main = "Variação por estado de fumador",
        xlab = "Estado de fumador", ylab = "Variação (kg)")
boxplot(VARIACAO ~ LIVING, data = bd, main = "Variação por situação de vida",
        xlab = "Situação de vida", ylab = "Variação (kg)", las = 2)

plot(bd$AGE, bd$VARIACAO, main = "Variação de peso vs idade",
     xlab = "Idade", ylab = "Variação (kg)", pch = 19)
abline(lm(VARIACAO ~ AGE, data = bd), col = "red")

# ---- 3. Modelo inicial (apenas efeitos principais) -----------
# Variáveis categóricas entram como dummy automaticamente (factor)
m0 <- lm(VARIACAO ~ TREAT + AGE + LIVING + SMOKE + GENDER, data = bd)
summary(m0)

# ---- 4. Modelo completo com todas as interações de 1.ª ordem -
m_cheio <- lm(VARIACAO ~ (TREAT + AGE + LIVING + SMOKE + GENDER)^2, data = bd)
summary(m_cheio)
anova(m_cheio)

# ---- 5. Seleção automática de variáveis/interações -----------
# 5.1 stepAIC (MASS) - stepwise nos dois sentidos
m_step <- stepAIC(m0,
                  scope = list(lower = ~ 1,
                               upper = ~ (TREAT + AGE + LIVING + SMOKE + GENDER)^2),
                  direction = "both", trace = TRUE)
summary(m_step)

# 5.2 stepGAIC (gamlss) - alternativa, com critério AIC (k = 2)
m_gamlss0 <- gamlss(VARIACAO ~ TREAT + AGE + LIVING + SMOKE + GENDER,
                    data = bd, family = NO, trace = FALSE)
m_gamlss  <- stepGAIC(m_gamlss0,
                      scope = list(lower = ~ 1,
                                   upper = ~ (TREAT + AGE + LIVING + SMOKE + GENDER)^2),
                      direction = "both", k = 2)
summary(m_gamlss)

# 5.3 Modelo final (parcimonioso) - o selecionado pelo stepAIC
m_final <- m_step
formula(m_final)

# ---- 6. Significância global e dos coeficientes --------------
summary(m_final)    # teste F global, R2, R2 ajustado, testes t
anova(m_final)      # ANOVA sequencial (tipo I)
Anova(m_final, type = "II")   # testes por termo (tipo II)
confint(m_final)    # IC a 95% para os coeficientes

# Comparação com o modelo só de efeitos principais
anova(m0, m_final)
AIC(m0, m_final, m_cheio)
BIC(m0, m_final, m_cheio)

# ---- 7. Multicolinearidade -----------------------------------
vif(m0)
tryCatch(print(vif(m_final)),
         error = function(e) message("VIF do modelo final não calculável: ",
                                     conditionMessage(e)))
vif(m_final, type = "predictor")
# ---- 8. Verificação dos pressupostos -------------------------
par(mfrow = c(2, 2))
plot(m_final)
par(mfrow = c(1, 1))

res  <- residuals(m_final)
rpad <- rstandard(m_final)

# Normalidade dos resíduos
shapiro.test(res)
qqnorm(rpad, main = "Q-Q plot dos resíduos padronizados"); qqline(rpad, col = "red")
hist(rpad, main = "Histograma dos resíduos padronizados",
     xlab = "Resíduos padronizados", col = "lightgreen", breaks = 15)

# Homocedasticidade
bptest(m_final)                  # Breusch-Pagan
ncvTest(m_final)                 # teste da variância não constante

# Independência dos erros
dwtest(m_final)                  # Durbin-Watson

# Linearidade (resíduos vs ajustados)
plot(fitted(m_final), res, pch = 19,
     main = "Resíduos vs valores ajustados",
     xlab = "Valores ajustados", ylab = "Resíduos")
abline(h = 0, col = "red")

# ---- 9. Observações influentes e outliers --------------------
cook <- cooks.distance(m_final)
plot(cook, type = "h", main = "Distância de Cook", ylab = "Distância de Cook")
abline(h = 4 / nrow(bd), col = "red", lty = 2)
which(cook > 4 / nrow(bd))

hat <- hatvalues(m_final)
which(hat > 2 * length(coef(m_final)) / nrow(bd))   # leverage elevado

which(abs(rpad) > 3)                                 # outliers
outlierTest(m_final)
influencePlot(m_final)


# ---- 10. Validação do modelo (validação cruzada 10-fold) -----
set.seed(123)
k      <- 10
folds  <- sample(rep(1:k, length.out = nrow(bd)))
erro   <- numeric(k)
for (i in 1:k) {
  treino <- bd[folds != i, ]
  teste  <- bd[folds == i, ]
  mod    <- lm(formula(m_final), data = treino)
  pred   <- predict(mod, newdata = teste)
  erro[i] <- sqrt(mean((teste$VARIACAO - pred)^2))
}
cat("RMSE médio (validação cruzada):", mean(erro), "\n")

# ---- 11. Resultados finais -----------------------------------
cat("\nEquação do modelo final:\n")
print(formula(m_final))
print(round(coef(summary(m_final)), 4))
cat("R2 =", round(summary(m_final)$r.squared, 4),
    "| R2 ajustado =", round(summary(m_final)$adj.r.squared, 4), "\n")

