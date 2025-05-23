library(ggplot2)

################## Ajuste a una Base de Datos ######################

# Definimos las funciones para calcular los estimadores
estimador_mu_integral <- function(S, diferencia) {
  deltaS <- diff(S)
  S_lag <- S[-length(S)]
  mu_hat <- sum((1 / S_lag) * deltaS) / (length(deltaS) * diferencia)
  return(mu_hat)
}

estimarsigma_vc <- function(S, diferencia){
  riemann_approx <- sum(S[-length(S)]^2) * diferencia
  deltaS = diff(S)
  var_cuad = sqrt(sum(deltaS^2) / riemann_approx)
  return(var_cuad)
}

estimador_transiciones <- function(S, diferencia){
  Y <- diff(log(S))
  N <- length(S)
  # Calculamos los estimadores totales
  estimador_sigmac <- 1/(N*diferencia) * sum((Y - mean(Y))^2) 
  estimador_mu <- mean(Y) / diferencia + estimador_sigmac / 2
  estimador_sigma <- sqrt(estimador_sigmac)
  return(data.frame(mu= estimador_mu, sigma=estimador_sigma))
}

# Descargamos la base de datos
df <- datoshistoricos

df <- df[-1, ]  # eliminamos la primera fila

df[ , 1:ncol(df)] <- lapply(df[ , 1:ncol(df)], as.numeric) # convertimos las
# columnas a flotantes

# Definimos nuestro proceso
St <- df$Close[1:252]

######## Estadísticas de la base de datos ###################
min(St)
quantile(St, probs = c(0.25, 0.5, 0.75))
max(St)
mean(St)
sd(St)
var(St)

#############################################################

# Definimos la delta del tiempo
delta_t <- 1/252

# Estimamos sigma y mu por información completa
sigma_hat_vc <- estimarsigma_vc(St, delta_t)
mu_hat_integral <- estimador_mu_integral(St, delta_t)

# Estimamos sigma y mu por transiciones
transiciones_hat <- estimador_transiciones(St, delta_t)
sigma_hat_transiciones <- transiciones_hat$sigma
mu_hat_transiciones <- transiciones_hat$mu

# Mostramos resultados
cat("Estimación de mu por integral estocástica:\n")
print(mu_hat_integral)

cat("\nEstimación de sigma por variación cuadrática:\n")
print(sigma_hat_vc)

cat("\nEstimaciones por log-retornos (transiciones):\n")
print(transiciones_hat)

################### Simulaciones con los parámetros estimados ##############

# La función para generar trayectorias
milsteinMBG <- function(S0, mu, sigma, T, N) {
  dt <- T / N  # Tamaño del paso
  t <- seq(0, T, length.out = N + 1)  # Vector de tiempo
  
  # Generamos el MBE
  dW <- rnorm(N, mean = 0, sd = sqrt(dt))
  W <- c(0, cumsum(dW))  # Proceso de Wiener
  transiciones_milstein <- dW  # Guardamos los incrementos
  
  # Método de Milstein
  S_milstein <- numeric(N + 1)
  S_milstein[1] <- S0  # Condición inicial
  
  for (i in 1:N) {
    S_milstein[i + 1] <- S_milstein[i] + mu * S_milstein[i] * dt + 
      sigma * S_milstein[i] * dW[i] + 
      0.5 * sigma^2 * S_milstein[i] * (dW[i]^2 - dt)
  }
  
  return(S_milstein)
}

set.seed(123)
S0 <- St[1]
iter <- length(St)
T <- 1
sim1_ic <- milsteinMBG(S0, mu_hat_integral, sigma_hat_vc, T, iter-1)
sim2_ic <- milsteinMBG(S0, mu_hat_integral, sigma_hat_vc, T, iter-1)
sim3_ic <- milsteinMBG(S0, mu_hat_integral, sigma_hat_vc, T, iter-1)
sim4_ic <- milsteinMBG(S0, mu_hat_integral, sigma_hat_vc, T, iter-1)
simulaciones_ic <- data.frame(t = 1:iter, St = St, s1 = sim1_ic, s2 = sim2_ic, s3 = sim3_ic, s4 = sim4_ic)
ggplot(simulaciones_ic, aes(x = t)) +
  geom_line(aes(y = St, color = "Datos Reales", linetype = "Datos Reales"), size = 1) +
  geom_line(aes(y = s1, color = "Simulación 1", linetype = "Simulación 1")) +
  geom_line(aes(y = s2, color = "Simulación 2", linetype = "Simulación 2")) +
  geom_line(aes(y = s3, color = "Simulación 3", linetype = "Simulación 3")) +
  geom_line(aes(y = s4, color = "Simulación 4", linetype = "Simulación 4")) +
  scale_color_manual(values = c("Datos Reales" = "red",
                                "Simulación 1" = "blue",
                                "Simulación 2" = "green",
                                "Simulación 3" = "yellow",
                                "Simulación 4" = "deeppink")) +
  scale_linetype_manual(values = c("Datos Reales" = "solid",
                                   "Simulación 1" = "dashed",
                                   "Simulación 2" = "dashed",
                                   "Simulación 3" = "dashed",
                                   "Simulación 4" = "dashed")) +
  labs(title = "Simulaciones del MBG",
       x = "Tiempo (días)",
       y = "Precio simulado",
       color = "Trayectoria",
       linetype = "Trayectoria") +
  theme_minimal()

# Caso con transiciones

sim1_t <- milsteinMBG(S0, mu_hat_transiciones, sigma_hat_transiciones, T, iter-1)
sim2_t <- milsteinMBG(S0, mu_hat_transiciones, sigma_hat_transiciones, T, iter-1)
sim3_t <- milsteinMBG(S0, mu_hat_transiciones, sigma_hat_transiciones, T, iter-1)
sim4_t <- milsteinMBG(S0, mu_hat_transiciones, sigma_hat_transiciones, T, iter-1)
simulaciones_t <- data.frame(t = 1:iter, St = St, s1 = sim1_t, s2 = sim2_t, s3 = sim3_t, s4 = sim4_t)
ggplot(simulaciones_t, aes(x = t)) +
  geom_line(aes(y = St, color = "Datos Reales", linetype = "Datos Reales"), size = 1) +
  geom_line(aes(y = s1, color = "Simulación 1", linetype = "Simulación 1")) +
  geom_line(aes(y = s2, color = "Simulación 2", linetype = "Simulación 2")) +
  geom_line(aes(y = s3, color = "Simulación 3", linetype = "Simulación 3")) +
  geom_line(aes(y = s4, color = "Simulación 4", linetype = "Simulación 4")) +
  scale_color_manual(values = c("Datos Reales" = "red",
                                "Simulación 1" = "blue",
                                "Simulación 2" = "green",
                                "Simulación 3" = "yellow",
                                "Simulación 4" = "deeppink")) +
  scale_linetype_manual(values = c("Datos Reales" = "solid",
                                   "Simulación 1" = "dashed",
                                   "Simulación 2" = "dashed",
                                   "Simulación 3" = "dashed",
                                   "Simulación 4" = "dashed")) +
  labs(title = "Simulaciones del MBG",
       x = "Tiempo (días)",
       y = "Precio simulado",
       color = "Trayectoria",
       linetype = "Trayectoria") +
  theme_minimal()



################### Validación del modelo ##################################

library(nortest)

# Caso con información completa
Yt1 <- log(St) / sigma_hat_vc
delta_Y <- diff(Yt1)
residuales <- (1 / sigma_hat_vc) * diff(log(St)) - (mu_hat_integral - 0.5 * sigma_hat_vc^2) * delta_t
qqnorm(residuales)
qqline(residuales)
hist(residuales)
lillie.test(residuales)
ad.test(residuales)
cvm.test(residuales)

# Caso con transiciones
Yt2 <- log(St) / sigma_hat_transiciones
delta_Y2 <- diff(Yt2)
residuales2 <- delta_Y2 - (mu_hat_transiciones / sigma_hat_transiciones - sigma_hat_transiciones / 2) * delta_t
qqnorm(residuales2)
qqline(residuales2)
hist(residuales2)
lillie.test(residuales2)
ad.test(residuales2)
cvm.test(residuales2)

# Validación del modelo del caso con información completa

set.seed(1)
n <- length(residuales)        # tamaño de muestra
B <- 1000                        # número de muestras simuladas 
ks_pvalues <- numeric(B)
muestras_normales <- matrix(NA, nrow = B, ncol = n)

for (b in 1:B) {
  muestra <- rnorm(n, mean = 0, sd = sqrt(delta_t)) # Varianza delta_t
  muestras_normales[b, ] <- muestra
  ks_pvalues[b] <- ks.test(residuales, muestra)$p.value
}

# Promedio de p-values
mean_pval <- mean(ks_pvalues)
cat("Promedio de p-values de las pruebas KS:", mean_pval, "\n")

# Hay información estadística suficiente para decir que los errores son normales

# Validación del modelo del caso con transiciones
Yt2 <- log(St) / sigma_hat_transiciones
delta_Y2 <- diff(Yt2)
residuales2 <- delta_Y2 - (mu_hat_transiciones / sigma_hat_transiciones - sigma_hat_transiciones / 2) * delta_t
ks_pvalues_transiciones <- numeric(B)
muestras_normales_transiciones <- matrix(NA, nrow = B, ncol = n)

for (b in 1:B) {
  muestra <- rnorm(n, mean = 0, sd = sqrt(delta_t)) # Varianza delta_t
  muestras_normales_transiciones[b, ] <- muestra
  ks_pvalues_transiciones[b] <- ks.test(residuales2, muestra)$p.value
}

# Promedio de p-values
mean_pval_transiciones <- mean(ks_pvalues_transiciones)
cat("Promedio de p-values de las pruebas KS:", mean_pval_transiciones, "\n")

