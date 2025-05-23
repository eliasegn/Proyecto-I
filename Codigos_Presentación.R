library(ggplot2)
library(sde)
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
  
  return(list(S = S_milstein, transiciones = transiciones_milstein))
}

##############################################################################

# Parámetros
mu_real <- 0.5
sigma_real <- 0.1
S0 <- 1
T <- 1
N <- 1000
delta <- T / N
tray <- 1000 # número de trayectorias
estimadores_trayectorias <- numeric(tray + 1)

estimador_mu <- function(S){
  I <- 0
  eme <- length(S) - 1
  deltaS <- diff(S)
  for(i in 1:eme){
    I <- I + 1/T * (1 / S[i]) * deltaS[i]
  }
  return(I)
}

for(i in 1:tray){
  St <- milsteinMBG(S0, mu_real, sigma_real, T, N)$S
  mu <- estimador_mu(St)
  estimadores_trayectorias[i] <- mu
}

# Estadísticas
mean(estimadores_trayectorias)
quantile(estimadores_trayectorias, probs = c(0.25, 0.75))
quantile(estimadores_trayectorias, probs = c(0.025, 0.975))

S = milsteinMBG(S0, mu_real, sigma_real, T, N)$S
eme <- length(S) - 1
deltaS <- diff(S)
I_vals <- numeric(eme)  # vector para guardar los estimadores en cada paso
I <- 0

for(i in 1:eme){
  I <- I + 1/T * (1 / S[i]) * deltaS[i]
  I_vals[i] <- I
}

data_completa <- data.frame(
  t = seq(1:length(I_vals)),
  integral = I_vals
)  

ggplot(data_completa, aes(x = t)) + 
  geom_line(aes(y=integral), color = 'steelblue', size = 1, alpha= 0.7) +
  geom_hline(yintercept = mu_real, linetype = 'dashed', color = 'blue') +
  labs(x = 'Tiempo', y= 'Estimadores', title = 'Evolución de mu con información completa') + 
  theme_minimal()
data_completa$integral[length(data_completa$integral)]

######## Variación Cuadrática ##############################

estimadores_sigma <- numeric(tray + 1)

estimarsigma <- function(S){
  riemann_approx <- sum(S^2 * delta)
  riemann_approx
  deltaS = diff(S)
  var_cuad = sqrt(sum(deltaS^2) / riemann_approx)
  return(var_cuad)
}

for(j in 1:length(estimadores_sigma)){
  S = milsteinMBG(S0, mu_real, sigma_real, T, N)$S
  sigma <- estimarsigma(S)
  estimadores_sigma[j] <- sigma
}

# Estadísticas
mean(estimadores_sigma)
quantile(estimadores_sigma, probs = c(0.25, 0.75))
quantile(estimadores_sigma, probs = c(0.025, 0.975))

S = milsteinMBG(S0, mu_real, sigma_real, T, N)$S
var_cuad_vals <- numeric(N - 2)

for (n in 2:(eme - 1)) {
  deltaS_n <- diff(S[1:(n + 1)])
  riemann_n <- sum(S[1:n]^2) * delta
  var_cuad_vals[n - 1] <- sqrt(sum(deltaS_n^2) / riemann_n)
}

# Creamos el data.frame
df_var_cuad <- data.frame(
  t = 1:(eme - 2),
  var_cuadr = var_cuad_vals
)

ggplot(df_var_cuad, aes(x = t)) +
  geom_line(aes(y=var_cuadr), color = "steelblue", size = 1) +
  geom_hline(yintercept = sigma_real, linetype = 'dashed', color = 'blue')+
  labs(x = "Tiempo", y = "Estimador de sigma",
       title = "Evolución del estimador de sigma por variación cuadrática") +
  theme_minimal()


############### ESTIMACIÓN POR TRANSICIONES #############################33

estimador_mu2 <- function(S){
  Y <- diff(log(S))
  
  # Calculamos los estimadores totales
  estimador_sigmac <- 1/(N*delta) * sum((Y - mean(Y))^2) 
  estimador_mu <- mean(Y) / delta + estimador_sigmac / 2
  estimador_sigma <- sqrt(estimador_sigmac)
  return(data.frame(mu= estimador_mu, sigma=estimador_sigma))
}

estimadores_mu3 <- numeric(tray)
estimadores_sigma3 <- numeric(tray)
for(j in 1:tray){
  S = milsteinMBG(S0, mu_real, sigma_real, T, N)$S
  muu <- estimador_mu2(S)$mu
  sigmaa <- estimador_mu2(S)$sigma
  estimadores_mu3[j] <- muu
  estimadores_sigma3[j] <- sigmaa
}

# Estadísticos para mu
mean(estimadores_mu3)
quantile(estimadores_mu3, probs = c(0.25, 0.75))
quantile(estimadores_mu3, probs = c(0.025, 0.975))

# Estadísticos para sigma
mean(estimadores_sigma3)
quantile(estimadores_sigma3, probs = c(0.25, 0.75))
quantile(estimadores_sigma3, probs = c(0.025, 0.975))

