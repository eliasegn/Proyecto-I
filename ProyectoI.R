library(ggplot2)
library(sde)

x <- 1
r <- 0.2
sigma <- 1
T <- 1
N <- 1000
delta <- T/N

set.seed(123)
X <- GBM(x=1, r=0.2, sigma=1, T=1, N=1000)
Y <- diff(log(X))

plot(X)

estimador_sigmac <- 1/(N*delta) * sum((Y - mean(Y))^2) 
estimador_sigmac

estimador_mu <- mean(Y) / delta + estimador_sigmac / 2

##############

eulerMBG <- function(S0, mu, sigma, T, N) {
  # S0     Precio inicial del activo
  # mu     Tasa de retorno esperada
  # sigma  Volatilidad del activo
  # T      Tiempo total de simulación (en años)
  # N      Número de pasos en la simulación
  
  dt <- T / N                            # Tamaño del paso de tiempo
  t <- seq(0, T, length.out = N + 1)     # Vector de tiempo
  
  # Generamos un MBE
  W <- rnorm(N, mean = 0, sd = sqrt(dt))
  W <- c(0, cumsum(W))                   # Insertamos W[0] = 0 al inicio
  
  # Aplicamos método de Euler
  S_euler <- numeric(N + 1)
  S_euler[1] <- S0                       # Condición inicial
  transiciones_euler <- numeric(N)
  
  for (i in 1:N) {
    dW <- W[i + 1] - W[i]
    transiciones_euler[i] <- dW
    S_euler[i + 1] <- S_euler[i] + mu * S_euler[i] * dt + sigma * S_euler[i] * dW
  }
  
  return(list(S_euler = S_euler, transiciones_euler = transiciones_euler))
}

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

x0 <- 1
m <- 0.2
s <- 1
Ti <- 1
N <- 1000
delta <- T/N

prueba1 <- milsteinMBG(x0, m, s, Ti, N)
plot(prueba1$S, type='l')
prueba2 <- GBM(x0, m, s, Ti, N)

procesos <- 100

libreria <- numeric()
nolibreria <- numeric()
for(i in 1:procesos){
  libreria <- c(libreria, list(GBM(x0, m, s, Ti, N)))
}

for(i in 1:procesos){
  nolibreria <- c(nolibreria, list(milsteinMBG(x0, m, s, Ti, N)$S))
}

quinlib <- numeric()
for(j in 1:length(libreria)){
  quinlib <- c(quinlib, libreria[[j]][500])
}

quinnolib <- numeric()
for(j in 1:lenth(nolibreria)){
  quinnolib <- c(quinnolib, nolibreria[[j]][500])
}

ks.test(quinlib, quinnolib)

# jaja, pues sí son iguales



############# EXPOSICIÓN ########################
# SIMULACIÓN Y PARÁMETROS DEL MBG

set.seed(103)
# Parámetros reales
mu_real <- 0.5
sigma_real <- 0.1
S0 <- 1
T <- 1
N <- 10000
delta <- T / N

# Simulamos el Browniano Geométrico
S <- milsteinMBG(S0, mu_real, sigma_real, T, N)$S

#############################################################################

grafico <- data.frame(t = seq(1:length(S)),
                      S_sim = S)

# Graficamos nuestro MBG
ggplot(grafico, aes(x = t)) + 
  geom_line(aes(y = S_sim), color = "black", size = 1, linetype = "solid") + 
  #geom_line(aes(y = S_exp), color = "red", size = 1, linetype = "dashed") +
  labs(y = "Valor", title = "Simulación") +
  theme_minimal()

# Diferencias de logaritmos
Y <- diff(log(S))

# Calculamos los estimadores totales
estimador_sigmac <- 1/(N*delta) * sum((Y - mean(Y))^2) 
estimador_mu <- mean(Y) / delta + estimador_sigmac / 2
estimador_mu
estimador_sigmac
estimador_sigma <- sqrt(estimador_sigmac)
estimador_sigma

# Inicializamos vectores para los estimadores
estim_mu <- numeric(length(Y))
estim_sigmac <- numeric(length(Y))

# Calculamos los estimadores a lo largo del tiempo
for (i in 10:length(Y)) {
  Yi <- Y[1:i]
  estim_sigmac[i] <- sqrt( 1 / (N* delta) * sum((Yi - mean(Yi))^2) )
  estim_mu[i] <- mean(Yi) / delta + estim_sigmac[i] / 2
}

# Graficamos
df_estim <- data.frame(
  t = 1:N,
  mu = estim_mu,
  sigma = (estim_sigmac)^2
)

ggplot(df_estim, aes(x = t * delta)) +
  geom_line(aes(y = mu), color = "#FF3030", size = 1, alpha = 0.7) +
  geom_hline(yintercept = mu_real, linetype = "dashed", color = "#8B1A1A") +
  geom_line(aes(y = sigma), color = "#1E90FF", size = 1, alpha = 0.7) +
  geom_hline(yintercept = sigma_real, linetype = "dashed", color = "dodgerblue4") +
  labs(y = "Estimadores", x = "Tiempo",
       title = "Evolución de estimadores de mu (rojo) y sigma (azul)") +
  coord_cartesian(ylim = c(-10 , 10)) +
  theme_minimal()

############### Información Completa con Integrales Estocásticas ##############3

I <- 0
M <- length(S) - 1
deltaS <- diff(S)
for(i in 1:M){
  I <- I + 1/T * (1 / S[i]) * deltaS[i]
}
I

#¡Es un buen estimador!

M <- length(S) - 1
deltaS <- diff(S)
I_vals <- numeric(M)  # vector para guardar los estimadores en cada paso
I <- 0

for(i in 1:M){
  I <- I + 1/T * (1 / S[i]) * deltaS[i]
  I_vals[i] <- I
}

data_completa <- data.frame(
  t = seq(1:length(I_vals)),
  integral = I_vals
)  

ggplot(data_completa, aes(x = t)) + 
  geom_line(aes(y=integral), color = 'coral4', size = 1, alpha= 0.7) +
  geom_hline(yintercept = mu_real, linetype = 'dashed', color = 'coral') +
  labs(x = 'Tiempo', y= 'Estimadores', title = 'Evolución de estimador con información completa') + 
  theme_minimal()
data_completa$integral[length(data_completa$integral)]

# Ahora intentamos hacer una estimación para sigma con variación cuadrática

# Suponiendo que ya tenemos S y T
M <- length(S)
dt <- T / M


int=function(path,delta)
{
  n=length(path)-1
  int=numeric(n-1)
  for(i in 1:n)
  {
    int[i]=(path[i])*(delta)
  }
  
  int=cumsum(int)
  return(int)
}

# Aproximaciones
riemann_approx <- sum(S^2 * dt)
riemann_approx
deltaS = diff(S)
var_cuad = sqrt(sum(deltaS^2) / riemann_approx)
var_cuad

# Es mejor estimador

var_cuad_vals <- numeric(N - 2)

for (n in 2:(M - 1)) {
  deltaS_n <- diff(S[1:(n + 1)])
  riemann_n <- sum(S[1:n]^2) * dt
  var_cuad_vals[n - 1] <- sqrt(sum(deltaS_n^2) / riemann_n)
}

# Creamos el data.frame
df_var_cuad <- data.frame(
  t = 1:(M - 2),
  var_cuadr = var_cuad_vals
)

ggplot(df_var_cuad, aes(x = t)) +
  geom_line(aes(y=var_cuadr), color = "steelblue", size = 1) +
  geom_hline(yintercept = sigma_real, linetype = 'dashed', color = 'blue')+
  labs(x = "Tiempo", y = "Estimador de sigma",
       title = "Evolución del estimador de sigma por variación cuadrática") +
  theme_minimal()

