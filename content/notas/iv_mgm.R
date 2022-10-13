# Espacio de trabajo ----
rm(list = ls())
options(scipen=999) # Prevenir notación científica

library(tidyverse)
library(janitor)
library(sandwich)
library(clubSandwich)
library(lfe)
library(AER)
library(gmm)

##Datos----

data.ingresos <- read_csv("ingresos_iv.csv",
                          locale = locale(encoding = "latin1"))

#Modelo exactamente identificado----

#Para tener una referencia, veamos lo que obtenemos con ivreg
#Nuestro modelo tiene cinco regresores más una constante
iv_ei <- ivreg(lwage ~ educ + exper + expersq + black + south |
                 . - educ + nearc4, data = data.ingresos)

stargazer(iv_ei,
          type="text",
          digits = 4)


#Estimadores con matrices
  
#Construimos las matrices X, Y y Z
data.ingresos <- data.ingresos %>% 
  mutate(constant=1)

#Cinco regresores y la constante
X <- data.matrix(select(data.ingresos, constant, educ, exper, expersq, black,
              south),
       rownames.force = T)

Y <- data.matrix(select(data.ingresos,lwage),
       rownames.force = T)

#Se excluye educación y se incluye en su posición nearc
Z <- data.matrix(select(data.ingresos, constant, nearc4, exper, expersq, black,
              south),
       rownames.force = T)
      
      
N <- nrow(X)
k <- ncol(X) # incluyendo la constante
      
#Estimamos beta
b <- solve(t(Z) %*% X) %*% t(Z) %*% Y
      
#Obtenemos el vector de coeficientes
b

#La matriz de varianzas, asumiendo homocedasticidad
u_hat <- Y-X%*%b

sigma2 <- as.numeric((1/N)*t(u_hat)%*%u_hat)

#Matriz de proyección
P <- Z%*%(solve(t(Z)%*%Z))%*%t(Z)

#Matriz de varianzas (con correción por muestras finitas (R por default multiplica por N/n-k)
V=sigma2*solve(t(X)%*%P%*%X)*(N/(N-k))

#Comparamos el coeficiente y el de educación con lo obtenido con ivreg
sqrt(diag(V))
stargazer(iv_ei,
          type="text",
          digits = 4)


#Modelo exactamente identificado con posible heterocedasticidad----
#En este caso los errores son más pequeños. Con VI pueden suceder cosas raras
stargazer(iv_ei, iv_ei, iv_ei,
          type="text",
          se = list(NULL,
                    sqrt(diag(vcovHC(iv_ei, type = "const"))),
                    sqrt(diag(vcovHC(iv_ei, type = "HC0")))),
          digits = 4,
          column.labels = c("VI, hom., ivreg",
                            "VI, hom., ivreg",
                            "VI, het., ivreg"))
    
#Con matrices
D <- diag(as.vector((Y-X%*%b)^2))
S_hat <- (1/(N)) * t(Z) %*% D %*% Z 

#Noten que HC0 no hace corrección por muestras pequeñas
Vr= N*solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)%*%(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%S_hat%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)%*%solve(t(X)%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X)

sqrt(diag(Vr))



#Modelo sobreidentificado con homocedasticidad----
iv_si <- ivreg(lwage ~ educ + exper + expersq + black + south |
                 . - educ + nearc4 + nearc2, data = data.ingresos)

stargazer(iv_si,
          type="text",
          digits = 4)


#Se excluye educación y se incluye en su posición nearc
Z <- data.matrix(select(data.ingresos, constant, nearc4, nearc2, exper, expersq, black,
                        south),
                 rownames.force = T)

P <- Z%*%(solve(t(Z)%*%Z))%*%t(Z)

b <- solve(t(X)%*%P%*%X) %*% t(X)%*%P%*%Y
b

#La matriz de varianzas es la misma que en el caso exactamente identificado
u_hat <- Y-X%*%b
sigma2 <- as.numeric((1/N)*t(u_hat)%*%u_hat)

#Noten que R hace correción de muestras finitas
V=sigma2*solve(t(X)%*%P%*%X)*(N/(N-k))
sqrt(diag(V))

stargazer(iv_si,
          digits = 4,
          type = 'text')




#Estimador de MGM óptimo----

#Usamos el paquete gmm

gmm_opt <- gmm(lwage ~ educ + exper + expersq + black + south,
               ~ nearc4 + nearc2 + exper + expersq + black + south,
               vcov = "HAC",
               wmatrix = "optimal",
               type = "twoStep",
               data = data.ingresos)

#Con matrices

#Primer paso: estimar beta con una matriz W subóptima
r <- k -1 + 2 # 1 endógena y 2 instrumentos
I <- data.matrix(diag(r))

#b1 <- solve(t(X)%*%P%*%X) %*% t(X)%*%P%*%Y
b1 <- solve(t(X)%*%Z%*%I%*%t(Z)%*%X)%*%t(X)%*%Z%*%I%*%t(Z)%*%Y

#Usemos b del punto anterior para obtener S_hat
D <- diag(as.vector((Y-X%*%b1)^2))
S_hat <- (1/N) * t(Z) %*% D %*% Z 






#Obtenemos bo (óptimo)
bo <- solve(t(X)%*%Z%*%solve(S_hat)%*%t(Z)%*%X)%*%
  t(X)%*%Z%*%solve(S_hat)%*%t(Z)%*%Y
bo

D <- diag(as.vector((Y-X%*%bo)^2))
S_tilde <- (1/N) * t(Z) %*% D %*% Z 

Vr <- (N)*solve(t(X) %*% Z %*% solve(S_tilde) %*% t(Z) %*% X)
sqrt(diag(Vr))

stargazer(gmm_opt,
          type="text",
          digits = 4)




##IV es el estimador de GMM para cualquier W----

gmm_iv_opt <- gmm(lwage ~ educ + exper + expersq + black + south,
                  ~ nearc4 + exper + expersq + black + south,
                  vcov = "iid",
                  wmatrix = "optimal",
                  type = "twoStep",
                  data = data.ingresos)

gmm_iv_ident <- gmm(lwage ~ educ + exper + expersq + black + south,
                    ~ nearc4 + exper + expersq + black + south,
                    vcov = "iid",
                    wmatrix = "ident",
                    type = "twoStep",
                    data = data.ingresos)

stargazer(gmm_iv_opt, gmm_iv_ident,
          type="text",
          digits = 4)

#Solo un instrumento
Z <- data.matrix(select(data.ingresos, constant, nearc4, exper, expersq, black,
                        south),
                 rownames.force = T)

#Estimamos beta
b <- solve(t(Z) %*% X) %*% t(Z) %*% Y
b
