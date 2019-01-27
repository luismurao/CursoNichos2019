## ----eval=TRUE, echo=FALSE, cache=FALSE,  message=FALSE, warning=FALSE, fig.align='center',fig.width=8----
library(dygraphs)
data_pkgs <- read.csv("Arichivos_01Intro/paquetesTS.csv")
pkgs <- data_pkgs[,2:3]
row.names(pkgs) <- as.Date(pkgs[,2])
dygraph(pkgs) %>% dyRangeSelector() 
library(knitr)


## ------------------------------------------------------------------------
set.seed(112)

## ------------------------------------------------------------------------
libro <- "En las montañas de la locura"
# Existe en la memoria de R el objeto libro?
exists("libro")
# y que hay sobre el objeto Libro?
exists("Libro")


## ------------------------------------------------------------------------
# Asigna a uno el valor numérico 1 usando el comando "<-" 
uno <- 1
# Alternativamente 
1 -> uno
# Asigna a dos el valor numérico 2 usando el comando "=" 
dos = 2
# Asignar el valor 3 a la variable 3 usando el comando assing
assign("tres", 3)
print(uno)
print(dos)
print(tres)


## ------------------------------------------------------------------------
Olor <- "Dulce"
olor <- "Agradable"
# Esto produce un error
# 0lor <- "Fétido" # la variable comienza con el digito cero
print(Olor)
print(olor)


## ------------------------------------------------------------------------
ol0r <- "fresco"
o1or <- "fresco"
print(ol0r)
print(o1or)

## ------------------------------------------------------------------------
.ayuda <- "SOS"
print(.ayuda)

## ------------------------------------------------------------------------
5 + 89

## ------------------------------------------------------------------------
18 - 167

## ------------------------------------------------------------------------
2 * 98

## ------------------------------------------------------------------------
6/5

## ------------------------------------------------------------------------
9 ^ 9

## ------------------------------------------------------------------------
9 %% 7

## ------------------------------------------------------------------------
9 %/% 7

## ------------------------------------------------------------------------
# Residuo + 7 * la parte entera de la fracción
(9 %% 7) + 7 * (9 %/% 7) 

## ------------------------------------------------------------------------
(1 + 1/300)^300

## ------------------------------------------------------------------------
exp(1)
sin(1)
cos(1)
log(1)
sqrt(1)

## ------------------------------------------------------------------------
pi

## ------------------------------------------------------------------------
options(digits = 16)
pi

## ------------------------------------------------------------------------
# Redondea dos digitos de precisión
round(10/3.4, 2)
# Entero mayor más cercano
ceiling(1/3)
# El entero menor más cercano
floor(1/3)

## ------------------------------------------------------------------------
# Auto Audi RS5 4.2 FSI 450 hp precio de lista
audi_rs5 <- 1484900.00
# Tasa anual de 12.9%
tasa_anual <- 0.129
# Cuanto pagaremos de interes en un año
(interes_anual <- audi_rs5*tasa_anual)
# Y en los cuatro años
(interes_4 <- audi_rs5*(tasa_anual*4))
# Total a pagar 
(total_audi_rs5 <- audi_rs5 + interes_4)

## ------------------------------------------------------------------------
# Generemos el vector f (el del ejemplo anterior)
f <- c(5,2,4,3)

## ------------------------------------------------------------------------
f[3]

## ------------------------------------------------------------------------
f[1]

## ------------------------------------------------------------------------
(f1 <- f[c(1,3,4)])

## ------------------------------------------------------------------------
1:10

## ------------------------------------------------------------------------
seq(from =0,to = 10,by = 2)

## ------------------------------------------------------------------------
seq(0,1,by=0.01)

## ------------------------------------------------------------------------
# longitud de la secuencia que queremos generar
tamano <- 40
seq(13,25,length.out=tamano)

## ------------------------------------------------------------------------
vec <- c("A","B")
rep(vec, times=2)
rep(vec,each=2)

## ------------------------------------------------------------------------
vec <- 1:1000
sample(x = vec,size = 100,replace = FALSE)


## ------------------------------------------------------------------------
l <- seq(0,1,by=0.0025)
print(l)
# La longitud de l
length(l)

## ------------------------------------------------------------------------
elementos <- seq(1:16)
mat <- matrix(elementos,nrow = 4,ncol = 4)
print(mat)

## ------------------------------------------------------------------------
mat_reng <- matrix(elementos,nrow = 4,ncol = 4,byrow = T)
print(mat_reng)

## ------------------------------------------------------------------------
# Matriz tipo cadena
letras <- LETTERS[1:9]
mat_carac <- matrix(letras,ncol=3,nrow = 3,byrow = F)
print(mat_carac)

## ------------------------------------------------------------------------
letras_b <- LETTERS[1:8]
numero_b <- 1:8
letras_numeros <- c(letras_b,numero_b)
mat_num_car <- matrix(letras_numeros,ncol=4,nrow = 4)
print(mat_num_car)

## ------------------------------------------------------------------------
dats <- seq(2,4,length.out = 48)
M <- matrix(dats,nrow = 6,ncol = 8,byrow = TRUE)
# Extraigo el elemento a_{3,6}
print(M[3,6])

## ------------------------------------------------------------------------
print(M[2, ])

## ------------------------------------------------------------------------
print(M[,3])

## ------------------------------------------------------------------------
d_mat <- 1:25
M <- matrix(d_mat,ncol = 5,nrow = 5, byrow = TRUE)
print(M[1:3,c(1,2)])

## ------------------------------------------------------------------------
print(M[c(1,2,3,4,5),c(1,5)])
print(M[1:5,c(1,5)])

## ------------------------------------------------------------------------
pollito_id <- 1:75
lote <- paste0("L",rep(1:5,each=15))
peso_inicial <- runif(75,min = 8,max = 12)
pollitos <- data.frame(pollito_id,lote,peso_inicial)

## ----eval=F--------------------------------------------------------------
## print(pollitos[c(1:3,21:23,31:33,51:53,71:73),])

## ----eval=T,echo=F-------------------------------------------------------
kable(pollitos[c(1:3,21:23,31:33,51:53,71:73),])

## ------------------------------------------------------------------------
dieta <- paste0("D",rep(1:5,each=15))
pollitos[,"Dieta"] <- dieta

## ------------------------------------------------------------------------
peso_d1 <- rnorm(mean = 2250, n = 15,sd = 1)
peso_d2 <- rnorm(mean = 2600, n = 15,sd = 1.4)
peso_d3 <- rnorm(mean = 1800, n = 15,sd = 3)
peso_d4 <- rnorm(mean = 2100, n = 15,sd = .6)
peso_d5 <- rnorm(mean = 2000, n = 15,sd = 1.7)
pesos_dietas <- c(peso_d1,peso_d2,peso_d3,peso_d4,peso_d5)
pollitos$peso_dieta <- pesos_dietas 

## ----eval=FALSE----------------------------------------------------------
## # Imprimamos una parte del data.frame
## print(pollitos[c(1:3,21:23,31:33,51:53,71:73),])

## ----results='asis',echo=FALSE-------------------------------------------
library(knitr)
kable(pollitos[c(1:3,21:23,31:33,51:53,71:73),])

print(pollitos[c(1:3,21:23,31:33,51:53,71:73),c(1,4)])

dim(pollitos)
## ------------------------------------------------------------------------
pollitos[1:10,c("Dieta")]
names(pollitos)

pollitos[5:75,c("peso_inicial","Dieta","peso_dieta")]

## ------------------------------------------------------------------------
pollitos$Dieta
pollitos$Dieta[1:10]

## ------------------------------------------------------------------------
# un vector de numeros aleatorio del 0 a la 1
vector <- runif(10)
# Un subconjunto del DF pollitos
data_frame_pollito <- pollitos[1:10,]
# uan lista
sublista <- list(a=1:10)

lista <- list(vector=vector, pollitos_df = data_frame_pollito, sublista=sublista)
print(lista)

## ------------------------------------------------------------------------
# El segundo elemento de la lista es el DF pollitos
lista[[2]]

## ------------------------------------------------------------------------
lista$pollitos_df

## ------------------------------------------------------------------------
lista[["pollitos_df"]]

## ------------------------------------------------------------------------
# El elemento de lista que a la sublista que contiene al vector objetivo es el 3
# Esta sublista esta contituida por un solo elemento el cual es el vector a
# El elemento del vector a que nos interesa es el 10
indice_nivel_1 <- 2
indice_nivel_2 <- 1
indice_nivel_3 <- 5
lista[[indice_nivel_1]][[indice_nivel_2]][[indice_nivel_3]]

## ------------------------------------------------------------------------
lista[["sublista"]][["a"]][[indice_nivel_3]]


## ------------------------------------------------------------------------
lista[["nuevo_elem"]] <- rnorm(30)
lista

## ------------------------------------------------------------------------
lista$nuevo_elem_2 <- rbinom(30,size = 1,prob = 0.5)
lista

setwd("~/Dropbox/CursoNichos2019/")

## ------------------------------------------------------------------------
df_csv <- read.csv("Arichivos_01Intro/data_Clean.csv")
head(df_csv)

## ------------------------------------------------------------------------
df_tab <- read.table("Arichivos_01Intro/ambys_data_dynamic_map.txt",header = T)
head(df_tab[,1:4])

## ------------------------------------------------------------------------
df_tab_to_csv <- read.csv("Arichivos_01Intro/ambys_data_dynamic_map.txt",header = T,sep=" ")
head(df_tab_to_csv[,1:4])

## ------------------------------------------------------------------------
d <- read.csv("Arichivos_01Intro/ambys_data_dynamic_map.txt",header = T,sep=" ",fileEncoding = "UTF-8")
head(d[,5:7])

## ------------------------------------------------------------------------
# install.packages("gdata")
library(gdata)
d_xls <- read.xls("Arichivos_01Intro/data_Clean.xlsx",sheet = 1)
head(d_xls)
d_occs <- read.xls("Arichivos_01Intro/data_Clean.xlsx",sheet = 2)


## ------------------------------------------------------------------------
library(XLConnect)
lib_xls <- loadWorkbook("Arichivos_01Intro/data_Clean.xlsx")
# Mostrar las hojas dispobles
getSheets(lib_xls)
# Cargar una hoja
d_xlsC <- readWorksheet(object = lib_xls,sheet = "ROC")
head(d_xlsC)
d_xlsOCCs <- readWorksheet(object = lib_xls,sheet = "OCCS")
head(d_xlsOCCs)

## ----eval=FALSE----------------------------------------------------------
## if (condicion) {
##   comandos o operaciones cuando TRUE
##   }
## if (condicion) {
##   comandos o operaciones cuando TRUE
## } else {
##     comandos o operaciones cuando FALSE
##   }
## 

## ------------------------------------------------------------------------
n1 <- 1
n2 <- 2
n3 <- 3
n4 <- 4

if(n1 %% 2 == 0){
  cat(n1,"Es par")
} else {
  cat(n1, "No es par")
}

if(n2 %% 2 == 0){
  cat(n2,"Es par")
} else {
  cat(n2, "No es par")
}


if(n3 %% 2 == 0){
  cat(n3,"Es par")
} else {
  cat(n3, "No es par")
}

if(n4 %% 2 == 0){
  cat(n4,"Es par")
} else {
  cat(n4, "No es par")
}



## ----eval=FALSE----------------------------------------------------------
## for (iterador in vector) {comandos o operaciones}
## 

## ------------------------------------------------------------------------
# Genero una secuencia de números enteros (del 150 al 550) aleatoria
vec_num <- sample(150:550,size = 100,replace = F)
# Nuestro vector es 
print(vec_num)

for(n in vec_num){
  if(n %% 2 == 0)
    cat(n,"Es par\n") 
  else 
    cat(n, "No es par\n")
}


## ----eval=FALSE----------------------------------------------------------
## # Iteramos sobre los elementos de vec_num
 vec_size <- length(vec_num)
 for(iterador in 1:vec_size){
   n <- vec_num[iterador]
    if(n %% 2 == 0)
     cat(n,"Es par\n")
   else
     cat(n, "No es par\n")
 }

## ------------------------------------------------------------------------
vec_pares <- NULL
vec_impares <- NULL
# Iteramos sobre los elementos de vec_num
vec_size <- length(vec_num)

for(iterador in 1:vec_size){
  n <- vec_num[iterador]
   if(n %% 2 == 0){
     vec_pares[iterador] <-  n
   }
  else{
    vec_impares[iterador] <- n
  } 
}

## ------------------------------------------------------------------------
vec_pares_limpio <- as.vector(na.omit(vec_pares))
print(vec_pares_limpio)
vec_impares_limpio <-as.vector( na.omit(vec_impares))
print(vec_impares_limpio)

## ------------------------------------------------------------------------
vec_pares_limpio <- vec_pares[which(!is.na(vec_pares))]
vec_impares_limpio <- vec_impares[which(!is.na(vec_impares))]

## ------------------------------------------------------------------------
cat("El numero de enteros pares en vec_num fue de",length(vec_pares_limpio),"\n")
cat("El numero de enteros impares en vec_num fue de",length(vec_impares_limpio),"\n")

## ------------------------------------------------------------------------
# Con n=2
fib_n_menos2 <- 0
fib_n_menos1 <- 1
fib_n <- fib_n_menos1 + fib_n_menos2 
fibonacci <- c(fib_n_menos2,fib_n_menos1)


while(fib_n < 400){
  
  fibonacci <- c(fibonacci,fib_n)
  fib_n_menos2 <- fib_n_menos1
  fib_n_menos1  <- fib_n
  fib_n <- fib_n_menos2 + fib_n_menos1

}

print(fibonacci)

## ------------------------------------------------------------------------
print(length(fibonacci))

## ------------------------------------------------------------------------
fib_n_menos2 <- 0
fib_n_menos1 <- 1
fib_n <- fib_n_menos1 + fib_n_menos2 
n=3

while(fib_n < 4000000000000){
  assign(paste0("fibonacci_",n),value = fib_n) 
  fib_n_menos2 <- fib_n_menos1
  fib_n_menos1  <- fib_n
  n <- n+1
  fib_n <- fib_n_menos2 + fib_n_menos1
}

print(c(fibonacci_3,fibonacci_4,fibonacci_5,fibonacci_6,fibonacci_9,fibonacci_15))

## ------------------------------------------------------------------------
es_par <- function(x){
  # Cheamos que x es entero
  if(x%%1==0){
    # Checamos si es par
    if(x%%2==0){
      return(TRUE)
    }
    else
      return(FALSE)
  }
  # Mensaje notificando que x no es entero
  warning("El numero tiene que ser un numero entero")
}

## ------------------------------------------------------------------------
# Primero con un numero que no es entero
es_par(2.2)
# Ahora con un impar
es_par(3)
# Finalmente con un par
es_par(18)

## ------------------------------------------------------------------------
occs_bd <- read.csv("./Arichivos_01Intro/occs_rep_amf.csv")
# Dimensiones del DF.
dim(occs_bd )
# Una explo de la BD (mostrando solo 6 registros)
head(occs_bd,6)
# Los nombres de las variables de la BD.
names(occs_bd)
# Los nombres de las especies
print(unique(occs_bd$name))

## ------------------------------------------------------------------------
# Creamos una lista vacia
occs_df_list <- list()
# Guardamos los nombres de las especies en un vector
sp_names <- unique(occs_bd$name)

## ------------------------------------------------------------------------
# Vector donde guardare los indices (vacio)
sp_index <- c()
# Iteraremos 
tiempo_novec <- system.time({
  # Recorremos cada especie en sp_names
  for(spi in 1:length(sp_names)){
    # Especie i en la iteracion i
    sp_name_i <- sp_names[spi]
    # Recorremos cada registro en BD 
  for(registro in 1:dim(occs_bd)[1]){
    # Preguntamos si coincide el nombre sp_name_i con el registro occs_bd$name (condicion 1)
    if(sp_name_i == occs_bd$name[registro]){
      sp_index <- c(sp_index,registro)
    }
  }
  # Los registros que cumplieron la condicion 1
  # en sp_index son ocupados para hacer el subset
  # Guardamos estos en la lista
  occs_df_list[[spi]] <- occs_bd[sp_index,]
  # reinializamos el vector sp_index
  sp_index <- c()
  }
  # Le ponemos nombre de la sp a cada elemento de la lista
  names(occs_df_list) <- sp_names
})

print(head(occs_df_list$`Pseudacris regilla`))
dim(occs_df_list$`Sceloporus occidentalis`)

## ------------------------------------------------------------------------
print(tiempo_novec)

## ------------------------------------------------------------------------
# Creamos una lista vacia
occs_df_list_vec <- list()
# nuevamente definamos una variable con los nombres de sps
sp_names <- unique(occs_bd$name)


## ------------------------------------------------------------------------
tiempo_vec <- system.time({
  for(spi in 1:length(sp_names)){
    # Especie i en la iteracion i
    sp_name_i <- sp_names[spi]
    sp_index_vec <- which(sp_name_i == occs_bd$name)
    occs_df_list_vec[[spi]] <- occs_bd[sp_index_vec,]
  }
})
names(occs_df_list_vec) <- sp_names
dim(occs_df_list_vec$`Sceloporus occidentalis`)

## ------------------------------------------------------------------------
print(tiempo_vec)

## ------------------------------------------------------------------------
tiempo_ratio <- tiempo_novec[[3]]/tiempo_vec[[3]]
print(tiempo_ratio)

## ------------------------------------------------------------------------
# Funcion de busqueda
# sp_name: Es el nombre de la especie
# bd_datos: data.frame con los datos de presencia
# col_name: Columna en el data.frame sobre la que se hara la busqueda
find_sp <- function(sp_name, bd_datos,col_name){
  # Preguntamos si coincide el nombre sp_name con el registro occs_bd[,col_name] (condicion 1)
  sp_index <- which(bd_datos[,col_name] == sp_name)
  # Subset con los datos que cumplen la condicion 1
  sp_df <- bd_datos[sp_index,]
  return(sp_df)
}

## ------------------------------------------------------------------------
sp_names <- unique(occs_bd$name)
tiempo_lapply <- system.time({
  occs_df_list_apply <- lapply(sp_names, 
                               function(sp) 
                                 find_sp(sp,bd_datos = occs_bd,col_name = "name") )
  
})


## ------------------------------------------------------------------------



tiempo_split <- system.time({
  occs_bdList <- occs_bd %>% split(.$name)
})

## ------------------------------------------------------------------------
cat("Tiempo de ejecucion:", tiempo_novec[[3]],"algoritmo no vectorizado\n")
cat("Tiempo de ejecucion:", tiempo_vec[[3]],"algoritmo vectorizado\n")
cat("Tiempo de ejecucion:", tiempo_lapply[[3]],"algoritmo vectorizado y lappy\n")
cat("Tiempo de ejecucion:", tiempo_split[[3]],"split\n")

## ------------------------------------------------------------------------
create_convex <- function(sp_db,sp_name,long,lat){
  ch <- chull(sp_db[,c(long,lat)])
  coords <- sp_db[c(ch, ch[1]), c(long,lat)]  
  coords <- cbind(sp_name,coords)
  return(coords)
}

## ------------------------------------------------------------------------
occs_bdList <- occs_bd %>% split(.$name)
poligonos <- occs_bdList %>% 
  purrr::map(~create_convex(sp_db = .x,
                           sp_name = .x$name[1],
                           long = "longitude",
                           lat = "latitude"))

## ------------------------------------------------------------------------
pol_subm <- poligonos[sample(length(poligonos),5)]
pol_subm %>% purrr::map(~ plot(.x$longitude,.x$latitude,type="l"))

## ------------------------------------------------------------------------
sessionInfo()

