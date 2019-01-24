# Script: 
# Instala los paquetes que usaremos en el curso.
# Autor: Luis Osorio-Olvera
# Fecha: 28/01/2019
# Curso: Curso Online de Nichos

# -----------------------------------------------------
# Nota muy importante a los usarios de Windows 
# INSTALEN RTools en
# https://cran.r-project.org/bin/windows/Rtools/

# -----------------------------------------------------
# Nota muy importante a los usarios de Mac instalen 
# las herramientas de desarrollador desde su terminal
# "sudo xcode-select --install"
# Escribir el pwd del usuario


# Paquetes que se intalar desde CRAN

# Instrucciones: copiar, pegar y correr el siguiente texto de instrucciones 
# en la consola de R o bien desde R-Studio como un script.

pkgs_curso <- c("devtools", # herramientas para desarrolladores
                "raster","leaflet", # Analisis espacial
                "rgdal","rasterVis",# Analisis espacial
                "maptools","sp", # Analisis espacial
                "spocc","dismo", # Distribucion de especies
                "ggplot2","rgl" , # Graficacion
                "dplyr","stringr","zoo", # manejo de datos
                "gdata", "XLConnect", # Lectura de datos
                "purrr", # Programacion funcional
                "tidyverse") # Paquetes para ciencia de datos (purrr,ggplot2, dplyr,readr...)


# Algoritmo que revisa que paquetes ya estan instalados en la computadora

pkgs_miss <- pkgs_curso[!(pkgs_curso %in% installed.packages())]

# Instalando paquetes 

if(length(pkgs_miss)>0L)
  install.packages(pkgs_miss, repos = "https://cloud.r-project.org/")

# Intalacion de paquetes que no estan en CRAN

devtools::install_github("luismurao/ntbox")
devtools::install_github("marlonecobos/kuenm")

install.packages("http://download.r-forge.r-project.org/src/contrib/demoniche_1.0.tar.gz",
                 repos=NULL,type="source")






