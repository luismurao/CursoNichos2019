# Script de introduccion a R
# Curso de modelos de nicho 2019
# Autor: Luis Osorio-Olvera
# Correo: luismurao@gmail.com

library(knitr)
library(rgl)
library(maptools)
data("wrld_simpl")

#' 
#' ### Resumen
#' 
#' En este tutorial vermos como utilizar R como un Sistema de Información Geográfica (SIG). 
#' R utiliza un conjunto de librerías para procesar objetos espaciales; las liberías que vamos a utilizar en este tutorial son:
#' 
#' - [sp](https://cran.r-project.org/web/packages/sp/index.html): Classes and methods for spatial data; the classes document where the spatial location information resides, for 2D or 3D data. Utility functions are provided, e.g. for plotting data as maps, spatial selection, as well as methods for retrieving coordinates, for subsetting, print, summary, etc.
#' - [rgdal](https://cran.r-project.org/web/packages/rgdal/index.html): Provides bindings to Frank Warmerdam's Geospatial Data Abstraction Library (GDAL) (>= 1.6.3) and access to projection/transformation operations from the PROJ.4 library. The GDAL and PROJ.4 libraries are external to the package, and, when installing the package from source, must be correctly installed first. Both GDAL raster and OGR vector map data can be imported into R, and GDAL raster data and OGR vector data exported. Use is made of classes defined in the sp package. Windows and Mac Intel OS X binaries (including GDAL, PROJ.4 and Expat) are provided on CRAN.
#' - [raster](https://cran.r-project.org/web/packages/raster/index.html): Reading, writing, manipulating, analyzing and modeling of gridded spatial data. The package implements basic and high-level functions. Processing of very large files is supported.
#' - [rgeos](https://cran.r-project.org/web/packages/rgeos/index.html): Interface to Geometry Engine - Open Source (GEOS) using the C API for topology operations on geometries. The GEOS library is external to the package, and, when installing the package from source, must be correctly installed first. Windows and Mac Intel OS X binaries are provided on CRAN.
#' - [leaflet](https://cran.r-project.org/web/packages/leaflet/index.html): Create and customize interactive maps using the 'Leaflet' JavaScript library and the 'htmlwidgets' package. These maps can be used directly from the R console, from 'RStudio', in Shiny apps and R Markdown documents.
#' - [rasterVis](https://cran.r-project.org/web/packages/rasterVis/index.html): Methods for enhanced visualization and interaction with raster data. It implements visualization methods for quantitative data and categorical data, both for univariate and multivariate rasters. It also provides methods to display spatiotemporal rasters, and vector fields. See the website for examples.
#' 
#' ![](figurasRSIG/mde.jpg)
#' 
#' En general trabajaremos sobre tareas específicas y vermos como ejecutarlas desde R. Finalmente se asume que los estudiantes ya tiene conocimientos básicos de SIG como:
#' 
#' ##### Vector
#' 
#' - Qué es un vector
#' - Cuantos tipos de vectores hay?
#'     - Puntos
#'     - Líneas
#'     - Polígonos
#' - Conocer la estructura interna de un vector
#'     - Número de vectores
#'     - Tabla de atributos
#'     
#' ##### Raster
#' 
#' - qué es un raster?
#' - Qué es un extent?
#' - Cuáles son sus dimensiones?
#' - Saber qué es una proyección?
#' - Saber si un raster está proyectado o no?
#' 
#' ### Asignar el directorio de trabajo
#' 
## ----eval=FALSE----------------------------------------------------------
## setwd("~/Dropbox/CursoNichos2019")

#' 
#' 
#' 
#' ## Trabajando con vectores desde R
#' 
#' Como se mencionó existen tre tipos de vectores espaciales; estos son puntos, líneas y polígonos. R utiliza un conjunto de librerías para procesar objetos espaciales. Las librerías que utilizamos para leer y hacer operaciones sobre vectores son: [`rgdal`](https://cran.r-project.org/web/packages/rgdal/index.html) y [`sp`](https://cran.r-project.org/web/packages/sp/index.html)
#' 
#' ### Cómo leo un vector desde R.
#' 
#' El comando para leer un shapefile es `readOGR`; este se encuentra en la librería [`rgdal`](https://cran.r-project.org/web/packages/rgdal/index.html). Los argumentos del comando son:
#' 
#' - `dsn`: Es el folder donde tenemos nuestro archivo vector (shapefile)
#' - `layer`: El nombre de la capa que queremos leer
#' 
## ------------------------------------------------------------------------
library(rgdal)
library(sp)
map_dir <- "destdv1gw"
map_vec <-  readOGR(dsn=map_dir ,layer = "destdv1gw")

#' 
#' Grafiquemos el mapa
#' 
## ----mpaita,cache=TRUE---------------------------------------------------
plot(map_vec)

#' 
#' ### Explorando los atributos del vector
#' 
#' Uno de los comandos de más útilies para explorar los elementos que conforman a un objerto, es el comando `str`.
#' 
## ----eval=FALSE----------------------------------------------------------
## str(map_vec)

#' 
#' o bien para darnos solo una idea de su estrutura podemos escribir el nombre del objeto
#' 
## ---- eval=FALSE---------------------------------------------------------
## map_vec

#' 
#' Los objetos espaciales en R son de la clase [**`S4`**](http://adv-r.had.co.nz/S4.html); para extraer los atributos de estos objetos (indexar)  se utiliza el símbolo de **@**. Veamos como extraer la tabla de atributos del mapa `map_vec`
#' 
## ------------------------------------------------------------------------
head(map_vec@data)

#' 
#' ### Número de vectores o poligonos
#' 
## ------------------------------------------------------------------------
length(map_vec@plotOrder)

#' 
#' ### Extraer un polígono del shapefile
#' 
## ------------------------------------------------------------------------
puebla <- map_vec[map_vec@data$ENTIDAD=="PUEBLA",]
puebla <-  map_vec[map_vec@data$ENTIDAD %in% "PUEBLA",]
plot(puebla)

#' 
#' ### Guardar un shapefile
#' 
#' Ahora guardamos como un shapefile separado, el polígono
#' 
## ----eval=FALSE----------------------------------------------------------
## writeOGR(obj = puebla,layer = "puebla",
##          dsn = "puebla",driver = "ESRI Shapefile")

#' 
#' 
#' ### Filtrar los puntos que caen dentro de un polígono
#' 
#' Generemos unos puntos aleatorios definidos en el extent del mapa; el comando `nombre_objeto@bbox` nos muestra el extent de nuestro shapefile.
#' 
## ------------------------------------------------------------------------
map_vec@bbox
# Coordenadas aleatorias
rand_x <- runif(n =300, map_vec@bbox[1,1],max = map_vec@bbox[1,2])
rand_y <- runif(n =300, map_vec@bbox[2,1],max = map_vec@bbox[2,2])
coords_rand <- data.frame(longitude=rand_x,latitude=rand_y)

#' 
#' #### Transformamar puntos a coordenadas espaciales
#' 
## ------------------------------------------------------------------------
# Cual es el CRS de map_vec
crs_mexico <- map_vec@proj4string@projargs
sp_df_rand <- SpatialPoints(coords = coords_rand,proj4string =CRS(crs_mexico))

#' 
#' 
#' 
#' #### Graficar mapa
#' 
## ----mapita,cache=T------------------------------------------------------
plot(map_vec)
plot(sp_df_rand,add=T)

#' 
#' #### Puntos en el poligono
#' 
## ----mapp, cache=TRUE----------------------------------------------------
data_poly_all <- sp::over( sp_df_rand , map_vec , fn = NULL)
en_poligono_index <- which(!is.na(data_poly_all[,1]))
p_en_poligono <- sp_df_rand[en_poligono_index ,]
plot(map_vec)
plot(p_en_poligono,add=T,col="blue")

#' 
#' #### Extraer los atributos en los puntos del poligono
#' 
## ------------------------------------------------------------------------
data_en_poli_extract <- data_poly_all[en_poligono_index,]
coords_en_poli <- sp_df_rand@coords[en_poligono_index,]
data_en_poli_extract_df <- cbind(coords_en_poli,data_en_poli_extract)
sp_DF_en_poli_extract <- SpatialPointsDataFrame(data_en_poli_extract_df[,c("longitude","latitude")],
                                                data_en_poli_extract_df)
#writeOGR(sp_DF_en_poli_extract,
#         dsn = map_dir,layer = "puntos_aleatorios_mexico",
#         driver = "ESRI Shapefile")

#' 
## ----mapll, cache=TRUE---------------------------------------------------
knitr::kable(head(sp_DF_en_poli_extract@data))
#plot(map_vec)
#cols <- unique(as.numeric(sp_DF_en_poli_extract@data$ENTIDAD))
#color_fun <- function(x) return()
#plot(p_en_poligono,add=T,col=sp_DF_en_poli_extract@data$ENTIDAD)

#' 
#' ### Reproyectar un polígono 
#' 
#' Con R podemos hacer reproyección de coordenadas, podemos pasar de coordenadas geográficas (logitude, latitud) a coordenadas utm (metrosal nivel del mar).
#' 
#' Veamo como proyectar nuesto mapa de Puebla a coordenadas utm
#' 
## ------------------------------------------------------------------------
puebla_reproj <- spTransform(puebla,CRS("+proj=utm +zone=10+datum=WGS84"))
plot(puebla_reproj,axes=TRUE)

#' 
#' Regresemos a coordenadas geográficas
#' 
## ------------------------------------------------------------------------
puebla_reproj_geo <- spTransform(puebla_reproj,CRS("+proj=longlat +ellps=WGS84"))
plot(puebla_reproj_geo,axes=TRUE)

#' 
#' 
#' 
#' ## Descargar datos espaciales desde R
#' 
#' Ahora veremos como descargar datos espaciales desde R; para ello utlizaremos la función [`getData`](https://www.rdocumentation.org/packages/raster/versions/2.5-8/topics/getData) del paquete `raster`. La función permite descargar datos de diferentes bases de datos como:
#' 
#' ### Límites territoriales
#' 
#' - [GADM](http://www.gadm.org/): Datos de límites territoriales de países
#' 
#' ### Datos climáticos
#' 
#' - [worldclim](http://www.worldclim.org/): Datos climáticos
#' 
#' ### Datos Digitales de Elevación
#' 
#' - [SRTM](http://srtm.csi.cgiar.org/): Datos Digitales de Elevación SRTM 90m.
#' - [mapzen](https://mapzen.com/):  Datos Digitales de Elevación aprox ~ 9m.
#' 
## ------------------------------------------------------------------------
library(raster)
wc10 <- raster::getData('worldclim', var='bio', res=10)
plot(wc10)

#' 
#' ## Trabajando con rasters desde R
#' 
#' La librería estrella para trabjar con capas raster es [`raster`]((https://cran.r-project.org/web/packages/raster/index.html)). Esta es un wrapper que se comunica con funciones de alto nivel escritas en el lenguaje `C` para procesar objetos raster. Es importante notar que los raster no son nada más que arreglos bidimensionales de tipo numérico o categórico. 
#' 
#' Al igual que los archivos de tipo vector tienen un extent, sin embargo tienen otros atributos como dimensión y resolución espacial. Entre mayor sea la resolución del raster, mayor será su dimensión (mayor número de filas y columnas).
#' 
#' ### Cómo leer un raster
#' 
#' El comando para leer 1 solo raster es `raster`
#' 
## ------------------------------------------------------------------------
library(raster)
bio2 <- raster("wc10/bio2.bil")
plot(bio2)

#' 
#' ### Cómo leer un conjunto (stack) de capas
#' 
#' A veces (casi siempre) será necesario leer un conjunto de capas que tienen la misma resolución y extent (***i.e.*** **[WorldClim layers](http://www.worldclim.org/)**). El comando `stack` de `raster` permite crear stack de capas; para ello es necesario indicar con en un vector los path o rutas a la carpeta que contiene nuestras capas.
#' 
#' Veamos como leer o almacenar en la memoria de R el conjunto de capas de **[WorldClim layers](http://www.worldclim.org/)**
#' 
#' Primero en guardamos en una variable los path a nustras capas
#' 
## ------------------------------------------------------------------------
paths_capas <- list.files("wc10/",
                          pattern = "*.bil$",full.names = TRUE)


#' 
#' Usamos el comando `stack`
#' 
## ------------------------------------------------------------------------
bios_wc <- raster::stack(paths_capas)

#' 
#' Graficamos las primer 5 capas
#' 
## ------------------------------------------------------------------------
plot(bios_wc[[1:5]])

#' 
#' notemos que utilice la misma sintaxis que utilizamos para indexa listas...
#' 
#' ### Hacer un crop
#' 
#' La función `crop` utiliza el extent de un shapfile para hacer recortar un raster o un stack. Veamos un recorte de las capas de [WorldClim](http://www.worldclim.org/) con el poligono de puebla.
#' 
#' 
## ------------------------------------------------------------------------
puebla_ras <- raster::crop(x = bios_wc,y = puebla)
plot(puebla_ras[[1]])
plot(puebla,add=T)

#' 
#' 
#' ### Hacer una mascara
#' 
#' Al igual que `crop`, `mask` hace un recorte pero conserva la forma del poligono.
#' 
## ------------------------------------------------------------------------
puebla_ras_forma <- raster::mask(puebla_ras, puebla)
plot(puebla_ras_forma,add=T)

#' 
#' ### Guardar un raster 
#' 
#' El comando para gurardar un raster es `wirteRaster` y está en el paquete `raster`
## ----eval=F--------------------------------------------------------------
## # Guardar el crop de puebla
## writeRaster(puebla_ras[[1]],"puebla_bio2.tif",overwrite=TRUE)

#' 
#' ### Guardar un stack completo
#' 
#' Mostraremos como guardar un stack completo en una carpeta
#' 
## ----eval=F--------------------------------------------------------------
## # Creamos la carpeta
## dir.create("Puebla_bios")
## lapply(names(puebla_ras_forma), function(x){
##   writeRaster(puebla_ras_forma[[x]], paste0("Puebla_bios/",
##                         x,".asc"),overwrite=TRUE)})

#' 
#' 
#' ### Convertir a puntos un stack
#' 
## ------------------------------------------------------------------------
stack_df <- data.frame(rasterToPoints(puebla_ras_forma))
head(stack_df)

#' 
#' 
#' ### Extraer valores de un raster o un stack
#' 
#' Podemos extraer los valores de un stack utilizando coordenadas.
#' 
## ------------------------------------------------------------------------
plot(bios_wc[[1]])
points(coords_rand,pch=".",col="red")
coord_rand_extract <- data.frame(extract(bios_wc,coords_rand))
coord_rand_extract <- na.omit(coord_rand_extract)
knitr::kable(head(coord_rand_extract))

#' 
#' 
#' Tabla de correlaciones...
#' 
## ------------------------------------------------------------------------
knitr::kable(round(cor(coord_rand_extract),2))

#' 
#' ### Reproyectar un raster
#' 
#' Al igual que con los vectores, podemos reproyectar a los rasters. Lo anterior se puede hacer mediante el comando `projectRaster` de `raster`
#' 
## ------------------------------------------------------------------------
puebla_bio1 <- puebla_ras[[1]]
puebla_bio1_reproj <- projectRaster(puebla_bio1,
                                    crs="+proj=utm +zone=10+datum=WGS84")
plot(puebla_bio1_reproj)
plot(puebla_reproj,add=T)

#' 
#' Regresemos a geográficas
#' 
## ------------------------------------------------------------------------
puebla_bio1_geo <- projectRaster(puebla_bio1, 
                                 crs="+proj=longlat +ellps=WGS84")
plot(puebla_bio1_geo)
plot(puebla,add=T)

#' 
#' 
#' ## La belleza de lo interactivo
#' 
#' En este ejemplo veremos como limpiar datos de precias de una especie con un mapa interactivo
#' 
## ------------------------------------------------------------------------
library(leaflet)
map <- leaflet() %>%
    addTiles()

map

#' 
#' 
#' 
#' Leemos los datos
#' 
## ----gbif1---------------------------------------------------------------
library(leaflet)
puma_occsDF <- read.csv("puma_concolor.csv")
head(puma_occsDF)

#' 
#' 1. Quitamos datos duplicados 
#' 
## ------------------------------------------------------------------------
library(ntbox)
puma_occsDF_L1 <- ntbox::clean_dup(puma_occsDF,
                                   longitude = "longitude",
                                   latitude = "latitude",
                                   threshold = 0) 

#' 
#' 
#' Visualizamos los datos 
#' 
## ------------------------------------------------------------------------

puma_occsDF_L1$rcID<- as.character(1:nrow(puma_occsDF_L1))
map %>% leaflet::addMarkers(lng =  puma_occsDF_L1$longitude,
                            lat = puma_occsDF_L1$latitude,
                            popup = puma_occsDF_L1$rcID)

#' 
#' 
#' 
#' Limpiamos datos 
#' 
## ------------------------------------------------------------------------
d_raros <- c(791,1367,966)
puma_occsDF_L2 <- puma_occsDF_L1[-d_raros ,]
map %>% leaflet::addMarkers(lng =  puma_occsDF_L2$longitude,
                            lat = puma_occsDF_L2$latitude,
                            popup = puma_occsDF_L2$rcID)

#' 
#' Escribimos nuestra tabla de datos 
#' 
## ------------------------------------------------------------------------
write.csv(puma_occsDF_L2, "puma_concolor_limpio.csv",row.names = F)

#' 
#' 
#' ### Descarga de datos de presencia `spocc`
#' 
#' 
## ----gbif----------------------------------------------------------------
library(spocc)

jaguar_occs <- occ("Panthera onca",
                    from = "gbif",
                    limit = 1500,
                    has_coords = T)
jaguar_occsDF <- occ2df(jaguar_occs)
map %>% leaflet::addCircleMarkers(lng =  jaguar_occsDF$longitude,
                                  lat = jaguar_occsDF$latitude)

#' 
#' Un poco más profesional 
#' 
## ------------------------------------------------------------------------
# install.packages("wicket")
library(wicket)
mex <- wrld_simpl[wrld_simpl$NAME %in% c("Mexico"),]
raster::extent(mex)
jaguar_mexExtent <- occ(query = 'Panthera onca',
                        geometry =c(-118.40417,14.55055,
                                    -86.70140,32.71846))
jaguar_mexExtentDF <- occ2df(jaguar_mexExtent)
map %>% leaflet::addCircleMarkers(lng =  jaguar_mexExtentDF$longitude,
                                  lat = jaguar_mexExtentDF$latitude)

#' 
#' Aún más pro
#' 
## ------------------------------------------------------------------------
jaguar_mex <- occ(query = 'Panthera onca', 
                  geometry = wicket::sp_convert(mex))
jaguar_mexDF <- occ2df(jaguar_mex)
map %>% leaflet::addCircleMarkers(lng =  jaguar_mexDF$longitude,
                                  lat = jaguar_mexDF$latitude)

#' 
#' 
#' ## Un ejemplo de programación funcional con `purrr`
#' 
#' 
#' - Leer la base de datos de reptiles y anfibios
#'     - Partir el `data.frame` por especie.
#'     - Crear polígnos mínimos convexos de los registros. 
#'     - Cortar las capas de `WorldClim` con el polígono.
#'     - Escribir los rasters (capas de modelación).
#'     - Buscar los datos de gbif usando polígonos.
#'     - Limpiar los duplicados espaciales.
#'     - Regresar el raster y los puntos de gbif.
#'     
## ----eval=TRUE-----------------------------------------------------------
library(magrittr)
d_sps <- read.csv("Arichivos_01Intro/occs_rep_amf.csv")
d_spsL <- d_sps %>% split(.$name)

#' 
#' Definimos la función para crear y escribir polígonos convexos en el ambiente global
#' 
## ------------------------------------------------------------------------
library(rgdal)
library(raster)
library(spocc)
# Objeto (funcion) en el ambiente global
create_convex <- function(sp_db,sp_name,long,lat,dir_sv){
  
  ch <- chull(sp_db[,c(long,lat)])
  coords <- sp_db[c(ch, ch[1]),
                  c(long,lat)]  
  sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(coords)),
                                           ID=1)))
  sp_poly_df <- SpatialPolygonsDataFrame(sp_poly,
                                         data=data.frame(ID=1))
  writeOGR(sp_poly_df, dir_sv, 
           layer=sp_name, 
           driver="ESRI Shapefile",
           overwrite_layer = TRUE,delete_dsn=TRUE)
  return(sp_poly_df)
}


#' 
#' Crando la función
#' 
## ------------------------------------------------------------------------
mi_funcion <- function(sp_db,sp_name,long,lat,dir_sv,env_layers,n_occs=500){
  dir_sv <- as.character(dir_sv)
  sp_name <- as.character(sp_name)
  if(!dir.exists(dir_sv)) dir.create(dir_sv)
  mi_ch <- create_convex(sp_db,sp_name,long,lat,dir_sv)
  sp_Layerscortadas <- crop(env_layers, mi_ch)
  sp_Layerscortadas <- mask(sp_Layerscortadas, mi_ch)
  
  lapply(names(sp_Layerscortadas), function(x){
  writeRaster(sp_Layerscortadas[[x]], paste0(dir_sv,"/",
                        x,".asc"),overwrite=TRUE)})
  
  sp_occs <- occ(sp_name,
                 from = "gbif",
                 limit = n_occs,
                 has_coords = T,
                 geometry = wicket::sp_convert(mi_ch))
  
  sp_occsDF <- occ2df(sp_occs)
  
  sp_occsLimpio <- ntbox::clean_dup(sp_occsDF,
                               longitude = "longitude",
                               latitude = "latitude",
                               threshold = raster::res(sp_Layerscortadas)[1])
  
  listaRes <- list(sp_occs=sp_occsLimpio,
                   layers= sp_Layerscortadas,
                   poligono=mi_ch)
  
  return(listaRes)
  
}

#' 
#' Veamos con una especie
#' 
## ------------------------------------------------------------------------
d_sps_DF <- d_spsL[[1]]
sp_ejemplo <- mi_funcion(sp_db = d_sps_DF,
                         sp_name = d_sps_DF$name[1],
                         long = "longitude",
                         lat = "latitude",
                         dir_sv = d_sps_DF$name[1],
                         env_layers = wc10,n_occs = 500)


#' 
#' Grafiquemos
#' 
## ------------------------------------------------------------------------
plot(sp_ejemplo$layers[[1]])
plot(sp_ejemplo$poligono,add=T)
points(sp_ejemplo$sp_occs[,c("longitude","latitude")],
       col="blue",pch=20)

#' 
#' Ahora hagámoslo con algunas (ustedes pueden con todas las sps)
#' 
## ------------------------------------------------------------------------
spsLista <- d_spsL[1:3]
spsResL <- spsLista %>% purrr::map(~mi_funcion(sp_db = .x,
                                    sp_name = .x$name[1],
                                    long = "longitude",
                                    lat = "latitude",
                                    dir_sv = .x$name[1],
                                    env_layers = wc10,
                                    n_occs = 500))


#' 
#' Exploremos la lista
#' 
## ------------------------------------------------------------------------
plot(spsResL$`Anaxyrus cognatus`$layers[[1]])
plot(spsResL$`Anaxyrus cognatus`$poligono,add=T)
points(spsResL$`Anaxyrus cognatus`$sp_occs[,c("longitude","latitude")],
       col="red",pch=20)

#' 
#' #### Revisen su directorio de trabajo
#' 
