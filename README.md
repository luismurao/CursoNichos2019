# CursoNichos2019

### Introducción a R

El presente repositorio tiene notas de programación en el lenguaje [R](https://www.r-project.org/). Se exploran desde las funciones más básicas del lenguaje (como leer tablas de datos, hacer aritmética, algunas operaciones algebraicas) hasta cómo utilizar funciones de paquetes especializados para el análisis espacial (*i.e.* [`raster`](https://cran.r-project.org/web/packages/raster/index.html), [`rgdal`](https://cran.r-project.org/web/packages/rgdal/index.html), [`dismo`](https://cran.r-project.org/web/packages/dismo/index.html)). 

También se hace énfansis en cómo programar funciones propias o funciones definidas por el usuario (**helper functions**) las cuales tienen el objetivo de facilitarnos tareas complejas. Es importante notar que operaciones básicas como la indexación suelen ser vitales para poder realizar operaciones mucho más complicadas y complejas dentro de los loops y estructuras de control.

### R como SIG

En este tutorial vermos como utilizar R como un Sistema de Información Geográfica (SIG). 
R utiliza un conjunto de librerías para procesar objetos espaciales; las liberías que vamos a utilizar en este tutorial son:

- [sp](https://cran.r-project.org/web/packages/sp/index.html): Classes and methods for spatial data; the classes document where the spatial location information resides, for 2D or 3D data. Utility functions are provided, e.g. for plotting data as maps, spatial selection, as well as methods for retrieving coordinates, for subsetting, print, summary, etc.
- [rgdal](https://cran.r-project.org/web/packages/rgdal/index.html): Provides bindings to Frank Warmerdam's Geospatial Data Abstraction Library (GDAL) (>= 1.6.3) and access to projection/transformation operations from the PROJ.4 library. The GDAL and PROJ.4 libraries are external to the package, and, when installing the package from source, must be correctly installed first. Both GDAL raster and OGR vector map data can be imported into R, and GDAL raster data and OGR vector data exported. Use is made of classes defined in the sp package. Windows and Mac Intel OS X binaries (including GDAL, PROJ.4 and Expat) are provided on CRAN.
- [raster](https://cran.r-project.org/web/packages/raster/index.html): Reading, writing, manipulating, analyzing and modeling of gridded spatial data. The package implements basic and high-level functions. Processing of very large files is supported.
- [rgeos](https://cran.r-project.org/web/packages/rgeos/index.html): Interface to Geometry Engine - Open Source (GEOS) using the C API for topology operations on geometries. The GEOS library is external to the package, and, when installing the package from source, must be correctly installed first. Windows and Mac Intel OS X binaries are provided on CRAN.
- [leaflet](https://cran.r-project.org/web/packages/leaflet/index.html): Create and customize interactive maps using the 'Leaflet' JavaScript library and the 'htmlwidgets' package. These maps can be used directly from the R console, from 'RStudio', in Shiny apps and R Markdown documents.
- [rasterVis](https://cran.r-project.org/web/packages/rasterVis/index.html): Methods for enhanced visualization and interaction with raster data. It implements visualization methods for quantitative data and categorical data, both for univariate and multivariate rasters. It also provides methods to display spatiotemporal rasters, and vector fields. See the website for examples.

### Literatura recomendada para aprender `R`

#### R como SIG

- Brunsdon,C. and Comber,L. (2019) An introduction to R for spatial analysis and mapping 2nd ed. SAGE Publishing.

- https://www.rspatial.org/

#### Programación 

##### Loops

- Braun,J. and Murdoch,D.J. (2007) A first course in statistical programming with R Cambridge University Press.

- Jones,O. (Owen D. et al. (2014) Introduction to scientific programming and simulation using R Chapman and Hall/CRC.

##### Funcional

- http://adv-r.had.co.nz/Functional-programming.html
- https://www.r-bloggers.com/functional-programming-in-r/
- https://www.slideshare.net/DASpringate/functional-programming-in-r

##### Vectorización

- http://alyssafrazee.com/2014/01/29/vectorization.html

#### Filogeografía, genética de poblaciones ...

- http://dyerlab.github.io/applied_population_genetics/mapping-populations.html
- https://arxiv.org/pdf/1604.01617.pdf
- https://cran.r-project.org/web/packages/phylin/vignettes/phylin_tutorial.pdf
- http://www.phytools.org/eqg/Exercise_3.2/
- https://www.youtube.com/channel/UCFXy3RKrbu595t8DyGjGoxg/videos
- http://dyerlab.github.io/popgraph/

#### Correr `python` desde `R`

- https://github.com/rstudio/reticulate

#### Guardar objetos de R

- https://www.r-bloggers.com/load-save-and-rda-files/
