### Seteo el directorio de trabajo
```R
setwd("G:/Mi unidad/INC") 
```
### Carga de librerías
```R
library(TCGAbiolinks) # para acceder a los datos de TCGA  
library(SummarizedExperiment) # para trabajar las matrices de diseños de experimentos  
library(edgeR) # para la transformación de cuentas crudas a CPM y el análisis de DEGs  
library(tidyverse) # para manipular los datos, incluye a dplyr, tidyr y ggplot2  
library(org.Hs.eg.db) # para convertir la notación de genes entre ESNSEMBL, SYMBOL, ENTREZ, etc  
library(ReactomePA) # para el ORA  
library("enrichplot") # para el gráfico de enriquecimiento  
```
### Generación de la búsqueda, descarga y creación de la tabla de expresión  
Aquí solamente utilizamos el paquete _TCGAbiolinks_. En caso de no tenerla previamente descargada, se puede hacer con el siguiente código:
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
```
Para más información sobre este paquete se puede acceder a la siguiente página de [Bioconductor](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html), o bien a través de R:
```
browseVignettes("TCGAbiolinks")
```
#### Búsqueda 
```R
query.exp <- GDCquery(
  project = "TCGA-BRCA", # datos pertenecientes al proyecto TCGA-BRCA (Breast Cancer)
  data.category = "Transcriptome Profiling", # transcriptos
  data.type = "Gene Expression Quantification", # cuantificación de expresión
  workflow.type = "STAR - Counts") # cuentas crudas
```
#### Descarga 
```R
GDCdownload(
  query = query.exp,
  files.per.chunk = 100)
```
#### Guardado del archivo
```R
BRCA.exp <- GDCprepare(
  query = query.exp,
  save = TRUE,
  save.filename = "BRCAExp.rda") # nombre del archivo, la extensión ".rda" corresponde a "datos de R"
```
### Acceso y exploración de la matriz de expresión (RNAseq)
En este paso utilizaremos la función 'assay' del paquete _SummarizedExperiment_:
```R
rnaseq <- assay(BRCA.exp) # para guardar la matriz de expresión en un objeto
```
Al explorar la matriz podemos ver cómo se asignan los códigos de cada muestra (nombre de columna) y cada gen (nombre de fila)
```R
rnaseq[0:2,0:2] # para poder ver los nombres de las filas y columnas, incluimos el 0  
                   TCGA-D8-A146-01A-31R-A115-07 TCGA-AQ-A0Y5-01A-11R-A14M-07
ENSG00000000003.15                         3414                          879
ENSG00000000005.6                           210                            9
```
El código de la muestra (nombre de la columna) está compuesto por letras y números que brindan información sobre el proyecto, el sitio de extracción de la muestra, etc. En el siguiente [link](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/) de la documentación del GDC se encuentra disponible la información para comprender cómo es la clasificación y qué información puedo obtener del código identificador. La descripción completa del código se encuentra en disponible en [TCGA Code Tables](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables). La correcta interpretación del código, nos permitirá seleccionar aquellas muestras de interés y generar una asignación correcta de las mismas para los análisis de contrastes entre grupos.  
El identificador del gen es el acordad por el Consorcio de Anotación de Genomas (Gencode) y está compuesto por una primera parte, previa al punto, que es el indentificador principal y luego corresponde a versiones o sub-identificaciones. Para nuestro análisis necesitamos únicamente identificador principal.  
A modo exploratorio podemos construir una tabla que contabilice la cantidad de tipos de muestra, según el criterio "sample type":
```R
table(substr(colnames(rnaseq[,-1]), 14, 15)) # para construir una tabla de conteo según el "sample type" indicado en los lugares 14 y 15 del código identificador
  01   06   11 
1110    7  113 
```
En nuestra matriz de expresión contamos con 1110 datos provenientes de tumores sólidos (01), 7 de tumores metastásicos (06) y 113 de tejido no tumoral (11).  

### Limpieza de datos
Aquí utilizaremos distintas funciones del paquete _tydiverse_ para acortar nombres y reasignar etiquetas
```R

```
```R

```
```R

```
