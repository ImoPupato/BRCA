### Seteo el directorio de trabajo
setwd("G:/Mi unidad/INC") 

### Carga de librerías
library(TCGAbiolinks) # para acceder a los datos de TCGA
library(SummarizedExperiment) # para trabajar las matrices de diseños de experimentos
library(edgeR) # para la transformación de cuentas crudas a CPM y el análisis de DEGs
library(tidyverse) # para manipular los datos, incluye a dplyr, tidyr y ggplot2
library(org.Hs.eg.db) # para convertir la notación de genes entre ESNSEMBL, SYMBOL, ENTREZ, etc
library(ReactomePA) # para el ORA
library("enrichplot") # para el gráfico de enriquecimiento

### Generación de la búsqueda, descarga y creación de la tabla de expresión
#### Búsqueda 
```R
query.exp <- GDCquery(
  project = "TCGA-BRCA", # datos pertenecientes al proyecto TCGA-BRCA (Breast Cancer)
  data.category = "Transcriptome Profiling", # transcriptos
  data.type = "Gene Expression Quantification", # cuantificación de expresión
  workflow.type = "STAR - Counts") # cuentas crudas
```
GDCdownload( # para descargar la búsqueda
  query = query.exp,
  files.per.chunk = 100)

BRCA.exp <- GDCprepare( # para guardar la descarga
  query = query.exp,
  save = TRUE,
  save.filename = "BRCAExp.rda") 

rnaseq <- assay(BRCA.exp) # creo el objeto rnaseq con la información descargada
