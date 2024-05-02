### Seteo el directorio de trabajo
```R
setwd("G:/Mi unidad/INC") 
```
### Carga de librerías
```R
library("TCGAbiolinks") # para acceder a los datos de TCGA  
library(SummarizedExperiment) # para trabajar las matrices de diseños de experimentos  
library("RTCGA") # para acceder a los paquetes de TCGA para R
library("RTCGA.clinical") # para acceder a la información clínica
library(dplyr)
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
Una vez guardada la información del dataset, podemos acceder a ella con el cógido:
```R
load("G:/Mi unidad/INC/BRCAExp.rda")
```

### Acceso y exploración de la matriz de expresión (RNAseq) y la información clínica
- Matriz de expresión (RNASeq)  
En este paso utilizaremos la función 'assay' del paquete _SummarizedExperiment_:
```R
rnaseq <- assay(BRCA.exp) # para guardar la matriz de expresión en un objeto, si trabajamos desde BRCA.exp, sino utilizar el siguiente código
rnaseq <- `rownames<-`( # para asignar las etiquetas a las filas
  `colnames<-`( # para asignar las etiquetas a las columnas
    data@assays@data@listData[["unstranded"]], # contenido de la tabla
    data@colData@rownames), # etiquetas de las columnas
  data@rowRanges@ranges@NAMES) # etiquetas de las filas
write.table(rnaseq, "RNASeq(counts)_BRCA.txt") # para guardar la matriz en un archivo .txt
```
Al explorar la matriz podemos ver cómo se asignan los códigos de cada muestra (nombre de columna) y cada gen (nombre de fila)
```R
rnaseq[0:2,0:2] # para poder ver los nombres de las filas y columnas, incluimos el 0  
                   TCGA-D8-A146-01A-31R-A115-07 TCGA-AQ-A0Y5-01A-11R-A14M-07
ENSG00000000003.15                         3414                          879
ENSG00000000005.6                           210                            9
```
El código de la muestra (nombre de la columna) está compuesto por letras y números que brindan información sobre el proyecto, el sitio de extracción de la muestra, etc. En el siguiente [link](https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/) de la documentación del _National Cancer Institute_ se encuentra disponible la información para comprender cómo es la clasificación y qué información puedo obtener del código identificador. La descripción completa del código se encuentra en disponible en [TCGA Code Tables](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables).  
La correcta interpretación del código, nos permitirá seleccionar aquellas muestras de interés y generar una asignación correcta de las mismas para los análisis de contrastes entre grupos.  
El identificador del gen es el acordado por el Consorcio de Anotación de Genomas (Gencode) y está compuesto por una primera parte, previa al punto, que es el indentificador principal y luego corresponde a versiones o sub-identificaciones. Para nuestro análisis necesitamos únicamente identificador principal.  
A modo exploratorio podemos construir una tabla que contabilice la cantidad de tipos de muestra, según el criterio "sample type":
```R
table(substr(colnames(rnaseq), 14, 15)) # para construir una tabla de conteo según el "sample type" indicado en los lugares 14 y 15 del código identificador
  01   06   11 
1111    7  113 
```
En nuestra matriz de expresión contamos con 1111 datos provenientes de tumores sólidos (01), 7 de tumores metastásicos (06) y 113 de tejido no tumoral (11).  
Para este análisis vamos a utilizar únicamente las muestras que provengan de tumores sólidos.

- Información clínica
Aquí utilizaremos el paquete _RTCGA.clinical_ que nos permite descargar directamente la base:
```R
clinica<-as.data.frame(BRCA.clinical) # para guardar la información clínica en un objeto de tipo data frame
dim(clinica)
[1] 1098 3703
```
Al explorar el data frame podemos ver que está compuesto por 1098 filas y 3703 columnas. Variables de interés pueden ser el sexo asignado al nacer, estadio de menopausia, etnia y raza, etc.:
```R
table(clinica$patient.gender)
female   male 
  1086     12 
```
```R
table(clinica$patient.ethnicity)
    hispanic or latino not hispanic or latino
                    39                    885 
```
```R
table(clinica$patient.race)
american indian or alaska native                            asian        black or african american 
                               1                               61                              183 
                           white 
                             758 
```
```R
table(clinica$patient.menopause_status)

                                               indeterminate (neither pre or postmenopausal) 
                                                                                          34 
                                              peri (6-12 months since last menstrual period) 
                                                                                          39 
           post (prior bilateral ovariectomy or >12 mo since lmp with no prior hysterectomy) 
                                                                                         705 
pre (<6 months since lmp and no prior bilateral ovariectomy and not on estrogen replacement) 
                                                                                         230                  39                    885
```
### Limpieza de datos
Antes de proceder a la limpieza de los datos, debemos preguntarnos cuál es nuestro objetivo. En este caso queremos contrastar la expresión diferencial de RNASeq en muestras de tumores de mama sólidos de acuerdo al estadío pre y post menupáusico.  
Para cumplir nuestro objetivo necesitamos tener los genes correctamente asignados (sin variantes), eliminar los registros de muestras no tumorales y tumores metastásicos, asignar la situación menopáusica.
- Recorte de los nombres de genes, remoción de muestras de tejido no tumoral y tumores metastásicos, recorte de códigos en los nombres de las muestras y orden alfabético en la matriz de expresión.
Aquí utilizaremos distintas funciones del paquete _tydiverse_ para acortar nombres, reasignar etiquetas y seleccionar datos:
```R
rnaseq <- rnaseq %>%
  `rownames<-`(sub(   # la función sub sirve para reemplazar mediante el criterio "match and replace", 
    "\\..*",  # el argumento "\\..*" es el match, un punto seguido de cualquier cantidad de caracteres
    "", # el argumento "" es el replace 
    rownames(.))) %>% # objeto a ser reemplazado
  `[`(, !(substr(colnames(.), 14, 15) == "11")) %>% # debemos tomar como caracter al número 11 porque el código posee letras y numeros
  `colnames<-`(substr(colnames(.), 1, 12)) %>% # para acortar el nombre para que sea igual al utilizado en la información clínica
  `[`(, order(colnames(.))) # ordenar alfabéticamente
```
- Asignación de etiquetas para el contraste y orden alfabético en la tabla de información clínica.
Aquí utilizaremos las funciones 'mutate', 'case_when' y 'arrange' del paquete _dplyr_:
```R
clinica <- clinica %>%
  mutate(
    estadio.menop = case_when(
      patient.menopause_status == "post (prior bilateral ovariectomy or >12 mo since lmp with no prior hysterectomy)" ~ "Post",
      patient.menopause_status == "pre (<6 months since lmp and no prior bilateral ovariectomy and not on estrogen replacement)" ~ "Pre",
      TRUE ~ "Indeterminada"
    ),
    patient.bcr_patient_barcode = toupper(patient.bcr_patient_barcode) # para que coincida con las etiquetas de rnaseq
  ) %>%
  arrange(patient.bcr_patient_barcode) # para ordenar alfabeticamente
```
- Generación del _sub_ set de datos con información clínica
```R
rnaseq <- rnaseq[, # para seleccionar todas las filas
                 colnames(rnaseq) # para seleccionar las columnas cuyo nombre cumple con el siguiente requisito
                 %in% clinica$patient.bcr_patient_barcode # para que solo selecciones según la información clínica
                 & !duplicated(colnames(rnaseq)) # para que no incluya valores duplicados en caso que los haya
                 & clinica$estadio.menop != "Indeterminada"] # para que no incluya a quienes no tienen un estadio menopausico asignado
```
- Generación del vector de etiquetas
```R
etiquetas <- subset(clinica$estadio.menop, clinica$patient.bcr_patient_barcode %in% colnames(rnaseq))
```
### Análisis de la expresión diferencial
Para este análisis utilizaremos el paquete edgeR de Bioconductor, si bien hay una breve explicación de cada línea del código, para mayor comprensión se recomienda leer la [guía de usuario](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).
```R
y<-DGEList(rnaseq) # creación del objeto DGE para el análisis, utilizamos las cuentas crudas de RNASeq
y[["samples"]]$group<-etiquetas # asignamos las etiquetas de cada columna para hacer el contraste
keep <- filterByExpr(y) # utilizando el método propuesto por el paquete edgeR para decidir qué "probs" se quedan y cuales deben eliminarse, se genera un vector lógico
y <- y[keep,,keep.lib.sizes=FALSE] # nos quedamos con aquellas que cumplen los requisitos para ser incluidas en el análisis
y <- calcNormFactors(y) # se determinan los factores de normalización
design<-model.matrix(~etiquetas) # se construye el diseño del contraste
y<-estimateDisp(y,design) # estimación de la dispersión o análisis propiamente dicho
et<-exactTest(y) # creación del objeto con la salida del análisis estadístico 
is.de<-decideTests(et) # con la función decideTest aplicada al objeto et, se indica cuáles son los diferencialmente expresados
summary(is.de) # resumen del análisis
salida<-as.data.frame(topTags(et,n=Inf)) #guardar en un objeto la salida del análisis
salida$Expresion = ifelse(salida$FDR < 0.01 & abs(salida$logFC) >= 1, 
                       ifelse(salida> 1 ,'Up','Down'),
                       'Stable') # asignamos la expresión según nuestro criterio
table(salida$Expresion)
write.table("salida","Exact test BCRA vs estadio menop.txt")
```
### Analisis de sobrerrepresentación de vías
Para este análisis utilizaremos los paquetes org.Hs.eg.db; ReactomePA para el análisis de sobrerrepresentación de vías y enrichplot para el gráfico de enriquecimiento.
```R
genes.de<-c(subset(rownames(salida),salida$Expresion!="Stable"))

???vias.enriquecidas <- enrichPathway(gene=select(org.Hs.eg.db,genes.de,"ENTREZID",pvalueCutoff=0.05, readable=T))
vav3<-c(subset(x$row_names,x$Expression!="Stable"))
my.symbols <- vav3

# Seleccionamos lo que queremos "comparar", en este caso Reactome necesita el ENTREZ y tenemos en SYMBOL 
IDS<-select(hs, 
            keys = my.symbols,
            columns = c("ENTREZID", "SYMBOL"),
            keytype = "SYMBOL")

#Seleccionamos el ENTREZ que es lo que usa Reactome
names_ids<-IDS$ENTREZID

# Hacemos el Enrichment analysis
x <- enrichPathway(gene=names_ids,pvalueCutoff=0.05, readable=T)
head(as.data.frame(vias.enriquecidas))

#Graficamos el dotplot
p<-dotplot(x, showCategory=26, font.size = 10)
```
### Bibliografía
