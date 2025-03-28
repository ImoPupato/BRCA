# Comparación de la expresión de ciertos genes de interés y análisis de sobrevida utilizando un set de datos de TCGA.
En este documento trabajaremos con el proyecto BRCA que contiene datos clínicos y de expresión de genes, entre otros, de muestras de cancer de mama humano.  
## Librerías utilizadas
### Descarga desde Cran
Primero descargamos los paquetes que hagan falta utilizando la función _install.packages()_ y luego las cargarmos:
```R
library(tidyverse)
library(survival)
library(survminer)
```
### Descarga desde Bioconductor 
Debemos instalar primero el paquete 'BiocManager' y luego los paquetes 'RTCGA', 'RTCGA.clinical' y 'RTCGA.mRNA':
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RTCGA")
BiocManager::install("RTCGA.clinical")
BiocManager::install("RTCGA.mRNA")
library(RTCGA.mRNA)
```
## Descarga de datos de sobrevida
Con la función _survivalTCGA()_ descargamos los datos clínicos que hacen referencia a la sobrevida de pacientes en el set de datos.
```R
clinica <- survivalTCGA(BRCA.clinical)
```
## Acceso, descarga y transformación de datos de expresión
```R
expr <- BRCA.mRNA %>%
  as_tibble() %>%
  select(bcr_patient_barcode, BRCA1, BRCA2, TP53, ESR1) %>% #genes que nos interesa comparar
  mutate(bcr_patient_barcode = str_sub(bcr_patient_barcode, 1, 12)) %>%
  inner_join(clinica, by = "bcr_patient_barcode") %>%
  pivot_longer(cols = c(BRCA1, BRCA2, TP53, ESR1), 
               names_to = "Gen", 
               values_to = "Expresión")
```
## Boxplot comparativo de la expresión de genes
```R
ggplot(expr, aes(x = Gen, y = Expresión, fill = Gen)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Expresión de genes en pacientes con BRCA",
       x = "Gen", 
       y = "Nivel de Expresión") +
  theme(legend.position = "none")
```
![](https://github.com/ImoPupato/BRCA/blob/main/Boxplot.png){width='100px'}


