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
![](https://github.com/ImoPupato/BRCA/blob/main/Boxplot.png)
## Análisis de sobrevida respecto a la expresión de uno de los genes
### Selección del valor de corte
Hay diferentes criterios propuestos para elegir el valor de expresión a utilizar como corte para armar los grupos de expresión alta y baja, en este ejemplo utilizaremos la mediana.
```R
expr <- expr %>%
  group_by(Gen) %>% # calcula la mediana por cada grupo de genes
  mutate(Exp_group = ifelse(Expresión > median(Expresión, na.rm = TRUE), "Alto", "Bajo")) %>%
  ungroup() # desagrupamos para seguir con el mismo formato
```
### Creación del objeto y de la función necesarios para el gráfico de sobrevida
Utilizando la función _Surv()_ del paquete 'survival' indicamos el tiempo (_times_) hasta que se produce el evento (1) o pérdida de seguimiento (0), y el estado del paciente respecto de este evento (_patient.vital_status_).
```R
surv_obj <- Surv(time = expr$times, event = expr$patient.vital_status)
```
Con la función _survfit()_ del paquete 'survival' generamos la "curva de sobrevida" que nos permite indicar la probabilidad de sobrevida respecto del tiempo bajo las dos condiciones elegidas. Para este ejemplo utilizaresmo el gen "BRCA2".
```
fit <- survfit(surv_obj ~ Exp_group, data = expr, subset = (Gen == "BRCA2"))
El paquete 'survminer' nos permite graficar la curva de sobrevida.
```R
ggsurvplot(fit, data = expr, pval = TRUE,
           legend.labs = c("Baja Expresión", "Alta Expresión"),
           xlab = "Tiempo (dias)", ylab = "Probabilidad de Sobrevida")
```
![](https://github.com/ImoPupato/BRCA/blob/main/Sobrevida.png)
