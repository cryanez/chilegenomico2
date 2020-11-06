##########################################################################################

 Documented by Cristian Yanez,  Engineer in bioinformatics, ChileGenomico Laboratory, 
 
 Faculty of Medicine, University of Chile.
 
 2020-07-13
 
##########################################################################################



# Descripción de GTSEQ_Ancestry_Estimate

GTSEQ_Ancestry_Estimate es una herramienta diseñada en el laboratorio Chilegenomico de la Facultad de Medicina de la Universidad de Chile. Estima la ancestría global de individuos chilenos mestizos utilizando ADMIXTURE. La estimación se realiza en muestras a las que se le genotipificaron 150 Marcadores informativos de la Ancestría (AIMS) chilena mediante la técnica Genotyping-in-Thousands by sequencing (GT-seq) de Nathan R Campbell (2015). En cada corrida se pueden genotipificar miles de muestras por lo que al estimar las ancestrías de un set en formato plink con esta cantidad de muestras chilenas, estas generan su propio componente ancestral (utilizando el algoritmo no supervisado de ADMIXTURE). Como solución al problema y así estimar de forma adecuada las ancestrías de muestras chilenas en el proyecto Chilegenomico2, la herramienta divide el set de miles de muestras mestizas en sets de 50 individuos chilenos seleccionados al azar o ventanas de 50 individuos (sin repetir muestras en ventanas). A cada ventana se le agregan 224 muestras de poblaciones que se utilizan como referencia para la estimación de las ancestrías en muestras chilenas.  Existen 4 componentes principales de ancestrías en Chile, estos son el componente africano, europeo, Aymara y del sur de Chile o Mapuche. Se utilizan como muestras referencias 30 africanos, 30 europeos, 110 Aymaras y 54 Mapuches (sumando 224 muestras en total). Por lo tanto se generan N sets compuestos por 224 muestras referencias + 50 muestras mestizas y se estima la ancestría global con ADMIXTURE para cada uno de forma independiente. 

El resultado principal de esta herramienta es la ancestría estimada de cada muestra chilena, la cual se reporta en una tabla ordenada por ancestría europea decreciente.

# Requerimientos

**Requirements:**

Bash

Python2.7

ADMIXTURE Version 1.3.0 

R

# Estructura de directorios

When downloading the code, you must create the missing "ancestry" and "input" directory, so that the directory structure remains as shown below:
```sh
working_directory
├── ancestry
│   └── CLG2
│       ├── input
│       ├── output
│       └── temporal
├── code
│   ├── 1_ancestry_windows_iterator.sh
│   ├── 2_ancestry_windows_calculation.sh
│   ├── 3_ancestry_windows_merge.sh
│   ├── get_Ancestry.sh
│   └── New_admix.R
├── input
│   ├── CLG2_Ancestry.tsv
│   ├── List_AIMS.csv
│   ├── list_REF224.list
│   ├── REF224_SORT.bed
│   ├── REF224_SORT.bim
│   ├── REF224_SORT.fam
│   └── REF224_popinfo.csv
└── README.txt
```

# Descripción de archivos de entrada y de salida
## Descripción de archivos de la carpeta "input"
Los datos de entrada para la estimación de ancestrías por ventanas de la carpeta "input" se deben solicitar al Dr. Ricardo Verdugo, Director del proyecto Chilegenomico, email: raverdugo@uchile.cl o a Cristian Yáñez, Ingeniero Bioinformático, email: cristianyanez@med.uchile.cl. 

La carpeta input contiene los siguientes archivos:

**CLG2_Ancestry.tsv:** Archivo tabulado con las ancestrías de las muestras chilenas a las que se les estimó la ancestría en el proyecto chilegenomico 1. Contiene 1909 muestras. La primera fila es el header y contiene 6 columnas. 
Las columnas son las siguientes: "ID_Interno", contiene el identificador de las muestras; "AFR", contiene la fracción de la ancestría africana para cada individuo; "EUR", contiene la fracción de ancestría europea; "AYMARA", fracción de ancestría Aymara, "CHILE_SUR", fracción del componente ancestral del sur de Chile o Mapuche; "AMERINDIO", fracción del componente amerindio compuesto por la suma de AYMARA + CHILE_SUR.

**List_AIMS.csv:** Lista de una columna sin header del total de los identificadores (rs ids) de los marcadores informativos de ancestrías (AIMS) utilizados en el análisis por la herramienta. 

**list_REF224.list:** Lista de los identificadores de 224 muestras referencias en formato plink utilizados para el cálculo de ancestrías por la herramienta. Contiene dos columnas separadas con un tabulador. Las muestras se ordenan por ancestría: 30 africanos, 30 europeos, 110 Aymaras y 54 Mapuches. 


**REF224_SORT.{bed,bim.fam}:** Set en formato plink de las 224 muestras referencias utilizadas para la estimación de ancestrías. Mismas muestras de list_REF224.list y en el mismo orden. Contienen 628.007 variantes. Para cada set al que se le estima la ancestría por ventanas con esta herramienta, la herramienta extrae las variantes en común a los AIMS del set de muestras mestizas de forma automática.

## Descripción de archivos de la carpeta "ancestry"
Explicación del contenido de la carpeta "ancestry". Esta carpeta debe ser creada siguiendo las siguientes instrucciones:

En la carpeta ancestry debe existir una carpeta por proyecto o set de datos al que se le estimará la ancestría en ventanas con la herramienta. En este caso tenemos la carpeta "ancestry/CLG2" que hace mención al proyecto Chilegenomico2 (CLG2). Dentro de la carpeta del proyecto deben existir 3 carpetas:

### Descripción de archivos de la carpeta "ancestry/project/input"
**CLG2/input:** Contiene el set en formato plink de las muestras 224 muestras referencias (las mismas de inputs/list_REF224.list) más las muestras mestizas chilenas a las que se les estimará las ancestrías con la herramienta (set de miles de muestras mestizas más las referencias). El set contiene los AIMS de las muestras mestizas que fueron genotipados por la técnica GT-Seq (150 AIMS) y que están presentes en más de un 80% de las muestras. Para las muestras referencias se obtienen las que están en común con estos AIMS a través del set inputs/REF224_SORT.{bed,bim.fam}. Este set debe estar ordenado por muestras, dejando las muestras referencias en el mismo orden que inputs/list_REF224.list al inicio y luego las muestras mestizas a las que se les estimará la ancestría.

Esta carpeta CLG2/input también debe contener un archivo tabulado con la información de la población para las muestras (popinfo) denominado: **popinfo_REF224_"NOMBRE PROYECTO".csv**. En el caso del proyecto CLG2 el archivo se denomina:  popinfo_REF224_CLG2.csv. Contiene una fila por muestra en el mismo orden que el set en formato plink al que se le estimará la ancestría por ventanas. Contiene las siguientes 16 columnas. IndID, Sex, Source, Region, Population, Ancestry, AncestryOriginal, IndLab, Recruitment, Reps, AFR, EUR, AYMARA, CHILE_SUR, AMERINDIO. Las últimas 6 columnas son la fracción de ancestría para las muestras. Para las muestras mestizas estas columnas se deben mantener vacías. 

### Descripción de archivos de salida de la carpeta "ancestry/project/output"
**CLG2/output:** Carpeta de salida de los resultados de la estimación de ancestría con la herramienta. A los archivos se le agrega el nombre del proyecto, en este caso CLG2. Los resultados que se obtienen son los siguientes:

**CLG2-Ancestry.tsv:** Archivo tabulado con las ancestrías de las muestras chilenas a las que se les estimó la ancestría por ventanas con la herramienta, utilizando los AIMS genotipados por GT-Seq. En el caso de CLG2 son 1720 muestras las que fueron genotipificadas por GT-Seq y que se les estimó la ancestría. En este archivo la primera fila es el header y contiene 6 columnas. Las columnas son las siguientes: "ID_Interno", contiene el identificador de las muestras; "AFR", contiene la fracción de la ancestría africana para cada individuo; "EUR", contiene la fracción de ancestría europea; "AYMARA", fracción de ancestría Aymara, "CHILE_SUR", fracción del componente ancestral del sur de Chile o Mapuche; "AMERINDIO", fracción del componente amerindio compuesto por la suma de AYMARA + CHILE_SUR. Los datos los entrega ordenados por ancestría europea decreciente. 

**ADMIXTURE_K_4_REF224-RANDOM-SORT_CLG2.pdf:** Gráfico de barras con las estimaciones de ancestrías usando los AIMS genotipados por GT-Seq y estimada con la herramienta. A la izquierda del gráfico están las 224 muestras referencias ordenadas y a la derecha las muestras mestizas a las que se les estimó las ancestrías y ordenadas por ancestría europea de forma decreciente.

**list_REF_RANDOM_SORT_CLG2.list:** Lista de los identificadores de las muestras en formato plink de las muestras en el mismo orden que CLG2-Ancestry.tsv

**POP_REF_RANDOM_SORT_CLG2.csv:** Archivo tabulado con la información de la población para las muestras (popinfo). Contiene una fila por muestra en el mismo orden del archivo anterior. Contiene las siguientes 16 columnas. IndID, Sex, Source, Region, Population, Ancestry, AncestryOriginal, IndLab, Recruitment, Reps, AFR, EUR, AYMARA, CHILE_SUR, AMERINDIO. Las últimas 6 columnas son la fracción de ancestría para las muestras. Para las muestras mestizas estas columnas reportan la ancestría estimada con la herramienta. 

**REF224-RANDOM-SORT_CLG2.4.Q:** Archivo en formato Admixture que reporta las fracciones de ancestría AFR, EUR, AYMARA, CHILE_SUR (en ese orden) estimadas por ventanas con la herramienta. Son cuatro columnas separadas por un espacio y mantienen el orden de filas (muestras) del archivo anterior.  

**REF224-RANDOM-SORT_CLG2-AMR.5.Q:** Archivo en formato Admixture que reporta las fracciones de ancestría AFR, EUR, AYMARA, CHILE_SUR, AMERNDIO (En ese orden) estimadas por ventanas con la herramienta. Son cinco columnas separadas por un espacio y mantienen el orden de filas (muestras) del archivo anterior.  La columna AMERINDIO es la suma de AYMARA + CHILE_SUR.


**REF224-RANDOM-SORT_CLG2.{bed, bim,fam}:** Set en formato plink ordenad, con los AIMS y las mismas muestras que el archivo anterior, conservando el mismo orden.

**Summary_ancestry_REF224-RANDOM-SORT_CLG2.csv:** Métricas resumen de la estimación de ancestría por ventanas. Contiene las estimaciones de ancestría promedio, desviación y error estándar. Valores referencias otenidos del proyecto CLG1 y la significancia de la diferencia entre ambos resultados.

### Descripción de archivos de la carpeta "ancestry/project/temporal"
**CLG2/temporal:** Carpeta donde se almacenan datos temporales


# Ejecución de herramienta GTSEQ_Ancestry_Estimate para estimar ancestría en miles de individuos por ventanas

Asegurese de estar en el directorio de trabajo.

### Ejecute el programa como se indica a continuación:
```sh
bash code/get_Ancestry.sh Name_Set_plink input/CLG2_Ancestry.tsv Project_folder_name Path_set_plink
Example:
bash code/get_Ancestry.sh CLG2 input/CLG2_Ancestry.tsv CLG2 PATH/CLG2_samples_GT-Seq
```
El script en bash code/get_Ancestry.sh recibe los siguientes 4 argumentos de entrada:

**Name_Set_plink:**
**input/CLG2_Ancestry.tsv:**
**Project_folder_name:**
**Path_set_plink:**
