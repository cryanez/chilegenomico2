
##########################################################################################
# Documented by Cristian Yanez,  Engineer in bioinformatics, ChileGenomico Laboratory, 
# Faculty of Medicine, University of Chile.
#
# 2020-07-13
##########################################################################################

Requirements:
Bash
Python2.7
ADMIXTURE Version 1.3.0 
R

When downloading the code, you must create the missing "ancestry" directory, so that the directory structure remains as shown below:

├── ancestry
│   └── CLG2
│       ├── input
│       ├── output
│       └── temporal
├── code
│   ├── 1_ancestry_windows_iterator.sh
│   ├── 2_ancestry_windows_calculation.sh
│   ├── 3_ancestry_windows_merge.sh
│   ├── get_Ancestry.sh
│   └── New_admix.R
├── input
│   ├── CLG2_Ancestry.tsv
│   ├── List_AIMS.csv
│   ├── list_REF224.list
│   ├── REF224_SORT.bed
│   ├── REF224_SORT.bim
│   ├── REF224_SORT.fam
│   └── REF224_popinfo.csv
└── README.txt

Los datos de entrada para la estimación de ancestrías por ventanas de la carpeta "input" se deben solicitar a Cristian Yáñez, email: cristianyanez@med.uchile.cl. 
La carpeta input contiene los siguientes archivos:

*CLG2_Ancestry.tsv: Archivo tabulado con las ancestrías de las muestras chilenas a las que se les estimó la ancestría en el proyecto chilegenomico 1. Contiene 1909 muestras. La primera fila es el header y contiene 6 columnas. 
Las columnas son las siguientes: "ID_Interno", contiene el identificador de las muestras; "AFR", contiene la frección de la ancestría africana para cada individuo;	"EUR", contiene la fracción de ancestría europea;	"AYMARA", fracción de ancestría Aymara,	"CHILE_SUR", fracción del componente ancestral del sur de Chile o Mapuche; "AMERINDIO", fracción del componente amerindio compuesto por la suma de AYMARA + CHILE_SUR.

*List_AIMS.csv: Lista de una columna sin header del total de los identificadores (rs ids) de los marcadores informativos de ancestrías (AIMS) utilizados en el análisis por la herramienta. 

*list_REF224.list: Lista de los identificadores de 224 muestras referencias en formato plink utilizados para el calculo de ancestrías por la herramienta. Contiene dos columnas separadas con un tabulador. Las muestras se ordenan por ancestría: 30 africanos, 30 europeos, 110 Aymaras y 54 Mapuches. 


*REF224_SORT.{bed,bim.fam}: Set en formato plink de las 224 muestras referencias utilizadas para la estimación de ancestrías. Mismas muestras de list_REF224.list y en el mismo orden. Contienen 628.007 variantes. Para cada set al que se le estima la ancestría por ventanas con esta herramienta, se extraen las variantes en común.



Explicación del contenido de la carpeta "ancestry". Esta carpeta debe ser creada siguendo las siguientes instrucciones:

En la carpeta ancestry debe existir una carpeta por proyecto o set de datos al que se le estimará la ancestría en ventanas con la herramienta. En este caso tenemos la carpeta "ancestry/CLG2" que hace mencion al proyecto Chilegenomico2 (CLG2). Dentro de la carpeta del proyecto deben existir 3 carpetas:

*CLG2/input: Contiene el set en formato plink de las muestras 224 muestras referencias (las mismas de inputs/list_REF224.list) más las muestras mestizas chilenas a las que se les estimará las ancestrías con la herramienta. El set cotiene los AIMS de las muestras mestizas que fueron genotipados por la técnica GT-Seq (150 AIMS) y que esán presentes en más de un 80% de las muestras. Para las muestras referencias se obtienen las que están en común con estos AIMS a travez del set inputs/REF224_SORT.{bed,bim.fam}. Este set debe estar ordenado por muestras, dejando las muestras referencias en el mismo orden que inputs/list_REF224.list al inicio y luego las muestras mestizas a las que se les estimará la ancestría.

Esta carpeta CLG2/input también debe conter un archivo tabulado con la información de la población para las muestras (popinfo) denominado: popinfo_REF224_"NOMBRE PROYECTO".csv. En el caso del proyecto CLG2 el archivo se denomina:  popinfo_REF224_CLG2.csv. Contiene una fila por muestra en el mismo orden que el set en formato plink al que se le estimará la ancestría por ventanas. Contiene las siguientes 16 columnas. IndID, Sex, Source, Region, Population, Ancestry, AncestryOriginal, IndLab, Recruitment, Reps, AFR, EUR, AYMARA, CHILE_SUR, AMERINDIO. Las últimas 6 columnas es la fracción de ancestría para las muestras. Para las muestras mestizas estas columnas se deben mantener vacías. 

*CLG2/output:


*CLG2/temporal









