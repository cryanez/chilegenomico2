
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
│   └── REF224_popinfo.csv
└── README.txt

In the inputs folder, currently the set used as reference set (REF224_SORT in plink format) for la estimación de ancestría en ventanas en muestras genotipadas por GT-SEQ is missing. If you want to use the script with the samples references used in the chilegenomco2 project, contact Cristian Yáñez, email: cristianyanez@med.uchile.cl
