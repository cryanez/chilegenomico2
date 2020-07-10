#!/bin/bash
##########################################################################################
# Programmed by Cristian Yanez,  Engineer in bioinformatics, ChileGenomico Laboratory, 
# Faculty of Medicine, University of Chile.
#
# 2020-07-13
##########################################################################################

shopt -s expand_aliases
source ~/.bashrc
conda activate anaconda3_python2

nameSet=$1
CLG2_Ancestry=$2
name_folder_set=$3
ruta_set_plink=$4

#Execute as:
#bash code/get_Ancestry.sh Nombre_Set_plink input/CLG2_Ancestry.tsv Nombre_carpeta_etapa4_set ruta_set_plink  
#bash code/get_Ancestry.sh GTS5_BiLS input/CLG2_Ancestry.tsv BiLs ../1_crear_sets_por_proyectos/SETS_FINALES/BiLs/GTS5_BiLS

echo $nameSet
echo $CLG2_Ancestry
echo $name_folder_set
echo $ruta_set_plink

#Extract AIMS
plink1.9 --bfile $ruta_set_plink --extract "input/List_AIMS.csv" --make-bed --out "ancestry/"$name_folder_set"/"$nameSet"_AIMS"
plink1.9 --bfile "ancestry/"$name_folder_set"/"$nameSet"_AIMS" --missing --out "ancestry/"$name_folder_set"/"$nameSet"_AIMS_missing"
awk '$5 > 0.2 {print $2}' "ancestry/"$name_folder_set"/"$nameSet"_AIMS_missing.lmiss" | tail -n +2 > "ancestry/"$name_folder_set"/list_SNPs_missing_"$nameSet".list"
plink1.9 --bfile "ancestry/"$name_folder_set"/"$nameSet"_AIMS" --exclude "ancestry/"$name_folder_set"/list_SNPs_missing_"$nameSet".list" --make-bed --out "ancestry/"$name_folder_set"/"$nameSet"_AIMS_filter1"
comm -12 <(cut -f 2 "ancestry/"$name_folder_set"/"$nameSet"_AIMS_filter1.bim" | sort | uniq) <(cut -f 2 "input/REF224_SORT.bim" | sort | uniq) > "ancestry/"$name_folder_set"/list_SNPs_comm.list"
plink1.9 --bfile "ancestry/"$name_folder_set"/"$nameSet"_AIMS_filter1" --extract "ancestry/"$name_folder_set"/list_SNPs_comm.list" --make-bed --out "ancestry/"$name_folder_set"/"$nameSet"_AIMS_filter2"
plink1.9 --bfile "input/REF224_SORT" --extract "ancestry/"$name_folder_set"/list_SNPs_comm.list" --indiv-sort f "input/list_REF224.list" --make-bed --out "ancestry/"$name_folder_set"/REF224_SORT_AIMS"
#We create samples popinfo
cat "input/REF224_popinfo.csv" > "ancestry/"$name_folder_set"/popinfo_REF224_"$nameSet".csv"
awk '{print $1}' "ancestry/"$name_folder_set"/"$nameSet"_AIMS_filter2.fam" | xargs -I{} echo -e {}"\t\tunknown\tunknown\tunknown\tadmixed\tadmixed\t"{}"\tunknown\t1\t\t\t\t\t" >> "ancestry/"$name_folder_set"/popinfo_REF224_"$nameSet".csv"
cut -f 1 "ancestry/"$name_folder_set"/popinfo_REF224_"$nameSet".csv" | tail -n +2 | xargs -I{} echo -e {}"\t"{} > "ancestry/"$name_folder_set"/list_REF224_"$nameSet".list"
plink1.9 --bfile "ancestry/"$name_folder_set"/"$nameSet"_AIMS_filter2" --bmerge "ancestry/"$name_folder_set"/REF224_SORT_AIMS" --make-bed --out "ancestry/"$name_folder_set"/REF224_"$nameSet
plink1.9 --bfile "ancestry/"$name_folder_set"/REF224_"$nameSet --indiv-sort f "ancestry/"$name_folder_set"/list_REF224_"$nameSet".list" --make-bed --out "ancestry/"$name_folder_set"/REF224_"$nameSet"_SORT"
mkdir "ancestry/"$name_folder_set"/input"
cp "ancestry/"$name_folder_set"/REF224_"$nameSet"_SORT.bed" "ancestry/"$name_folder_set"/REF224_"$nameSet"_SORT.bim" "ancestry/"$name_folder_set"/REF224_"$nameSet"_SORT.fam" "ancestry/"$name_folder_set"/input"
cp "ancestry/"$name_folder_set"/popinfo_REF224_"$nameSet".csv" "ancestry/"$name_folder_set"/input"
#Ancestry calculation by windows
endAdmixed=$(wc -l "ancestry/"$name_folder_set"/input/REF224_"$nameSet"_SORT.fam" | awk '{print $1}')
bash code/1_ancestry_windows_iterator.sh -p "ancestry/"$name_folder_set"/input/REF224_"$nameSet"_SORT" -s 225 -e $endAdmixed -i "ancestry/"$name_folder_set"/input/popinfo_REF224_"$nameSet".csv" -t 8 -d 1 -w "ancestry/"$name_folder_set -r $nameSet -a $CLG2_Ancestry

