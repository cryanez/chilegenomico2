#!/bin/bash
##########################################################################################
# Programmed by Cristian Yanez,  Engineer in bioinformatics, ChileGenomico Laboratory, 
# Faculty of Medicine, University of Chile.
#
# 2020-07-13
#
# Script to joins the ancestry calculations by windows into a single set ordered by descending European ancestry
#
##########################################################################################

start_value=$1								#beginning of admixed samples in the fam (row number)
end_value=$2								#number where the admixed samples end in the fam (row number)
j=$3
workDirectory=$4
reference=$5
CLG2_Ancestry=$6

listRefRandomSort="list_REF_RANDOM_SORT_"$reference".list"  		#list name in plink format of reference individuals + randomly admixed individuals ordered by predicted European ancestry                                                                                                                        
popinfoAllRandomSort="POP_REF_RANDOM_SORT_"$reference".csv"			#name of the popinfo of the set of reference individuals + admixed individuals ordered by predicted European ancestry
outConcatenated="REF"$(($start_value-1))							#variable name used to generate outputs where references samples are included
setPlinkRandom="$outConcatenated-RANDOM"							#Plink set name of reference individuals + admixed individuals randomly
setPlinkRandomSort=$setPlinkRandom"-SORT_"$reference				#name of the plink set of reference individuals + admixed individuals ordered by predicted European ancestry
popinfoREF="REF$(($start_value-1))_popinfo.csv"						#popinfo name ordered only with samples references
end_value_cat=$end_value"p"											#variable used to extract rows using sed command

##We create popinfo with ordered samples Reference + admixed samples ordered by European ancestry in descending order
echo -e "***We create popinfo with ordered samples Reference + admixed samples ordered by European ancestry in descending order:	$workDirectory/output/$popinfoAllRandomSort***"
cat $workDirectory/temporal/$popinfoREF <( paste $workDirectory/temporal/POP_RANDOM.csv <(sed -n $start_value,$end_value_cat $workDirectory/temporal/$setPlinkRandom.4.Q) | awk '{print $1"\t\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"0}' | sort -k 11 -r ) > $workDirectory/output/$popinfoAllRandomSort
##We create list of individuals in order of the previous popinfo
echo -e "***We create list of individuals in order of the previous popinfo:		$workDirectory/output/$listRefRandomSort***"
awk '{print $1"\t"$1}' $workDirectory/output/$popinfoAllRandomSort | tail -n +2 > $workDirectory/output/$listRefRandomSort
##We create Q file in the above order
echo -e "***We create Q file in the above order:				$workDirectory/output/$setPlinkRandomSort.4.Q***"
cat <(head -n $(($start_value-1)) $workDirectory/temporal/$setPlinkRandom.4.Q )  <(tail -n +$(($start_value+1)) $workDirectory/output/$popinfoAllRandomSort | awk '{print $10" "$11" "$12" "$13}' ) > $workDirectory/output/$setPlinkRandomSort.4.Q
##We create set plink in previous order
echo -e "***We create set plink in previous order:					$workDirectory/output/$setPlinkRandomSort***"
plink1.9 --bfile $workDirectory/temporal/$setPlinkRandom --indiv-sort f $workDirectory/output/$listRefRandomSort --make-bed --out $workDirectory/output/$setPlinkRandomSort
##Create Q file with Amerindian ancestry (only for calculations, not graphed)
echo -e "***Create Q file with Amerindian ancestry (only for calculations, not graphed):		$workDirectory/output/$setPlinkRandomSort-AMR.5.Q***"
awk '{print $0" "$3+$4}' $workDirectory/output/$setPlinkRandomSort.4.Q > $workDirectory/output/$setPlinkRandomSort-AMR.5.Q

##We graph complete set joined with ancestries calculated by blocks and ordered
echo -e "***We graph complete set joined with ancestries calculated by blocks and ordered:	$workDirectory/output/ADMIXTURE_K_4_$setPlinkRandomSort.pdf***"
cp $workDirectory/temporal/CV_$outConcatenated-50-$j.txt $workDirectory/output/CV_$setPlinkRandomSort.txt
Rscript code/New_admix.R $setPlinkRandomSort $popinfoAllRandomSort 2 $setPlinkRandomSort-AMR.5.Q $start_value $workDirectory $CLG2_Ancestry

##We create final ancestry archive file with Amerindian Ancestry included
echo -e "***We create final ancestry archive file with Amerindian Ancestry included:		$workDirectory/output/$reference-Ancestry.tsv***"
echo -e "ID_Muestra\tAFR\tEUR\tAYMARA\tCHILE_SUR\tAMERINDIO" > $workDirectory/output/$reference-Ancestry.tsv
tail -n +$(($start_value+1)) $workDirectory/output/$popinfoAllRandomSort | awk -F "\t" '{print $1"\t"$11"\t"$12"\t"$13"\t"$14"\t"$13+$14}' >> $workDirectory/output/$reference-Ancestry.tsv
