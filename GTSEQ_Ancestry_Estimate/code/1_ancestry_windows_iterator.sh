#!/bin/bash
##########################################################################################
# Programmed by Cristian Yanez,  Engineer in bioinformatics, ChileGenomico Laboratory, 
# Faculty of Medicine, University of Chile.
#
# 2020-07-13
#
# Script to automate the calculation of ancestry by windows for samples sequenced with the GT-Seq protocol. Script that allows organizing input data.
# Run script 2_ancestry_windows_calculation.sh to iterate the ancestry calculation in windows of 50 samples + references.
# Then run 3_ancestry_windows_merge.sh script which joins the ancestry calculations into a single set ordered by descending European ancestry
#
##########################################################################################

#Execute as: 
#bash 1_ancestry_windows_iterator.sh -p ../input/REF224_829GTS4_SORT -s 225 -e 1053 -i ../input/popinfo_REF224_829GTS4.csv -t 8 -d 1

setPlink=""									#plink set name (ordered with reference to start)
start_value=225								#beginning of admixed samples in the fam (row number)
end_value=2000								#number where the admixed samples end in the fam (row number)
popinfoAll=""								#name of the ordered popinfo of the set plink
N_CPU_ADMIX=8								#Number of processors (CPU) to use in admixture
workDirectory=""
reference=""
CLG2_Ancestry=""

listRefRandom="list_REF_RANDOM.list"		#list name in plink format of reference individuals + randomly admixed individuals
popinfoAllRandom="POP_REF_RANDOM.csv"		#popinfo name of the plink set of reference individuals + randomly mixed individuals
outConcatenated="REF"$(($start_value-1))	#variable name used to generate outputs, where references samples are included 
setPlinkRandom="$outConcatenated-RANDOM"	#Plink set name of reference individuals + admixed individuals randomly
listRef="list_REF$(($start_value-1)).list"	#list name in plink format of reference individuals
popinfoREF="REF$(($start_value-1))_popinfo.csv"	#popinfo name ordered only with samples references
end_value_cat=$end_value"p"					#variable used to extract rows using the sed command

usage() {                                      # Function: Print a help message.
  echo "Usage: $0 [-p NAME_SET_PLINK] [-s START_ADMIXED] [-e END_ADMIXED] [-i NAME_POPINFO] [-t THREADS_ADMIXTURE] [-d DELETE_OPTION] [-w WORK_DIRECTORY] [-r REFERENCE_PROYECT] [-a CLG2_ANCESTRY]" 1>&2 
}
exit_abnormal() {                              # Function: Exit with error.
  usage
  exit 1
}
while getopts ":p:s:e:i:t:d:w:r:a:h" options; do
  case "${options}" in                         
    p)
      setPlink=${OPTARG}
      ;;        
    s)
      start_value=${OPTARG}
      re_isanum='^[0-9]+$'                     		# Regex: match whole numbers only
      if ! [[ $start_value =~ $re_isanum ]] ; then 	# if $TIMES not a whole number:
        echo "Error: START_ADMIXED must be a positive, whole number."
        exit_abnormally
        exit 1
      elif [ $start_value -eq "0" ]; then      		# If it's zero:
        echo "Error: START_ADMIXED must be greater than zero."
        exit_abnormal                          		# Exit abnormally.
      fi      
      ;;
    e)
      end_value=${OPTARG}
      if ! [[ $end_value =~ $re_isanum ]] ; then  	# if $TIMES not a whole number:
        echo "Error: END_ADMIXED must be a positive, whole number."
        exit_abnormally
        exit 1
      elif [ $end_value -eq "0" ]; then           	# If it's zero:
        echo "Error: END_ADMIXED must be greater than zero."
        exit_abnormal                          		# Exit abnormally.
      fi            
      ;;   
    i)
      popinfoAll=${OPTARG}
      ;;          
    t)
      N_CPU_ADMIX=${OPTARG}
      re_isanum='^[0-9]+$'                     		# Regex: match whole numbers only
      if ! [[ $N_CPU_ADMIX =~ $re_isanum ]] ; then  # if $TIMES not a whole number:
        echo "Error: THREADS_ADMIXTURE must be a positive, whole number."
        exit_abnormally
        exit 1
      elif [ $N_CPU_ADMIX -eq "0" ]; then         	# If it's zero:
        echo "Error: THREADS_ADMIXTURE must be greater than zero."
        exit_abnormal                          		# Exit abnormally.
      fi            
      ;;   
    d)
      delete=${OPTARG}
      re_isanum='^[1-2]+$'                     		# Regex: match whole numbers only
      if ! [[ $delete =~ $re_isanum ]] ; then   	# if $TIMES not a whole number:
        echo "Error: DELETE_OPTION must be 1 or 2, whole number."
        exit_abnormally
        exit 1
      fi            
      ;;  
    w)
      workDirectory=${OPTARG}
      ;;     
    r)
      reference=${OPTARG}
      ;;    
    a)
      CLG2_Ancestry=${OPTARG}
      ;;                                
    h)
      echo "Usage: $0 [-p NAME_SET_PLINK] [-s START_ADMIXED] [-e END_ADMIXED] [-i NAME_POPINFO] [-t THREADS_ADMIXTURE] [-d DELETE_OPTION]"
      echo "-p		NAME_SET_PLINK: Name of set plink"
      echo "-s		START_ADMIXED"
      echo "-d		DELETE_OPTION: 1 to delete data temporal, 2 no delete data"
      exit_abnormal                            		# Exit abnormally.
      ;; 
    :)                                         		# If expected argument omitted:
      echo "Error: -${OPTARG} requires an argument."
      exit_abnormal                            		# Exit abnormally.
      ;;
    *)                                         		# If unknown (any other) option:
      exit_abnormal                            		# Exit abnormally.
      ;;
  esac
done

echo "####### ARGUMENTS ########"
echo "Set plink:	$setPlink"
echo "Start value:	$start_value"
echo "End value:	$end_value"
echo "Popinfo:	$popinfoAll"
echo "Threads admixture:	$N_CPU_ADMIX"
echo "Delete set temporal:	$delete"
echo "Work directory:		$workDirectory"
echo "##########################"

#Create job directories
if [ -d $workDirectory/output/ ];
then
echo "'output' directory already exists"
else
echo "Creating 'output' directory"
mkdir $workDirectory/output/
fi
if [ -d $workDirectory/temporal/ ];
then
echo "'temporal' directory already exists"
else
echo "Creating 'temporal' directory"
mkdir $workDirectory/temporal/
fi
echo "##########################"

##Create reference individuals list
echo "***Create reference individuals list:	$workDirectory/temporal/$listRef ***"
head -n $(($start_value-1)) $setPlink.fam | awk '{print $1"\t"$2}' > $workDirectory/temporal/$listRef
##We create popinfo for reference only
echo -e "***We create popinfo for reference only:		$workDirectory/temporal/$popinfoREF***"
head -n $start_value $popinfoAll > $workDirectory/temporal/$popinfoREF
##We create popinfo complete set REFERENCE + random admixed individuals
echo -e "***We create popinfo complete set REFERENCE + random admixed individuals:				$workDirectory/temporal/$popinfoAllRandom***"
end_pop_cat=$(($end_value+1))"p"
sed -n $(($start_value+1)),$end_pop_cat $popinfoAll | sort -R > $workDirectory/temporal/POP_RANDOM.csv
cat $workDirectory/temporal/$popinfoREF $workDirectory/temporal/POP_RANDOM.csv > $workDirectory/temporal/$popinfoAllRandom
##We create list of individuals of the new pop info with REFERENCE + random admixed individuals
echo -e "***We create list of individuals of the new pop info with REFERENCE + random admixed individuals:	$workDirectory/temporal/$listRefRandom***"
awk '{print $1"\t"$1}' $workDirectory/temporal/$popinfoAllRandom | tail -n +2 > $workDirectory/temporal/$listRefRandom
##We create a new set ordered according to the REFERENCE list + random admixed individuals
echo -e "***We create a new set ordered according to the REFERENCE list + random admixed individuals:		$workDirectory/temporal/$setPlinkRandom***"
plink1.9 --bfile $setPlink --indiv-sort f $workDirectory/temporal/$listRefRandom --make-bed --out $workDirectory/temporal/$setPlinkRandom
echo -e "\n***Calculation of ancestry in windows of REF + 50 admixed individuals"
j=0
for i in $(seq $start_value 50 $end_value)
do
  j=$(($j+1))
  echo "bash code/2_ancestry_windows_calculation.sh $setPlink $start_value $end_value $popinfoAll $N_CPU_ADMIX $j $i $workDirectory"  
  bash code/2_ancestry_windows_calculation.sh $setPlink $start_value $end_value $popinfoAll $N_CPU_ADMIX $j $i $workDirectory
done

bash code/3_ancestry_windows_merge.sh $start_value $end_value $j $workDirectory $reference $CLG2_Ancestry
echo "*** Delete option: $delete"
if [ "$delete" -eq 1 ]
then
echo "Deleting temporal folder data"
rm $workDirectory/temporal/*
else
echo "Conserving temporal folder data"
fi
echo "##########################"
exit 0

