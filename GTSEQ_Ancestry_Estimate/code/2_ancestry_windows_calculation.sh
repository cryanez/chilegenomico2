#!/bin/bash
##########################################################################################
# Programmed by Cristian Yanez,  Engineer in bioinformatics, ChileGenomico Laboratory, 
# Faculty of Medicine, University of Chile.
#
# 2020-07-13
#
# Script to iterate the ancestry calculation in windows of 50 samples + references using a set plink
#
##########################################################################################

setPlink=$1									#plink set name (ordered with reference to start)
start_value=$2								#beginning of admixed samples in the fam (row number)
end_value=$3								#number where the admixed samples end in the fam (row number)
popinfoAll=$4								#name of the ordered popinfo of the set plink
N_CPU_ADMIX=$5								#Number of processors (CPU) to use in admixture
j=$6
i=$7
workDirectory=$8
listRefRandom="list_REF_RANDOM.list"		#list name in plink format of reference individuals + randomly admixed individuals
popinfoAllRandom="POP_REF_RANDOM.csv"		#popinfo name of the plink set of reference individuals + randomly mixed individuals
outConcatenated="REF"$(($start_value-1))	#variable name used to generate outputs, where references samples are included 
setPlinkRandom="$outConcatenated-RANDOM"	#Plink set name of reference individuals + admixed individuals randomly
listRef="list_REF$(($start_value-1)).list"	#list name in plink format of reference individuals
popinfoREF="REF$(($start_value-1))_popinfo.csv"	#popinfo name ordered only with samples references
end_value_cat=$end_value"p"					#variable used to extract rows using the sed command
StartAfr=1									#Start position of the african references samples in the fam (row number)
EndAfr=30"p"								#End position of African references samples in fam (row number)
StartEur=31									#Start position of the european references samples in the fam (row number)
EndEur=60"p"								#End position of european references samples in fam (row number)
StartAym=61									#Start position of the Aymara references samples in the fam (row number)
EndAym=170"p"								#End position of Aymara references samples in fam (row number)
StartMap=171								#Start position of the Mapuche references samples in the fam (row number)
EndMap=224"p"								#End position of Mapuche references samples in fam (row number)
start=$i
startPop=$(($start+1))
end=$(($i + 49))
concatenado_end=$end"p"
concatenado_endPop=$(($end+1))"p"

##We create a list with 50 admixed individuals
echo -e "\n***We create a list with 50 admixed individuals:				$workDirectory/temporal/list_$j.txt***"
sed -n $start,$concatenado_end $workDirectory/temporal/$setPlinkRandom.fam | awk '{print $1"\t"$2}' > $workDirectory/temporal/list_$j.txt
##We concatenate individuals references + 50 from the generated list
echo -e "***We concatenate individuals references + 50 from the generated list:	$workDirectory/temporal/$outConcatenated-50-$j.list***"
cat $workDirectory/temporal/$listRef $workDirectory/temporal/list_$j.txt > $workDirectory/temporal/$outConcatenated-50-$j.list
##We create popinfo with 50 individuals
echo -e "***We create popinfo with 50 individuals:					$workDirectory/temporal/popinfo_list_$j.csv***"
sed -n $startPop,$concatenado_endPop $workDirectory/temporal/$popinfoAllRandom > $workDirectory/temporal/popinfo_list_$j.csv
##We concatenate popinfo references + 50 from the generated list
echo -e "***We concatenate popinfo references + 50 from the generated list:		$workDirectory/temporal/popinfo_$outConcatenated-50-$j.csv***"
cat $workDirectory/temporal/$popinfoREF $workDirectory/temporal/popinfo_list_$j.csv > $workDirectory/temporal/popinfo_$outConcatenated-50-$j.csv
##We create reference subset + 50 from the generated list
echo -e "***We create reference subset + 50 from the generated list:		$workDirectory/temporal/$outConcatenated-50-$j***"
plink1.9 --bfile $workDirectory/temporal/$setPlinkRandom --keep $workDirectory/temporal/$outConcatenated-50-$j.list --indiv-sort f $workDirectory/temporal/$outConcatenated-50-$j.list --make-bed --out $workDirectory/temporal/$outConcatenated-50-$j
##We run admixture of the set with K 4
echo -e "***We run admixture of the set with K 4:				$workDirectory/temporal/$outConcatenated-50-$j.K4.log***"
admixture -j$N_CPU_ADMIX --cv $workDirectory/temporal/$outConcatenated-50-$j.bed 4 > $workDirectory/temporal/$outConcatenated-50-$j.K4.log
echo "mv *.P $workDirectory/temporal/"
mv *.P $workDirectory/temporal/
echo "mv *.Q $workDirectory/temporal/"
mv *.Q $workDirectory/temporal/
####### AFTER ADMIXTURE  #######
##We create CV file with the information obtained from admixture
echo -e "***We create CV file with the information obtained from admixture:	$workDirectory/temporal/CV_$outConcatenated-50-$j.txt***"
grep CV $workDirectory/temporal/$outConcatenated-50-$j.K4.log | awk '{ print 4, $4}' > $workDirectory/temporal/CV_$outConcatenated-50-$j.txt
##We graph ancestry by windows
echo -e "***We graph ancestry by windows:				$workDirectory/temporal/ADMIXTURE_K_4_$outConcatenated-50-$j.pdf***"
Rscript code/New_admix.R $outConcatenated-50-$j popinfo_$outConcatenated-50-$j.csv 1 $workDirectory

##Obtain columns of admxture Q files corresponding to each ancestry REFERENCE
echo -e "***Obtain columns of admxture Q files corresponding to each ancestry REFERENCE***"
sed -n $StartAfr,$EndAfr $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}'
ColsAfr=$(sed -n $StartAfr,$EndAfr $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}' | awk '{if( $1>$2 && $1>$3 && $1>$4){ print 1} else if($2>$1 && $2>$3 && $2>$4){print 2} else if($3>$1 && $3>$2 && $3>$4){print 3} else if($4>$1 && $4>$2 && $4>$3){print 4}  }')
echo "Column Afr in $workDirectory/temporal/$outConcatenated-50-$j.4.Q: $ColsAfr"
sed -n $StartEur,$EndEur $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}'
ColsEur=$(sed -n $StartEur,$EndEur $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}' | awk '{if( $1>$2 && $1>$3 && $1>$4){ print 1} else if($2>$1 && $2>$3 && $2>$4){print 2} else if($3>$1 && $3>$2 && $3>$4){print 3} else if($4>$1 && $4>$2 && $4>$3){print 4}  }')
echo "Column Eur in $workDirectory/temporal/$outConcatenated-50-$j.4.Q: $ColsEur"
sed -n $StartAym,$EndAym $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}'
ColsAym=$(sed -n $StartAym,$EndAym $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}' | awk '{if( $1>$2 && $1>$3 && $1>$4){ print 1} else if($2>$1 && $2>$3 && $2>$4){print 2} else if($3>$1 && $3>$2 && $3>$4){print 3} else if($4>$1 && $4>$2 && $4>$3){print 4}  }')
echo "Column Aymara in $workDirectory/temporal/$outConcatenated-50-$j.4.Q: $ColsAym"
sed -n $StartMap,$EndMap $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}'
ColsMap=$(sed -n $StartMap,$EndMap $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk '{sum1 += $1} {sum2 += $2} {sum3 += $3} {sum4 += $4} END {print sum1" "sum2" "sum3" "sum4}' | awk '{if( $1>$2 && $1>$3 && $1>$4){ print 1} else if($2>$1 && $2>$3 && $2>$4){print 2} else if($3>$1 && $3>$2 && $3>$4){print 3} else if($4>$1 && $4>$2 && $4>$3){print 4}  }')
echo "Column Mapuche in $workDirectory/temporal/$outConcatenated-50-$j.4.Q: $ColsMap"		
##Si es i=1 copiamos el archivo Q completo, si es 2 en adelante copiamos solo las muestras de 50 en 50

##We order Q files by Afr,Eur,Aym,Map columns
echo -e "***We order Q files by Afr,Eur,Aym,Map columns:		$workDirectory/temporal/$setPlinkRandom.4.Q***"
if [ "$j" -eq 1 ]
then
	awk -v c1=$ColsAfr -v c2=$ColsEur -v c3=$ColsAym -v c4=$ColsMap 'BEGIN{OFS=FS=" "} {print $c1,$c2,$c3,$c4}' $workDirectory/temporal/$outConcatenated-50-$j.4.Q > $workDirectory/temporal/$setPlinkRandom.4.Q
else
	concatenado_end_Q=$(($start_value+49))"p"
	sed -n $start_value,$concatenado_end_Q $workDirectory/temporal/$outConcatenated-50-$j.4.Q | awk -v c1=$ColsAfr -v c2=$ColsEur -v c3=$ColsAym -v c4=$ColsMap 'BEGIN{OFS=FS=" "} {print $c1,$c2,$c3,$c4}' >> $workDirectory/temporal/$setPlinkRandom.4.Q
fi
echo "*************************************"




