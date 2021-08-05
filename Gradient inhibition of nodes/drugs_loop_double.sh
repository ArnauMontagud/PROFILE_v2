#!/usr/bin/env bash
#title           :Build the running script for the double drug inhibition
#author      :Arnau Montagud, arnau.montagud@gmail.com
#date            :June 2021

model2="LNCAP_mutRNA_"
suffix=("EGF")
# suffix=("00" "AR" "EGF" "AR_EGF")
drugs=("AKT" "AR" "AR_ERG" "Caspase8" "cFLAR" "EGFR" "ERK" "GLUT1" "HIF1" "HSPs" "MEK1_2" "MYC_MAX" "p14ARF" "PI3K" "ROS" "SPOP" "TERT")
echo "# running script" >run_double.sh
for ((i = 0; i < ${#drugs[@]}; i++))
do
	for ((j = i; j < ${#drugs[@]}; j++)); do
		if [ ${drugs[i]} != ${drugs[j]} ]; then
				for k in ${suffix[@]}
				do
				model2_suffix=$model2$k
					echo "python ./MBSS_DrugSim5_mn4.py "$model2_suffix" -d \""${drugs[i]}, ${drugs[j]}"\" -c \"0, 0.2, 0.4, 0.6, 0.8, 1\"" >> run_double.sh
				done
			fi
	done
done

sbatch run_double.sh

