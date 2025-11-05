#!/bin/bash

# running script
# run first the WTs
if [ ! -d "WT_sims" ] || [ ! -f "WT_sims/LNCAP_mutRNA_EGF.bnd" ] || [ ! -f "WT_sims/LNCAP_mutRNA_EGF.cfg" ]; then
    mkdir -p WT_sims
    echo Creating the WT_sims folder
    cp -r LNCAP_mutRNA_EGF* ./WT_sims/
fi
cd WT_sims
echo Running the WT model
../MBSS_FormatTable.pl LNCAP_mutRNA_EGF.bnd LNCAP_mutRNA_EGF.cfg 
cd ..

python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "AKT" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "AR" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "AR_ERG" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "Caspase8" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "EGFR" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "ERK" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "GLUT1" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "HIF1" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "HSPs" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "MEK1_2" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "MYC_MAX" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "PI3K" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "ROS" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "SPOP" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "TERT" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "cFLAR" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
python ./PROFILE_DrugSim.py LNCAP_mutRNA_EGF -d "p14ARF" -c "0, 0.2, 0.4, 0.6, 0.8, 1"
