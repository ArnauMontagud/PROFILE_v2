# Personalization of patient-specific treatments based on logical models (PROFILE_v2)

This is a repository of code and analyses related to the paper "Patient-specific Boolean models of signaling networks guide personalized treatments". 

The paper can accessed here: [https://www.biorxiv.org/content/10.1101/2021.07.28.454126v1](https://www.biorxiv.org/content/10.1101/2021.07.28.454126v1)

Present code is an extension to use the [PROFILE tool](https://github.com/sysbio-curie/PROFILE), paper available [here](https://www.frontiersin.org/articles/10.3389/fphys.2018.01965), to simulate patient-specific drug inhibitions to find patient-specific treatments.

## Brief tutorial on performing the gradient inhibition of nodes
### Requirements
- Python version 3.0 or greater
- Python's package numpy
- Perl
- R
- MaBoSS requires: flex, bison, gcc and g++
- A MaBoSS Boolean model

You need to have a MaBoSS Boolean model to be able to use this tool, such as the one provided (`LNCAP_mutRNA_EGF` BND and CFG files). Additionally, you may want to have a personalised model by using the [PROFILE tool](https://github.com/sysbio-curie/PROFILE). A tutorial is available [here](https://github.com/sysbio-curie/PROFILE/blob/master/Tutorial_PROFILE.pdf).

### Steps
A. In the `Gradient inhibition of nodes` folder, first prepare the files for the gradual inhibition of the nodes of interest:

1. `./drugs_loop_single.sh`

Script that builds `run_single.sh` that has one line for the gradual inhibition of each single node of interest.

2. `./drugs_loop_double.sh`

Script that builds `run_double.sh` that has one line for the  gradual inhibition of each combination of nodes of interest.

B. Launch the simulations
1. `run_single.sh`
2. `run_double.sh`

C. Gather the simulation results
1. `data_gathering_single.R`
2. `data_gathering_double.R`

D. Analyse and plot results
1. `data_analysis.R`

A processed dataframe is available to perform the analysis

## Scripts to reproduce figures of the paper
- Figure 2: Available in Supplementary File 1.
- Figure 3: Available in Supplementary File 2.
- Figure 4: Available in the `TCGA_plot.Rmd` script in `Analysis of TCGA patients' simulations` folder.
- Figure 5: Available in the `data_analysis.R` script in `Gradient inhibition of nodes` folder.
- Figure 6: Available in the `CL_plot.Rmd` script in `Analysis of drug sensitivities across cell lines` folder.
- Figures 7, 8 and 9: Available in the `drug_assays.R` script in `Analysis of experimental validation` folder.

Processed datasets are available to obtain figures 4 through 9.
