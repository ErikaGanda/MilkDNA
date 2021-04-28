# MilkDNA
Analysis and data supporting:  
["DNA extraction and host depletion methods significantly impact and potentially bias bacterial detection in a biological fluid"](https://www.biorxiv.org/content/10.1101/2020.08.21.262337v1) - Under review

The following scripts are as described below. Each input and output data set is indicated in the respective \*.Rmd files. Final output can be found in the respective HTML files.

For questions or additional information, please contact:  
Erika Ganda,  DVM  PhD
Assistant Professor, Food Animal Microbiomes  
Pennsylvania State University  
ganda at psu dot edu  

# The figure below provides a study overview:
Study overview. Each panel depicts the design for each aspect evaluated, with respective numbers of biological and technical replicates. The qPCR assays were performed in duplicate for each experimental replicate. Copy numbers were calculated for Bovine, Total 16S, Bacillus, Listeria, Mycobacterium, and Salmonella in (A); Bovine, Total 16S, Listeria, and Salmonella in (C), and only Bovine and Total 16S were calculated in (B) and (D) because no inoculation was performed.


![Ganda2021Figures1](https://user-images.githubusercontent.com/47567236/106322591-0cfae700-6244-11eb-828a-73937cb17182.jpg)

## A. DNAResults.Rmd			
Overview:  
This script will load qPCR results on DNA extraction with various methods, provide descriptive statistics, and generate visualizations for the paper indicated above. Linear model including all interactions is employed in `JDS_analysis.R`.

Experimental design is described in panel A of figure 1.

## B. PMAconcentrations.Rmd
Overview:  
This script will load qPCR results from experiments evaluating different PMA concentrations, employ a linear model to identify dose dependent effects, and generate visualizations for the paper indicated above.

Experimental design is described in panel B of figure 1.

## C. SelectiveLysisEnzymePMA.Rmd
Overview:
This script will load qPCR results from experiments evaluating different PMA concentrations of a mild protease coupled with PMA, employ a linear model to identify dose dependent effects, and generate visualizations for the paper indicated above.

Experimental design is described in panel C of figure 1.

## D. HostDepletionKitComparison.Rmd
Overview:
This script will load qPCR results from experiments evaluating different host DNA depletion methods, employ models to compare estimated marginal mean log copy numbers between methods, and generate visualizations for the paper indicated above.

Experimental design is described in panel C of figure 1.

## DNAconcentrations.Rmd
Overview:
This script will load data and calculate descriptive statistics of DNA concentrations.

## SequencingResults.Rmd
Overview:This script will load data and calculate descriptive statistics of a subset of 9 samples that were sequenced.

