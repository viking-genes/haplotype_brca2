#!/bin/bash
#$ -N brca2_hapl_len
#$ -l h_vmem=20G
#$ -l h_rt=02:30:00


cd ${haplotype_analysis_project_path}

# file containing the list of all carrier haplotypes, eg COHORT0001_1. The file has no header
hapl_list=brca2_gsa_carrier_haplotypes.txt

# selecting the current carrier haplotype
carr_hapl=$(head -n ${SGE_TASK_ID} ${hapl_list} | tail -n 1)

echo $carr_hapl

# calculate the length between each haplotype in the hapl_list file and all other carrier haplotypes. 
R --no-save -f scripts/a02_haplotype_length.r --args ${carr_hapl} 

