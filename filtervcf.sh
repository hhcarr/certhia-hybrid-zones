#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=06:00:00
#SBATCH --partition=shas
#SBATCH --job-name=filtvcf
#SBATCH --output=filtvcf_%j.out
#SBATCH --array=1-26
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hazel.carr@ucdenver.edu

module purge
source /curc/sw/anaconda/default
conda activate birbenv

workdir=/scratch/summit/hhcarr@xsede.org/certhia21
cd ${workdir}
region=$( head -n${SLURM_ARRAY_TASK_ID} indices.txt | tail -n1 )

mkdir admix50kbpN
mkdir admix100kbpN
mkdir trees50kbpN

#gunzip -c variants/calls_filt.vcf.gz > variants/calls_filt.vcf

# run vcftools with SNP output spaced 50kbp 
vcftools --vcf variants/north/calls_${region}.vcf --remove-indv Certhia_familiaris_outgroup_KU92846 --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 50000 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out admix50kbpN/callsFilt${region}.vcf

# run vcftools with SNP output spaced 100kbp 
vcftools --vcf variants/north/calls_${region}.vcf --remove-indv Certhia_familiaris_outgroup_KU92846 --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 100000 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out admix100kbpN/callsFilt${region}.vcf

# run vcftools with SNP and invariant site output, 20% max missing data, no indels
vcftools --vcf variants/north/calls_${region}.vcf --max-missing 0.8 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out trees50kbpN/callsFilt${region}.vcf

# bgzip and tabix index files that will be subdivided into windows
# directory 1
bgzip trees50kbpN/callsFilt${region}.vcf.recode.vcf
#tabix
tabix -p vcf trees50kbpN/callsFilt${region}.vcf.recode.vcf.gz



#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )

# define main working directory
workdir=/lustre/scratch/jmanthey/01b_certhia_genomics2

# run vcftools with SNP output spaced 50kbp 
vcftools --vcf ${workdir}/03_vcf/${input_array}.g.vcf --remove-indv Certhia_familiaris_outgroup_KU92846 --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 50000 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/13_ds1_admix50kbp/${input_array}

# run vcftools with SNP output spaced 100kbp 
vcftools --vcf ${workdir}/03_vcf/${input_array}.g.vcf --remove-indv Certhia_familiaris_outgroup_KU92846 --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 100000 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/14_ds2_admix100kbp/${input_array}

# run vcftools with SNP and invariant site output, 20% max missing data, no indels
vcftools --vcf ${workdir}/03_vcf/${input_array}.g.vcf --max-missing 0.8 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out ${workdir}/15_ds3_trees50kbp/${input_array}

# run vcftools with SNP output with no missing data, no thinning
vcftools --vcf ${workdir}/03_vcf/${input_array}.g.vcf --remove-indv Certhia_familiaris_outgroup_KU92846 --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 6 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/17_ds5_lostruct/${input_array}

# bgzip and tabix index files that will be subdivided into windows
# directory 1
bgzip ${workdir}/15_ds3_trees50kbp/${input_array}.recode.vcf
# directory 2 while keeping original files
bgzip -c ${workdir}/17_ds5_lostruct/${input_array}.recode.vcf > ${workdir}/17_ds5_lostruct/${input_array}.recode.vcf.gz
#tabix
tabix -p vcf ${workdir}/15_ds3_trees50kbp/${input_array}.recode.vcf.gz
tabix -p vcf ${workdir}/17_ds5_lostruct/${input_array}.recode.vcf.gz
