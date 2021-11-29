#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=2:00:00
#SBATCH --partition=shas
#SBATCH --job-name=vcfconcat
#SBATCH --output=vcfconcat_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hazel.carr@ucdenver.edu

module purge
source /curc/sw/anaconda/default
conda activate birbenv

workdir=/scratch/summit/hhcarr@xsede.org/certhia21
cd ${workdir}

mkdir admix50kbpNTotal
mkdir admix100kbpNTotal
mkdir trees50kbpNTotal

bcftools concat --threads $SLURM_NTASKS -Ov -o admix50kbpNTotal/calls.vcf admix50kbpN/*.vcf

bcftools concat --threads $SLURM_NTASKS -Ov -o admix100kbpNTotal/calls.vcf admix100kbpN/*.vcf

bcftools concat --threads $SLURM_NTASKS -Ov -o trees50kbpNTotal/calls.vcf trees50kbpN/*.vcf