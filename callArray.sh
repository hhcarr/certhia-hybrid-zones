#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=24:00:00
#SBATCH --partition=shas
#SBATCH --job-name=callN
#SBATCH --output=callN_%j.out
#SBATCH --mail-type=ALL
#SBATCH --array=1-26
#SBATCH --mail-user=hazel.carr@ucdenver.edu

module purge
source /curc/sw/anaconda/default
conda activate birbenv

region=$( head -n${SLURM_ARRAY_TASK_ID} indices.txt | tail -n1 )
refgenome=/scratch/summit/hhcarr@xsede.org/certhia21/reference/06_certhia_reordered.fasta
workdir=/scratch/summit/hhcarr@xsede.org/certhia21
cd ${workdir}

mkdir variants/north

bcftools mpileup --threads $SLURM_NTASKS -r ${region} -d 0 -Ou -f ${refgenome} finalBams/north/*.bam | bcftools call --threads $SLURM_NTASKS -vmO u -o variants/north/calls_${region}.bcf

bcftools view -Ov variants/north/calls_${region}.bcf -o variants/north/calls_${region}.vcf
