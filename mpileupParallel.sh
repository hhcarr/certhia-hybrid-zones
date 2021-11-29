#!/bin/bash

#SBATCH --time 02:00:00
#SBATCH --partition=shas
#SBATCH --ntasks=500
#SBATCH --job-name gpDemo
#SBATCH --output gnuparallel.out

module purge
module load gnu_parallel

source /curc/sw/anaconda/default
conda activate birbenv

workdir=/scratch/summit/hhcarr@xsede.org/certhia21
refgenome=$workdir/reference/06_certhia_reordered.fasta
cd ${workdir}

mkdir pileTmp

my_parallel="parallel --delay 5 -j $SLURM_NTASKS --colsep '\t' --tmpdir ${workdir}/pileTmp"
my_srun="srun --export=all --exclusive -n1 --cpus-per-task=1 --cpu-bind=cores"
$my_parallel "$my_srun bcftools mpileup -Ou -f ${refgenome} -r {} ${workdir}/finalBams/*.bam | bcftools call -m -v -Oz -o ${workdir}/{}.vcf.gz" {1} :::: indices.txt


#parallel --colsep '\t' --tmpdir ${workdir}/pileTmp samtools mpileup -b bams.fofn -r {1} :::: ${refgenome} > variants/calls.bcf
#parallel 'bcftools mpileup -Ou -f ${refgenome} -r {} ${workdir}/finalBams/*.bam | bcftools call -m -v -Oz -o {}.vcf.gz' ::: {1..${refgenome}.fai}

rm -rf pileTmp
echo done
#awk 'NF =1' /scratch/summit/hhcarr@xsede.org/certhia21/reference/06_certhia_reordered.fasta.fai > indices.txt
