#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --time=06:30:00
#SBATCH --partition=shas
#SBATCH --job-name=qcAln
#SBATCH --output=qcAln_%j.out
#SBATCH --mail-type=ALL
#SBATCH --array=4-51
#SBATCH --mail-user=hazel.carr@ucdenver.edu

module purge
source /curc/sw/anaconda/default
conda activate birbenv

workdir=/scratch/summit/hhcarr@xsede.org/certhia21
cd ${workdir}

mkdir qc
mkdir aln
mkdir bams
mkdir reports
mkdir finalBams

basename=$( head -n${SLURM_ARRAY_TASK_ID} basenames.txt | tail -n1 )
refgenome=/scratch/summit/hhcarr@xsede.org/certhia21/reference/06_certhia_reordered.fasta
gatk=/scratch/summit/hhcarr@xsede.org/certhia21/gatk-*/gatk

mkdir tmp

fastp -w $SLURM_NTASKS -q 20 -i input/${basename}_R1.fastq.gz -I input/${basename}_R2.fastq.gz \
    -o qc/${basename}_R1.fastq.gz -O qc/${basename}_R2.fastq.gz \
    -h reports/${basename}_fastpQC.html -j reports/${basename}_fastpQC.json

#align reads to reference
bwa mem -t $SLURM_NTASKS ${refgenome} qc/${basename}_R1.fastq.gz qc/${basename}_R2.fastq.gz > aln/${basename}.sam

rm qc/${basename}*

#convert sam to bam
samtools view -b -o bams/${basename}.bam aln/${basename}.sam

rm aln/${basename}.sam

$gatk CleanSam -I bams/${basename}.bam -O bams/${basename}_cleaned.bam

rm bams/${basename}.bam

$gatk --java-options "-Djava.io.tmpdir=/scratch/summit/hhcarr@xsede.org/${SLURM_ARRAY_TASK_ID} -Xmx32G" SortSam --TMP_DIR /scratch/summit/hhcarr@xsede.org/certhia21/tmp/${SLURM_ARRAY_TASK_ID} -SO coordinate -I bams/${basename}_cleaned.bam -O bams/${basename}_cleaned_sorted.bam

rm bams/${basename}_cleaned.bam

$gatk --java-options "-Djava.io.tmpdir=/scratch/summit/hhcarr@xsede.org/${SLURM_ARRAY_TASK_ID} -Xmx32G" AddOrReplaceReadGroups -I bams/${basename}_cleaned_sorted.bam -O bams/${basename}_cleaned_sorted_rg.bam --TMP_DIR /scratch/summit/hhcarr@xsede.org/certhia21/tmp/${SLURM_ARRAY_TASK_ID} --RGLB 1 --RGPL ILLUMINA --RGPU unit1 --RGSM ${basename}

rm bams/${basename}_cleaned_sorted.bam

$gatk --java-options "-Djava.io.tmpdir=/scratch/summit/hhcarr@xsede.org/${SLURM_ARRAY_TASK_ID} -Xmx60G" MarkDuplicates -ASO coordinate --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 -M reports/${basename}_markups_metrics_file.txt \
    -I bams/${basename}_cleaned_sorted_rg.bam -O finalBams/${basename}.bam

samtools index finalBams/${basename}.bam

