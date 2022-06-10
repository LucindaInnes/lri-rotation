#!/bin/bash -l

#SBATCH --partition=panda
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=mapping
#SBATCH --time=24:00:00
#SBATCH --mem=32G

source ~/.bashrc

cd $1

mkdir output

pycoQC -f sequencing_summary* -o output/pycoQC_output.html


# porechop

conda activate /software/apps/porechop/0.2.3_seqan2.1.1

porechop -i fastq_pass > chopped.fastq

# minimap2 + samTools

conda activate myenv
conda install minimap2
conda install bedtools
module load samtools-1.13-gcc-8.2.0-nj2zx6n

minimap2 -ax map-ont $2 chopped.fastq | samtools view -bT $2 | samtools sort -o output/aln.bam
samtools index output/aln.bam output/aln.bai
genomeCoverageBed -ibam output/aln.bam > output/aln_coverage.txt


