#!/bin/bash
#SBATCH --partition={PARTITION_NAME}
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --account={ACCOUNT_NAME}
#SBATCH --job-name=STARa

# load modules
module load usc
module load star/2.7.0e
module load samtools/1.10
module load bedtools2/2.27.1

# your commands below
FILE=$1

### Filter reads ###
fastp -i fastq/${FILE}_R1.fastq.gz -I fastq/${FILE}_R2.fastq.gz \
                -o fastq_postQC/${FILE}_postQC_R1.fastq.gz -O fastq_postQC/${FILE}_postQC_R2.fastq.gz \
                --length_required 20 \
                --average_qual 20 \
                --detect_adapter_for_pe –correction \
                -h QCs/${FILE}_fastp.html \
                -j QCs/${FILE}_fastp.json

### Align reads ###
STAR --genomeDir ../STAR_index_mm39p \
        --runThreadN 8 \
        --readFilesCommand zcat \
        --readFilesIn fastq_postQC/${FILE}_postQC_R1.fastq.gz fastq_postQC/${FILE}_postQC_R2.fastq.gz \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --outFileNamePrefix STAR_alignments/${FILE}_ \
        --outSAMunmapped Within \
        --outSAMattributes Standard