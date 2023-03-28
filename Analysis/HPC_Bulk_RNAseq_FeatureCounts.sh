#!/bin/bash
#SBATCH --partition={PARTITION_NAME}
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=12:00:00
#SBATCH --account={ACCOUNT_NAME}
#SBATCH --job-name=FeatureC

# load modules
module load usc
module load star/2.7.0e
module load samtools/1.10
module load bedtools2/2.27.1

# your commands below
featureCounts -T 4 -p -t exon -g gene_name \
  -a gencode.vM29.primary_assembly.annotation.gtf \
  -o counts/${OUTPUT_NAME}_Read_Count_Table.txt \
  STAR_alignments/*_Aligned.sortedByCoord.out.bam