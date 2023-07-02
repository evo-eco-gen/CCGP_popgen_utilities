#!/bin/bash
### Evaluate, filter and map paired end Illumina reads to a Perognathus longimembris longimembris CCGP reference, using thw Savio HPC at the UC.
# Job name:
#SBATCH --job-name=read_map
#
# Account:
#SBATCH --account=fc_nachman    #ALLOWANCE
#
# Partition:
#SBATCH --partition=savio3
#
# Nodes to use
#SBATCH --nodes=1
#
# Number of tasks needed for use case (example):
#SBATCH --ntasks=1
#
# Processors per task:
#SBATCH --cpus-per-task=40
#
# Wall clock limit:
#SBATCH --time=3-0
#
#SBATCH --mail-user=kozakk@berkeley.edu
#SBATCH --mail-type=ALL
#
## Command(s) to run (example):
module load bio/samtools/1.4
module load fastqc/0.11.7

#check file input names, adjust as needed, e.g.
#rename '_001' '' *.fastq.gz
#rename '-' '.' *.fastq.gz
#rename '_' '.' *.fastq.gz

#path to the reference genome
reference='/global/scratch/users/kozakk/CCGP/ref/Perognathus'

#Evaluate and clean the PE reads
#fastQC with 20 cores b/c there are only 20 samples
fastqc -t 20 --noextract --nogroup *.fastq.gz

while read sample; do
  /global/scratch/users/kozakk/bin/fastp -w 40  -i $sample.R1.fastq.gz -I $sample.R2.fastq.gz -o $sample.R1.fastp.fastq.gz -O $sample.R2.fastp.fastq.gz --failed_out $sample.fastp_failed_reads \
  --dedup \
  --average_qual 10 \
  --length_required 50 \
  --detect_adapter_for_pe --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
  mv fastp.html $sample.fastp.html
done < amplus_new_samplus.txt

#Index the reference
/global/scratch/users/kozakk/bin/bwamem2/bwa-mem2 index $reference

#Map
while read sample; do
echo $sample
  /global/scratch/users/kozakk/bin/bwamem2/bwa-mem2 mem -t 40 \
   $reference \
   $sample.R1.fastp.fastq.gz \
   $sample.R2.fastp.fastq.gz | samtools sort -O BAM -o $sample.sorted.bam \
   && samtools flagstat $sample.sorted.bam > $sample.flagstat.log \
done < amplus_new_samplus.txt
