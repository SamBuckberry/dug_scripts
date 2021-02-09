#!/bin/bash

# Activate modules in team area
module use /p9/mcc_hg19ips/sw/modules
module add bashutils

# Add modules to current environment
module add parallel
module add bs-seeker2
conda deactivate && conda activate /p9/mcc_hg19ips/sw/bs-seeker2/bs-seeker2-2.1.7/bs-seeker2.env
module add samtools
module add fastp
module add bbmap
module add sambamba
module add bowtie2
module add java64/12.0.1
bowtie2Path=/p9/mcc_hg19ips/sw/bowtie2/bowtie2-2.4.2/bin

#immediately exit if any command has a non-zero exit status
set -euo pipefail

# Fastq input
readR1="$1"

# prefix
prefix="$2"

# Check required input files exist
fileExists "$readR1"

msg "=== Setting up directories ==="

# Setup the job directory
pathJob="/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/map_se/$prefix"
[[ -d $pathJob ]] || mkdir -p "$pathJob"
cd "$pathJob"

# Setup up the working directory
pathWorking="$pathJob/working"
[[ -d $pathWorking ]] || mkdir -p "$pathWorking"

# Setup the scratch directory
pathScratch="$pathWorking/000scratch/$SLURM_JOB_ID"
[[ -d $pathScratch ]] || mkdir -p "$pathScratch"

# Path to the directory containing the indexed reference genome
indexBS="/p9/mcc_hg19ips/data/indexes/hg19_bsseeker_bt2_index/" # For human hg19

# Name of the genome file to use. This file must be in the "indexBS" path above
genome="hg19_L_PhiX.fa" # For human hg19 and chrL

# Adapter sequence
rightAdapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

# Log active variables
logvars readR1 prefix pathJob pathWorking pathScratch indexBS genome rightAdapter

msg "=== Adapter and quality trimming ==="
bash bbduk.sh \
in="$readR1" out="$pathWorking"/"$prefix".trimmed.fastq.gz \
literal="$rightAdapter" \
ktrim=r \
mink=3 \
qtrim=rl \
trimq=10 \
minlength=20 \
overwrite=true

##### Split the trimmed fastq file for parallel alignments
msg "=== Split fastq file ==="
fastp -A --disable_quality_filtering --html "$prefix"_fastp.html \
-w 6 -s 7 --in1 "$pathWorking"/"$prefix".trimmed.fastq.gz -o "$pathWorking"/$prefix.split.fastq.gz &&
rm "$pathWorking"/"$prefix".trimmed.fastq.gz

msg "=== Single-end alignment ==="
parallel -j7 bs_seeker2-align.py \
--aligner=bowtie2 --bt2--end-to-end --bt2-p 2 -e 300 \
--temp_dir="$pathScratch" \
-i {} \
-o {.}.bam \
-d "$indexBS" \
-g "$genome"  ::: "$pathWorking"/000[1-7].$prefix.split.fastq.gz

msg "=== Sort bam files ==="
for f in "$pathWorking"/000[1-7].$prefix.split.fastq.bam; do
	sambamba sort -t 24 $f
done

rm "$pathWorking"/000[1-7].$prefix.split.fastq.gz
rm "$pathWorking"/000[1-7].$prefix.split.fastq.bam

msg "=== Merge the bam files ==="
sambamba merge -t 24 "$prefix".bam "$pathWorking"/000[1-7].$prefix.split.fastq.sorted.bam &&
rm "$pathWorking"/000[1-7].$prefix.split.fastq.sorted.bam
rm "$pathWorking"/000[1-7].$prefix.split.fastq.sorted.bam.bai

msg "=== Create a merged log file ==="
cat "$pathWorking"/000[1-7].$prefix.split.fastq.bam.bs_seeker2_log > "$prefix".bs_seeker2.log &&
rm "$pathWorking"/000[1-7].$prefix.split.fastq.bam.bs_seeker2_log

msg "=== Create md5 sum for output bam file ==="
md5sum "$prefix".bam > "$prefix".bam.md5

msg "=== Clean up scratch ==="
ls "$pathScratch/"*/* >/dev/null 2>&1 && rm "$pathScratch/"*/*
