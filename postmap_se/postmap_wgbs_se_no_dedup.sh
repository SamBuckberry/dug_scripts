#!/bin/bash

# Activate modules in team area
module use /p9/mcc_hg19ips/sw/modules
module add bashutils
module add samtools
module add sambamba
module add parallel

#immediately exit if any command has a non-zero exit status
set -euo pipefail

# Argument 1: Path to a sorted BAM file
file="$1"

fileExists $file

# Allocated cores
cores="24"

# number of cores for sambamba markdup processing
sambamba_markdup_cores="16"

# Path to the directory containing the indexed reference genome
indexBS="/p9/mcc_hg19ips/data/indexes/hg19_bsseeker_bt2_index/" # For human hg19

# Name of the genome file to use. This file must be in the "indexBS" path above in fasta format (must be uncompressed)
genome="hg19_L_PhiX.fa" # For human hg19 and chrL

prefix=$(basename "$file" .bam)

# Setup the job directory
pathJob="/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/postmap_se/$prefix"
[[ -d $pathJob ]] || mkdir -p "$pathJob"
cd "$pathJob"

# Setup up the working directory
pathWorking="$pathJob/working"
[[ -d $pathWorking ]] || mkdir -p "$pathWorking"

# Setup the scratch directory
pathScratch="$pathWorking/000scratch/$SLURM_JOB_ID"
[[ -d $pathScratch ]] || mkdir -p "$pathScratch"

logvars file cores sambamba_markdup_cores indexBS genome prefix

#msg "==== Remove the PCR duplicates and save output file===="
#sambamba markdup -r -t $sambamba_markdup_cores -p --tmpdir=$pathScratch "$file" $pathWorking/"$prefix".dedup.bam

msg "==== Get the chromosome/contig id's for splitting. grep -V avoids bug with samtools output===="
#samtools idxstats $pathWorking/"$prefix".dedup.bam | cut -f 1 | grep -v \* > $pathWorking/"$prefix".chromNames

samtools idxstats "$file" | cut -f 1 | grep -v \* > $pathWorking/"$prefix".chromNames

msg "==== Split the bamfile into chromosome specific files ===="
# Split the bamfile into chromosome specific files
cat $pathWorking/"$prefix".chromNames | parallel -j$cores samtools view -b -o $pathWorking/"$prefix".{1}.temp.bam "$file" {1}

msg "==== Remove phiX bam ===="
# Remove lambda bam
rm $pathWorking/"$prefix".gi_9626372_ref_NC_001422_1_.temp.bam

msg "==== call mC levels===="
module add cgmaptools
parallel -j$cores CGmapFromBAM -b {} -g "$indexBS/$genome" -o {.} ::: $pathWorking/"$prefix".chr*.temp.bam

msg "==== Recombine the output files and clean up===="
# Recombine the output files and clean up

ls $pathWorking/"$prefix".*.temp*.ATCGmap.gz | xargs zcat > "$prefix".ATCGmap && pigz -f -p $cores "$prefix".ATCGmap &&

ls $pathWorking/"$prefix".*.temp*.CGmap.gz | xargs zcat > "$prefix".CGmap && pigz -f -p $cores "$prefix".CGmap &&

msg "==== Cleanup files ===="
rm $pathWorking/"$prefix".*wig.gz $pathWorking/"$prefix"*temp* $pathWorking/"$prefix".chromNames \

msg "=== Clean up scratch ==="
rm -r "$pathScratch/"
