#!/bin/bash
#rj name=map_pe_dug_pt3 queue=hg19ips features=knl,centos7 schema=input.schema

#MAPPING PIPELINE JOB: PART 3, merged mapped bams and clean up

# Set up environment and common parameters
. map_pe_dug_common.sh

# Set base name for output
baseR1=$(basename "$readR1" $suffixFastq)
baseR2=$(basename "$readR2" $suffixFastq)

# DUG Configuration
pathJob="$pathWorking/$baseR1"
[[ -d $pathJob ]] || mkdir -p "$pathJob"
cd "$pathJob"

# Log active variables
logvars cores readR1 readR2 pathJob baseR1

# Check required input files exist
fileExists "$baseR1"_merged_reads.bam
fileExists "$baseR1"_pairs.bam

msg "==== Merge all bam files ==="
sambamba merge -t "$cores" "$baseR1".bam "$baseR1"_merged_reads.bam "$baseR1"_pairs.bam 

msg "==== Create md5 sum for output bam file ===="
md5sum "$baseR1".bam > "$baseR1".bam.md5 

msg "==== Clean up the temp files ===="
rm 00*."$baseR1"* \
00*."$baseR2"* \
"$baseR1"_trimmed.fastq \
"$baseR2"_trimmed.fastq \
"$baseR1"_merged.fastq \
"$baseR1"_unmerged.fastq \
"$baseR2"_unmerged.fastq \
"$baseR1"_files \
"$baseR2"_files \
"$baseR1"_out \
"$baseR1"_map_manifest \
"$baseR1"_merged_read_bams \
"$baseR1"_merged_fq

msg "=== Moving results to: $pathData/$baseR1"
cd "$pathData" && 
mv "$pathJob" ./
ls "$pathData/$baseR1"

msg "==== part 3 completed ===="