#!/bin/bash
#rj name=map_pe_dug_pt1 queue=hg19ips features=knl,centos7 schema=input.schema

#MAPPING PIPELINE JOB: PART 1, merge overlapping read pairs

# Set up environment and common parameters
. map_pe_dug_common.sh

# Fastq input files from schema
fileExists $readR1
fileExists $readR2

# Set base name for output
baseR1=$(basename "$readR1" $suffixFastq)
baseR2=$(basename "$readR2" $suffixFastq)

# DUG Configuration
pathJob="$pathWorking/$baseR1"
[[ -d $pathJob ]] || mkdir -p "$pathJob"
cd "$pathJob"

# Log active variables
logvars cores readR1 readR2 baseR1 known_adapter pathJob

msg "==== Pre-trim reads report ===="
fastp -Q -A -G -L --in1="$readR1" --in2="$readR2" --thread=16 \
--html="$baseR1"_fastp_report.html --json="$baseR1"_fastp_report.json 

msg "==== Remove adapters using bbduk script in the bbmap toolkit ===="
bbduk.sh \
in="$readR1" in2="$readR2" \
out="$baseR1"_trimmed.fastq out2="$baseR2"_trimmed.fastq \
literal="$known_adapter" threads="$cores" overwrite=t -Xmx60g \
overwrite=true 

msg "==== Merge overlapping read pairs using the bbmerge script in the bbmap toolkit ===="
bbmerge.sh in1="$baseR1"_trimmed.fastq in2="$baseR2"_trimmed.fastq qtrim=r \
out="$baseR1"_merged.fastq outu1="$baseR1"_unmerged.fastq outu2="$baseR2"_unmerged.fastq \
-Xmx20g 

msg "==== part 1 completed ===="