#!/bin/bash
#rj name=bam_merge queue=hg19ips features=knl,centos7 schema=input.schema

# Activate modules in team area
module use /p9/mcc_hg19ips/sw/modules
module add bashutils
module add sambamba

## Set the output folder
pathWorking=/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/aligned

lib="$library"

bam_list="$file_list"

bam_files=$(cat "$bam_list")

# Create folder for output
mkdir "$out_folder"/"$lib"

out_file="$out_folder"/"$lib"/"$lib"_all_merged.bam

# Copy the bam list to the output folder for record keeping
cp "$bam_list" "$out_folder"/"$lib"

# Run the file merging

sambamba merge -p -t 26 "$out_file" $bam_files
