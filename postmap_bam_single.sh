#!/bin/bash

# Activate modules in team area
module add hpc
module use /p9/mcc_hg19ips/sw/modules
module add bashutils
module add sambamba

## Set the output folder
out_folder=/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/aligned

# input bam list
file_list=$1

base=$(basename "$file_list" | sed 's/_bam_file_list.txt//g')

bam_files=$(cat "$file_list")

# Create folder for output
mkdir "$out_folder"/"$base"

out_file="$out_folder"/"$base"/"$base"_all_merged.bam

# Copy the bam list to the output folder for record keeping
cp "$file_list" "$out_folder"/"$base"

# Run the file merging
sambamba merge -p -t 26 "$out_file" $bam_files
