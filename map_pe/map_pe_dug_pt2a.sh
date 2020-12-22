#!/bin/bash
#rj name=map_pe_dug_pt2a queue=hg19ips features=knl,centos7 schema=input.schema

#MAPPING PIPELINE JOB: PART 2A, map unmerged read pairs and convert to bam

# Set up environment and common parameters
. map_pe_dug_common.sh

# Set base name for output
baseR1=$(basename "$readR1" $suffixFastq)
baseR2=$(basename "$readR2" $suffixFastq)

# DUG Configuration
pathJob="$pathWorking/$baseR1"
[[ -d $pathJob ]] || mkdir -p "$pathJob"
cd "$pathJob"

pathScratch="$pathJob/000scratch/$SLURM_JOB_ID"
[[ -d $pathScratch ]] || mkdir -p "$pathScratch"

# Log active variables
logvars cores files bt_cores readR1 readR2 baseR1 indexBS genome pathJob

# Check required input files exist
fileExists "$baseR1"_unmerged.fastq
fileExists "$baseR2"_unmerged.fastq 

# Here is where we split the files for run time optimisation.

msg "==== File split pairs and generate reports ===="
fastp -A --split="$files" --split_prefix_digits=4 --thread="$files" \
--in1="$baseR1"_unmerged.fastq --in2="$baseR2"_unmerged.fastq \
--out1="$baseR1"_unmerged_split.fastq --out2="$baseR2"_unmerged_split.fastq \
--json "$baseR1"_unmerged_report.json --html "$baseR1"_unmerged_report.html 

msg "====== Map the unmerged paired-end reads ======"

msg "==== Setup the trimmed file manifest for mapping ===="
ls 00*."$baseR1"_unmerged_split.fastq > "$baseR1"_files 
ls 00*."$baseR2"_unmerged_split.fastq > "$baseR2"_files 
sed 's/.fastq/.bam/g' "$baseR1"_files > "$baseR1"_out 
paste "$baseR1"_files "$baseR2"_files "$baseR1"_out > "$baseR1"_map_manifest 

msg "==== Paired-end alignment ===="
# Paired-end alignment. The --split_line argument has not undergone optimisation testing. Larger number may be more efficient!
cat "$baseR1"_map_manifest | parallel --colsep="\t" -j"$files" bs_seeker2-align.py \
--path=$bowtie2Path \
--aligner=bowtie2 --bt2--end-to-end --bt2-p "$bt_cores" -e 300 -X 2000 \
--temp_dir="$pathScratch" \
--split_line=4000000 \
-1 {1} -2 {2} -o {3} \
-d "$indexBS" \
-g "$genome" 

msg "==== Sort the output bam files ===="
while read i; do sambamba sort -t "$cores" $i; done < "$baseR1"_out 

msg "==== Clip the overlaps in the sorted BAM files using bamUtil ===="
parallel -j"$files" bam clipOverlap --stats --in {} --out {.}_clipped.bam ::: 00*."$baseR1"_unmerged_split.sorted.bam 

msg "==== Merge the bam files ===="
sambamba merge -t "$cores" "$baseR1"_pairs.bam 00*."$baseR1"_unmerged_split.sorted_clipped.bam 

msg "==== Create a merged log file ===="
cat 00*."$baseR1"_unmerged_split.bam.bs_seeker2_log > "$baseR1"_pairs_bs_seeker2.log 

msg "=== Clean up scratch ==="
ls "$pathScratch/"*/* >/dev/null 2>&1 && rm "$pathScratch/"*/*

msg "==== part 2a unmerged workflow completed ===="