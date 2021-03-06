#!/bin/bash
#rj name=map_pe_dug_pt2b queue=hg19ips features=knl,centos7 schema=input.schema

#MAPPING PIPELINE JOB: PART 2B, map merged read pairs and convert to bam

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
fileExists "$baseR1"_merged.fastq

# Here is where we split the files for run time optimisation.

msg "==== File split pairs and generate reports ===="
fastp -A --split="$files" --split_prefix_digits=4 --thread="$files" \
--in1="$baseR1"_merged.fastq \
--out1="$baseR1"_merged_split.fastq \
--json "$baseR1"_merged_report.json --html "$baseR1"_merged_report.html 

msg "====== Map the merged pairs ======"
msg "==== Setup the trimmed file manifest for mapping ===="
ls 00*."$baseR1"_merged_split.fastq > "$baseR1"_merged_fq 
ls 00*."$baseR1"_merged_split.fastq | sed 's/.fastq/.bam/g' > "$baseR1"_merged_read_bams 

msg "==== Paired-end alignment ===="
paste "$baseR1"_merged_fq "$baseR1"_merged_read_bams | parallel -j"$files" --colsep="\t" bs_seeker2-align.py \
--path=$bowtie2Path \
--aligner=bowtie2 --bt2--end-to-end --bt2-p "$bt_cores" -e 400 \
--temp_dir="$pathScratch" \
-i {1} -o {2} \
-d "$indexBS" \
-g "$genome" 

msg "====  Create a merged log file ===="
cat 00*."$baseR1"_merged_split.bam.bs_seeker2_log > "$baseR1"_merged_pairs_bs_seeker2.log 

msg "==== Sort the output bam files ===="
while read i; do sambamba sort -t "$cores" $i; done < "$baseR1"_merged_read_bams 

msg "==== Merge the single end map bam files ===="
sambamba merge -t "$cores" "$baseR1"_merged_reads.bam 00*."$baseR1"_merged_split.sorted.bam 

msg "=== Clean up scratch ==="
ls "$pathScratch/"*/* >/dev/null 2>&1 && rm "$pathScratch/"*/*

msg "==== part 2b merged workflow completed ===="
