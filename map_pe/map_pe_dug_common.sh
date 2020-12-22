#!/bin/bash

#MAPPING PIPELINE: SHARED / COMMON SCRIPT CONFIGURATION
#CALLED BY: map_pe_dug_pt*.sh

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
module add bamUtil
module add bowtie2
module add java64/12.0.1
bowtie2Path=/p9/mcc_hg19ips/sw/bowtie2/bowtie2-2.4.2/bin

#immediately exit if any command has a non-zero exit status
set -euo pipefail

# Number of cores.
cores="272"

# Number of files split into during mapping stage. 6-10 is safest. 
files="10" 

# Magic number for best bowtie thread parameter: 
#   (cores / (files * constant)) rounded down to nearest 4
#   4 -> number of hardware threads per core
constant=2
bt_cores=$(( $cores / $files / $constant / 4 * 4 )) #integer math nearest mult of 4
[[ $bt_cores -gt 0 ]] || bt_cores=1

# Suffix to trim from input files for naming
suffixFastq=".fastq.gz"

. ./config.sh
