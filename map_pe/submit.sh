#!/bin/bash

#SUBMIT SCRIPT FOR MAPPING PIPELINE
#Usage: ./submit.sh
#1. Update variables in config.sh
#2. Update R1/R2 read pairs in input.schema. Each line contains one pair of read files
#3. Run ./submit.sh
#   Optional arguments: array=(sample#),(sample#),(sample#)... Only run these samples from the schema file
#   Optional arguments: dep=SLURM_JOB_ID           Do not start this job until job# SLURM_JOB_ID is complete

#Expected schema file: (sample#) depTasks=(sample#) readR1=/path/to/R1_reads.fastq.gz readR2=/path/to/R2_reads.fastq.gz
#1 depTasks=1 readR1=/p9/mcc_hg19ips/data/test_1M_reads_R1.fastq.gz readR2=/p9/mcc_hg19ips/data/test_1M_reads_R2.fastq.gz 
#2 depTasks=2 readR1=/p9/mcc_hg19ips/data/test_10M_reads_R1.fastq.gz readR2=/p9/mcc_hg19ips/data/test_10M_reads_R2.fastq.gz

#immediately exit if any command has a non-zero exit status
set -e 
set -o pipefail

module add hpc
module use /p9/mcc_hg19ips/sw/modules
module add bashutils

export PATH=/d/sw/Insight/latest/scripts:$PATH

export fileScriptCommon=map_pe_dug_common.sh
export fileScriptOne=map_pe_dug_pt1.sh
export fileScriptTwoA=map_pe_dug_pt2a.sh
export fileScriptTwoB=map_pe_dug_pt2b.sh
export fileScriptThree=map_pe_dug_pt3.sh

CL_PARAMS="array dep"
. command-line

[[ $array ]] && RJSArray="array=$array"

fileExists $fileScriptCommon
fileExists $fileScriptOne
fileExists $fileScriptTwoA
fileExists $fileScriptTwoB
fileExists $fileScriptThree
fileExists input.schema
fileExists config.sh

msg "=== Submitting map jobs to cluster ==="

depJidsOne=$dep
jidsOne=$(rjs $@ $RJSArray dep=$depJidsOne "$fileScriptOne")
if [[ ! $jidsOne ]]; then
    msg failed to submit $fileScriptOne >&2
    exit 100
fi

msg $fileScriptOne job $jidsOne has been submitted

depJidsTwo=$jidsOne
jidsTwoA=$(rjs $@ $RJSArray dep=$depJidsTwo "$fileScriptTwoA")
if [[ ! $jidsTwoA ]]; then
    msg failed to submit $fileScriptTwoA >&2
    exit 100
fi

msg $fileScriptTwoA job $jidsTwoA has been submitted

jidsTwoB=$(rjs $@ $RJSArray dep=$depJidsTwo "$fileScriptTwoB")
if [[ ! $jidsTwoB ]]; then
    msg failed to submit $fileScriptTwoB >&2
    exit 100
fi

msg $fileScriptTwoB job $jidsTwoB has been submitted

depJidsThree="$jidsTwoA,$jidsTwoB"
jidsThree=$(rjs $@ $RJSArray dep=$depJidsThree "$fileScriptThree")
if [[ ! $jidsThree ]]; then
    msg failed to submit $fileScriptThree >&2
    exit 100
fi

msg $fileScriptThree job $jidsThree has been submitted

