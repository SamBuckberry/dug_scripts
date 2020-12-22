#!/bin/bash

#MAPPING PIPELINE: CONFIGURATION & SETTINGS
#CALLED BY: map_pe_dug_common.sh

# Path to store output final data / workflow results
pathData="/p9/mcc_hg19ips/data/map_pe"

# Path for jobs that have not completed
pathWorking="$pathData/working"

# Path to the directory containing the indexed reference genome
indexBS="/p9/mcc_hg19ips/data/indexes/hg19_bsseeker_bt2_index/" # For human hg19

# Name of the genome file to use. This file must be in the "indexBS" path above
genome="hg19_L_PhiX.fa" # For human hg19 and chrL

known_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC,AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

