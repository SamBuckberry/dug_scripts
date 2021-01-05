#!/bin/bash

#SBATCH --job-name=submit_dmr.job
#SBATCH --partition=hg19ips
#SBATCH --output=submit_dmr.out
#SBATCH --error=submit_dmr.err

# Add the modules needed for the R script
module add R

Rscript ./run_dmrseq_dug.R \
/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/workflow/call_dmr/dmr_manifest_test.csv \
/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/dmrs/ \
test_data



