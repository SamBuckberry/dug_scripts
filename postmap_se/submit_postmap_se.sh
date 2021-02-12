#!/bin/bash
#rj name=postmap_se queue=hg19ips features=knl,centos7,fastio logdir=logs schema=input.schema runtime=24

#immediately exit if any command has a non-zero exit status
set -e
set -o pipefail

module add hpc
module use /p9/mcc_hg19ips/sw/modules
module add bashutils

sh ./postmap_wgbs_se.sh $bam
