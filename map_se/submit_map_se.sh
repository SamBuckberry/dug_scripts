#!/bin/bash
#rj name=map_se queue=hg19ips features=knl,centos7,fastio logdir=logs schema=input.schema runtime=48

#immediately exit if any command has a non-zero exit status
set -e 
set -o pipefail

module add hpc
module use /p9/mcc_hg19ips/sw/modules
module add bashutils

sh ./map_se_bt2_hg19.sh $readR1 $prefix