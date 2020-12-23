library(bsseq)
library(dmrseq)
library(magrittr)
library(stringr)
library(BiocParallel)
library(data.table)


#--------Call DMRs

# This is where we need to the parallel processing. 
# See page 12 of https://www.bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
# More in depth documentation at https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf
# and here https://hpc.nih.gov/apps/R.html#biocparallel

library(BiocParallel)
library(Rmpi)
library(batchtools)

tmpl <- system.file(package="batchtools", "templates", "slurm-simple.tmpl")
noquote(readLines(tmpl))




FUN <- function(i) system("hostname", intern=TRUE)
library(BiocParallel)

## register SLURM cluster instructions from the template file
param <- BatchtoolsParam(workers=5, cluster="slurm", template=tmpl)
register(param)
## do work
xx <- bplapply(1:100, FUN)
table(unlist(xx))




BiocParallel::register(BPPARAM = BatchtoolsParam(workers=100, cluster="slurm"))

FUN <- function(i) system("hostname", intern=TRUE)

xx <- bplapply(1:100, FUN)



library(BiocParallel)
library(Rmpi)
FUN <- function(i) system("hostname", intern=TRUE)

param <- SnowParam(mpi.universe.size() - 1, "MPI")
register(param)



# Load the data
obj_list <- lapply(X = obj_fls, read_bs_obj)
obj_list <- bsseq::combineList(x = obj_list)

# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

# Setup the groups
pData(obj_list)$CellType <- factor(c(rep("ESC", times=5), rep("Naive2Primed", times=5)))
pData(obj_list)$Replicate <- c(1:5,1:5)




obj_list <- readRDS("/home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/processed_data/dmrseq_dmrs/esc_vs_n2p_BSseq_obj.rds")

primed_n2p_dmrs <- dmrseq(obj_list, testCovariate = "CellType",
                          bpSpan = 500,
                          maxGap = 500,
                          maxPerms = 10)

saveRDS(object = primed_n2p_dmrs,
        file = "/home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/processed_data/dmrseq_dmrs/esc_vs_n2p_dmrseq_DMRs.rds")

