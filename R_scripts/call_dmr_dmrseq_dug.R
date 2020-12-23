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

tmpl <- "/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/workflow/R_scripts/slurm-dmr-template.tmpl"

FUN <- function(i) system("hostname", intern=TRUE)

## register SLURM cluster instructions from the template file
param <- BatchtoolsParam(workers=50, cluster="slurm", template=tmpl)
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


md <- read.csv(file = "/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/workflow/resources/sample_metadata.csv")
colnames(md)

obj_list <- readRDS("/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/processed_data/bsseq_cg/all_samples_BSseq_obj.rds")

library(stringr)

# Ensure metadata and data are in same order
stopifnot(all(str_sub(md$cg_bsseq, start = 1, end = 6) == str_sub(colnames(obj_list), start = 1, end = 6)))

colnames(obj_list) <- md$library


# Subset to samples for testing
keep_ids <- c("Primed", "TNT")
obj_list <- obj_list[ ,md$group %in% keep_ids]


md <- md[md$group %in% keep_ids, ]
stopifnot(all(colnames(obj_list) == md$library))



# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

# Setup the groups
pData(obj_list)$CellType <- md$group

primed_tnt_dmrs <- dmrseq(obj_list, testCovariate = "CellType",
                          bpSpan = 500,
                          maxGap = 500,
                          maxPerms = 10)

saveRDS(object = primed_tnt_dmrs,
        file = "/home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/processed_data/dmrseq_dmrs/esc_vs_n2p_dmrseq_DMRs.rds")

