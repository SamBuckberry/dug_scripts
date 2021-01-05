library(bsseq)
library(dmrseq)
library(magrittr)
library(stringr)
library(BiocParallel)
library(data.table)
#library(Rmpi)
#library(batchtools)

args = commandArgs(TRUE)

# Path to DMR manifest file. This is a csv file with the column names id, group, rds_path 
manifest <- args[1]
out_path <- args[2]
prefix <- args[3]

#-------- Read the DMR manifest file

# Check that the manifest file exists

(!file.exists(manifest)){
    cat('Cannot find', manifest, 'exiting!\n')
    stop()
}

# Read the manifest file
man_dat <- read.table(manifest, header = TRUE, sep = ",",
                      colClasses = "character")

# Check that Rds files exist
rds_check <- lapply(man_dat$rds_path, file.exists) %>% unlist()
stopifnot(all(rds_check) == TRUE)

# Get groups and test that only two groups exist
uniq_groups <- man_dat$group %>% unique()
stopifnot(length(uniq_groups) == 2)

# Check that library ids are unique
stopifnot(length(unique(man_dat$id)) == nrow(man_dat))

# Chromosomes to test for DMRs
chrom_list <- str_c("chr", 1:22)

# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path, chroms){

   bs_obj <- readRDS(file = rds_path)
   bs_obj <- strandCollapse(bs_obj)
   bs_obj <- chrSelectBSseq(BSseq = bs_obj,
                            seqnames = chroms,
                            order = TRUE)
   return(bs_obj)

}

# Load the data, and sub-select targeted chromosomes
obj_list <- lapply(X = obj_fls, read_bs_obj, chroms = str_c("chr", 1:22))
obj_list <- bsseq::combineList(x = obj_list)

# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

# Setup the groups
pData(obj_list)$CellType <- factor(man_dat$group)

# Ensure metadata and data are in same order
stopifnot(all(str_sub(md$cg_bsseq, start = 1, end = 6) == str_sub(colnames(obj_list), start = 1, end = 6)))

colnames(obj_list) <- man_dat$id

# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

#-------- Call DMRs

# This is where we need to the parallel processing.
# See page 12 of https://www.bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
# More in depth documentation at https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf
# and here https://hpc.nih.gov/apps/R.html#biocparallel

# The goal is to use one node per chromosome for parallel processing

BiocParallel::register(BPPARAM = BatchtoolsParam(workers=length(chrom_list), cluster="slurm"))

#!!!!!------- This is the function that will use the parallel processing
dmrs <- dmrseq(obj_list, testCovariate = "CellType")

out_list <- list(obj_list, dmrs)

# Output results
dmr_out <- str_c(out_path, prefix, "dmrseq_dmrs.Rds")
saveRDS(object = dmrs, file = dmr_out)
message(str_c("DMR file saved to ", dmr_out))


savehistory(file = str_c(out_path, prefix, "dmrseq_dmrs.log"))

sink(str_c(out_path, prefix, "dmrseq_dmrs.log"), split = TRUE)

