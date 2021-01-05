library(bsseq)
library(dmrseq)
library(magrittr)
library(stringr)
library(BiocParallel)
library(batchtools)

args = commandArgs(TRUE)

# Path to DMR manifest file. This is a csv file with the column names id, group, rds_path 
manifest <- args[1] # must be absolute path
out_path <- args[2] # must be folder path ending in /
prefix <- args[2] # Out file prefix

#-------- Read the DMR manifest file

# Check that the manifest file exists
stopifnot(file.exists(manifest) == TRUE)

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
   bs_obj <- chrSelectBSseq(BSseq = bs_obj,
                            seqnames = chroms,
                            order = TRUE)
   bs_obj <- strandCollapse(bs_obj)
   return(bs_obj)

}

# Load the data, and sub-select targeted chromosomes
obj_list <- lapply(X = man_dat$rds_path, read_bs_obj,
                   chroms = chrom_list)

# Combine all of the bsseq objects into one
obj_list <- bsseq::combineList(x = obj_list)

# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

# Setup the groups
pData(obj_list)$CellType <- factor(man_dat$group)
colnames(obj_list) <- man_dat$id

#-------- Call DMRs

# This is where we need to the parallel processing.
# See page 12 of https://www.bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
# More in depth documentation at https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf
# and here https://hpc.nih.gov/apps/R.html#biocparallel

# The goal is to setup the `dmrseq` function for parallel processing using multiple nodes

##  Something like this should get the `dmrseq` function below to run on multiple nodes.
## Will need to setup slurm template file for the "template" argument below

#BiocParallel::register(BPPARAM = BatchtoolsParam(workers=24,
#                                                 cluster="slurm",
#                                                 template = "path/to/template/file"))

# The following will use 22 core on one node... Might be too much for a full dataset
# The line below will need to be removed to test the above 
BiocParallel::register(BPPARAM = MulticoreParam(workers = 22))

#!!!!!------- This is the function that will use the parallel processing
# If we run into memory issues, the chrsPerChunk argument could be reduced to a smaller number, like 11

chunks <- length(levels(seqnames(obj_list)))

message(str_c("Chunks = ",chunks))

dmrs <- dmrseq(obj_list, testCovariate = "CellType",
               maxPerms = 10, cutoff = 0.05,
               chrsPerChunk = chunks)

# Output results
out_list <- list(manifest_data = man_dat, dmr_granges = dmrs)

out_file <- str_c(out_path, prefix, "_", uniq_groups[1],
                 "_vs_", uniq_groups[2],
                 "_dmrseq_dmrs.Rds")

saveRDS(object = out_list, file = out_file)

message(str_c("DMR file saved to ", out_file))

