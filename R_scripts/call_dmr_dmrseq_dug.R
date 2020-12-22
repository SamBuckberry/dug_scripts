library(bsseq)
library(dmrseq)
library(magrittr)
library(stringr)
library(BiocParallel)
library(data.table)


# This is where we need to the parallel processing. 
# See page 12 of https://www.bioconductor.org/packages/release/bioc/vignettes/BiocParallel/inst/doc/Introduction_To_BiocParallel.pdf
# More in depth documentation at https://bioconductor.org/packages/release/bioc/manuals/BiocParallel/man/BiocParallel.pdf
# and here https://hpc.nih.gov/apps/R.html#biocparallel
BiocParallel::register()

obj_fls <- list.files(path = "/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/processed_data/bsseq_cg",
                      pattern = ".Rds", full.names = TRUE)

# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path){
    bs_obj <- readRDS(file = rds_path)
    bs_obj <- strandCollapse(bs_obj)
    return(bs_obj)
}

#--------Call DMRs

# Load the data
obj_list <- lapply(X = obj_files, read_bs_obj)
obj_list <- bsseq::combineList(x = obj_fls)


pData(obj_list)$CellType <- factor(c(rep("ESC", times=5), rep("Naive2Primed", times=5)))
pData(obj_list)$Replicate <- c(1:5,1:5)

# Remove CpG with no coverage
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(obj_list, type="Cov")==0) == 0)
obj_list <- obj_list[loci.idx, ]

# Save object for easy retreval
saveRDS(object = obj_list,
        file = "/home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/processed_data/dmrseq_dmrs/esc_vs_n2p_BSseq_obj.rds")

obj_list <- readRDS("/home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/processed_data/dmrseq_dmrs/esc_vs_n2p_BSseq_obj.rds")

primed_n2p_dmrs <- dmrseq(obj_list, testCovariate = "CellType",
                          bpSpan = 500,
                          maxGap = 500,
                          maxPerms = 10,
                          chrsPerChunk = 1)

saveRDS(object = primed_n2p_dmrs,
        file = "/home/sbuckberry/working_data_02/polo_project/human_ips/methylCseq/processed_data/dmrseq_dmrs/esc_vs_n2p_dmrseq_DMRs.rds")

