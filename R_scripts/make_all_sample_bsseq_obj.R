library(bsseq)

obj_fls <- list.files(path = "/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/processed_data/bsseq_cg",
                      pattern = ".Rds", full.names = TRUE)

# Read a Bs_seq boject from .Rds file
read_bs_obj <- function(rds_path){
    bs_obj <- readRDS(file = rds_path)
    bs_obj <- strandCollapse(bs_obj)
    return(bs_obj)
}

# Save object for easy retreval
saveRDS(object = obj_list,
        file = "/d/home/hg19ips/hg19ips_samb/mcc_hg19ips/data/processed_data/bsseq_cg/all_samples_BSseq_obj.rds")
