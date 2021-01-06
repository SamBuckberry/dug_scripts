#!/usr/bin/env Rscript

library(bsseq)
library(magrittr)
library(stringr)
library(data.table)
library(R.utils)
library(rtracklayer)
library(GenomicRanges)

# Load library for hg19. This will need to change for other genome assemblies
library(BSgenome.Hsapiens.UCSC.hg19)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

# Path to a CGmap file from BSseeker2
path <- args[1]

# Dinucleotide context
context <- args[2]

# Out path
file_ext <- str_c("_", context, ".bigwig")
out_path <- str_replace(string = path, pattern = ".CGmap.gz$", replacement = file_ext)

#---- Check inputs
stopifnot(file.exists(path))

# Ensure context is upper case
context <- casefold(context, upper = TRUE)

# Check context is correct
stopifnot(context %in% c("CG", "CA", "CC", "CT", "CH"))

# Function to create bigwig file
make_bsseq_obj <- function(path, out_path, context="CG"){
    
    message(str_c("Reading ", path, "..."))
    
    # Read the data
    dat <- data.table::fread(path, header = FALSE,
                             select = c(1,2,3,5,7,8), sep = "\t",
                             col.names = c("chr", "base", "position",
                                           "diContext", "C_reads", "CT_reads"))
    
    # subset data by context
    message(str_c("Subsetting data for dinucleotide context ", context, "..."))
    
    keep <- NA
    
    if (context == "CH"){
        keep <- dat$diContext != "CG"
    } else {
        keep <- dat$diContext == context
    }

    dat <- subset(dat, keep)
    
    message("Calculating methylation levels...")
    
    # Add strand info
    dat$base <- ifelse(test = dat$base == "C", yes = "+", no = "-")
    
    #Load data into BSseq object
    bs_obj <- bsseq::BSseq(chr = dat$chr, pos=dat$position,
                           M = matrix(dat$C_reads),
                           Cov = matrix(dat$CT_reads))
    
    # Add strand to bsseq object
    bsseq::strand(bs_obj) <- dat$base
    
    # Collapse strand if context is CG.
    # This helps reduce file size through aggregating symmetrical CGs
    if (context=="CG"){
        message("Collapsing strands for CG context...")
        bs_obj <- bsseq::strandCollapse(bs_obj)
    }
    
    # Convert to Granges
    gr <- bs_obj@rowRanges
    
    gr$score <- bsseq::getMeth(BSseq = bs_obj, type = "raw", what = "perBase")
    
    gr <- gr[!seqnames(gr) %in% c("chrM", "chrL")]
    
    gr <- dropSeqlevels(gr, value = c("chrM", "chrL"))
    
    # Clean up to free up memory
    rm(bs_obj)
    
    #Add seqinfo

    hg19_info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
    info <- intersect(seqinfo(gr), hg19_info)    
    seqinfo(gr) <- info

    # Write the bigwig file
    message(str_c("Saving ", out_path))
    rtracklayer::export.bw(object = gr, con = out_path)

    message("Done!")
}

# Execute the function with specified arguments
make_bsseq_obj(path = path, out_path = out_path, context = context)