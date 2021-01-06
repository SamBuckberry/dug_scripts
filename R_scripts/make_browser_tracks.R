library(data.table)
library(R.utils)
library(magrittr)
library(stringr)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

# Arguments

# Path to a CGmap file from BSseeker2
path <- args[1]

# Prefix for output bigwig files
prefix <- args[2]

# Load the CGmap data
loadData <- function(path) {
        header <- c("chr", "base", "position", "triContext",
                    "diContext", "mC", "C_reads", "CT_reads")
        
        dat <- data.table::fread(path,
          header = TRUE, col.names = header)

        return(dat)
}
dat <- loadData(path)

# Drop contigs not suitable for browser
dat <- dat[!dat$chr %in% c("chrL", "gi_9626372_ref_NC_001422_1_"), ]

# Get the CG methylation levels
aggregate_cg <- function(dat){
        
        # Subset to CG context only
        dat <- dat[dat$diContext == "CG", ]
        
        # Set up the locus id
        dat$locus <- NA
        
        # Get a locus id relative to forward strand
        dat$locus <- ifelse(test = dat$base == "G", 
                            yes = paste(dat$chr, dat$position - 1, sep = ":"),
                            no = paste(dat$chr, dat$position, sep = ":"))

        # Sum the read counts for + and - strand
        combined <- dat[, lapply(.SD, sum), by=.(locus), .SDcols=c("C_reads", "CT_reads")]

        # Set the position
        pos <- str_split(string = combined$locus, pattern = ":", n = 2, simplify = TRUE)
        
        # Make into GRanges
        cg_gr <- GRanges(seqnames = pos[ ,1],
                         ranges = IRanges(start = as.numeric(pos[ ,2]),
                                          end = as.numeric(pos[ ,2])),
                         score=combined$C_reads/combined$CT_reads)

        #Add seqinfo
        info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
        info <- info[seqlevels(cg_gr)]                 
        seqinfo(cg_gr) <- info
        
        return(cg_gr)
        
}

message("Aggregating CG data...")
cg_gr <- aggregate_cg(dat)
message("Writing mCG levels bigwig file...")
rtracklayer::export.bw(object = cg_gr, con = str_c(prefix, "_mCG_levels.bigwig"))
rm(cg_gr)
gc()


# make GRanges of dat for next steps
gr <- GRanges(seqnames = dat$chr, ranges = IRanges(start = dat$position, end = dat$position),
              base=dat$base, context=dat$diContext, C_reads=dat$C_reads, CT_reads=dat$CT_reads,
              mC=dat$mC)
rm(dat)
gc()

# Get the CH methylation levels
# make_ch_gr <- function(gr){
        
#         # Subset to CG context only
#         #ch_dat <- dat[dat$diContext != "CG", ]
#         ch_gr <- gr[gr$context != "CG"]
        
#         ch_gr <- GRanges(seqnames = seqnames(ch_gr),
#                          ranges = IRanges(start = start(ch_gr),
#                                           end = end(ch_gr)),
#                          score=ch_gr$mC)
        
#         seqlevels(ch_gr) <- seqnames(BSgenome.Hsapiens.UCSC.hg19)
#         seqinfo(ch_gr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
        
#         return(ch_gr)
        
# }

# message("Calculating mCH levels...")
# ch_gr <- make_ch_gr(gr)
# message("Writing mCH levels bigwig file...")
# rtracklayer::export.bw(object = ch_gr, con = str_c(prefix, "_mCH_levels.bigwig"))
# rm(ch_gr)
# gc()


# Get the C coverage
message("Calculating coverage...")
gr <- GRanges(seqnames = seqnames(gr),
                  ranges = IRanges(start = start(gr), end = end(gr)),
                  score=gr$CT_reads)
seqlevels(gr) <- seqnames(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(gr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
message("Writing C coverage bigiwg file...")
rtracklayer::export.bw(object = gr, con = str_c(prefix, "_C_coverage.bigwig"))




