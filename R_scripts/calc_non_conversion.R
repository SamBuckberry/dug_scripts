library(data.table)
library(R.utils)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

CGmap <- args[1]

# Create the temp directory for uncompressed CGmap file
tmpFile <- tempfile(tmpdir = tempdir(), fileext = ".CGmap.tmp")

#gunzip(CGmap, remove=FALSE, destname = tmpFile)
command <- paste("zcat ", CGmap, " > ", tmpFile, sep = "")
system(command)

# Read the file
dat <- fread(input = tmpFile, sep = "\t", select = c(1,2,3,5,7,8),
             col.names = c("chr", "base", "position", "context", 
                           "C_reads", "CT_reads"))

dat <- dat[dat$chr == "chrL", ,drop=TRUE]


stat <- round(data.frame(non_conv_rate_pc = (sum(dat$C_reads)/sum(dat$CT_reads))*100,
          positions_covered = nrow(dat),
          median_x_coverage = median(dat$CT_reads),
          mC_min_pc = min((dat$C_reads / dat$CT_reads)*100),
          mC_max_pc = max((dat$C_reads / dat$CT_reads)*100),
          mean_mC_level_pc = mean((dat$C_reads / dat$CT_reads)*100), 
          mC_std_dev = sd((dat$C_reads / dat$CT_reads)*100)), digits = 2)

write.table(t(stat), file = paste(CGmap, ".non_conversion_rate.txt", sep = ""),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
