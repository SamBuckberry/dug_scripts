library(data.table)
library(R.utils)
library(Rsamtools)
library(BSgenome.Hsapiens.UCSC.hg19)
library(magrittr)
library(stringr)

# Do not output scientific notation
options(scipen=999)

# initial garbage collection
gc()

# Set command line arguments
args <- commandArgs(TRUE)

#fastq <- str_c(id, ".fastq.gz")
bam <- args[1]
bam_dedup <- args[2]
CGmap <- args[3]
id <- str_replace(bam, ".bam", "")

### Count the reads in the fastq file
#lines <- system(command = paste("zcat", fastq, "| wc -l"), intern = TRUE)
#reads <- as.numeric(lines) / 4

### Count the reads in the bam file
bam_counts <- countBam(bam)[1,6]

### Count the reads in the deduplicated bam file
bam_counts_dedup <- countBam(bam_dedup)[1,6]

# Get map stats
#map_pc <- (bam_counts / reads) * 100
dedup_pc <- (1 - (bam_counts_dedup / bam_counts)) * 100


# Calc non-conversion stats
dat <- data.table::fread(input = paste('gzip -dc ', CGmap), select = c(1,2,3,5,7,8),
             col.names = c("chr", "base", "position", "context", 
                           "C_reads", "CT_reads"))

chrL <- dat[dat$chr == "chrL", ,drop=TRUE]

# chrL_stat <- round(data.frame(non_conv_rate_pc = (sum(chrL$C_reads)/sum(chrL$CT_reads))*100,
#                          positions_covered = nrow(chrL),
#                          median_x_coverage = median(chrL$CT_reads),
#                          mC_min_pc = min((chrL$C_reads / chrL$CT_reads)*100),
#                          mC_max_pc = max((chrL$C_reads / chrL$CT_reads)*100),
#                          mean_mC_level_pc = mean((chrL$C_reads / chrL$CT_reads)*100), 
#                          mC_std_dev = sd((chrL$C_reads / chrL$CT_reads)*100)), digits = 2)

# Calc autosomal stats
regions <- GRanges(seqnames = seqnames(BSgenome.Hsapiens.UCSC.hg19),
                   ranges = IRanges(start = 1, end = seqlengths(BSgenome.Hsapiens.UCSC.hg19)))

chroms <- paste("chr", c(1:22), sep = "")
regions <- regions[seqnames(regions) %in% chroms]

nuc_counts <- getSeq(BSgenome.Hsapiens.UCSC.hg19, regions) %>%
        alphabetFrequency() %>% data.frame()

c_count <- sum(nuc_counts$C) + sum(nuc_counts$G)

dat <- dat[dat$chr %in% chroms, ,drop=TRUE]

covered_c <- nrow(dat) / c_count

c_summary <- summary(dat$CT_reads)
names(c_summary) <- str_c("Nuclear_genome_", names(c_summary)) %>%
        str_replace(pattern = " ", "_")

c_var <- var(dat$CT_reads)

df <- data.frame(Library = id,
    #Reads = reads,
    Mapped = bam_counts,
    #Map_pc = map_pc,
                 Mapped_dedup=bam_counts_dedup, Duplicate_pc = dedup_pc,
                 chrL_non_conv_pc = (sum(chrL$C_reads)/sum(chrL$CT_reads))*100,
                 chrL_positions_covered = nrow(chrL),
                 chrL_median_coverage = median(chrL$CT_reads),
                 Nuclear_genome_covered_c_pc = covered_c * 100,
                 Nuclear_genome_min = as.numeric(c_summary[1]),
                 Nuclear_genome_1stQu = as.numeric(c_summary[2]),
                 Nuclear_genome_median = as.numeric(c_summary[3]),
                 Nuclear_genome_mean = as.numeric(c_summary[4]),
                 Nuclear_genome_3rdQu = as.numeric(c_summary[5]),
                 Nuclear_genome_max = as.numeric(c_summary[6]),
                 Nuclear_genome_cov_variance = c_var)
                 
# transpose
df <- t(df)

write.table(x = df, file = str_c(id, "_map_stats.txt"), quote = FALSE, sep = "\t",
            row.names = TRUE, col.names = FALSE)
