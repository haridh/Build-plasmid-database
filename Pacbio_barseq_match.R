install.packages("magrittr")
install.packages("fuzzyjoin")
suppressMessages(library(fuzzyjoin, quietly = TRUE))
suppressMessages(library(tibble, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))
args <- commandArgs(trailingOnly = TRUE)

# Insert base path
base_path <- "/groups/banfield/users/haridh/Ben_Brady_collab/MagicPool/BARseq"

# Two arguments needed. The first is a list of barcodes identified from BARseq and 2nd is an output prefix
# This script works just like the Pacbio_ETseq_match.R script except for BARseq files.
barseq_file <- args[1]
barseq_pfx <- args[2]
barseq_file_path <- paste(base_path, barseq_file, sep = "/")
bar <- read.table(barseq_file_path, sep = "\t", stringsAsFactor = FALSE)
pacbio_files <- paste(base_path, "pacbio_split", sep = "/")
file_list <- list.files(pacbio_files)
for (file in file_list) {
	final_file <- data.frame("Promoter"="", "TPase"="", "pbio_bar"="", "barseq"="")
	file_path <- paste("./pacbio_split", file, sep = "/")
	out_pfx <- paste(base_path, barseq_pfx, sep = "/")
	out_file <- paste(out_pfx, file, sep = "_")
	pac <- read.table(file_path, sep = "\t", stringsAsFactor = FALSE)
	file_len <- nrow(pac)
	part_size <- round(file_len/5)
	for (part in c(1:4)) {
		part_stop <- part*part_size
		part_start <- part_stop - part_size
		cat("Merging with: ",file,part_start,"-",part_stop, "out of", file_len, "\n")
		pac1 <- pac[, c(1,2,4)]
		pac1 <- pac1[c(part_start:part_stop),]
		names(pac1)[1] <- paste("Promoter")
		names(pac1)[2] <- paste("TPase")
		names(pac1)[3] <- paste("pbio_bar")
		k <- as_tibble(bar$V1)
		k1 <- rename(k, barseq = value)
		pac_bar <- pac1 %>% stringdist_inner_join(k1, by = c(pbio_bar = "barseq"), max_dist = 2)
		pac_bar <- unique(pac_bar)
		new_lines <- nrow(pac_bar)
		final_file <- rbind(final_file, pac_bar)
		cat("Adding", new_lines, "to", out_file,"....\n")
	}
	part_start <- part_stop
	cat("Merging with: ",file,part_start,"-",file_len, "out of", file_len, "\n")
	pac1 <- pac[, c(1,2,4)]
	pac1 <- pac1[c(part_start:file_len),]
	names(pac1)[1] <- paste("Promoter")
	names(pac1)[2] <- paste("TPase")
	names(pac1)[3] <- paste("pbio_bar")
	k <- as_tibble(bar$V1)
	k1 <- rename(k, barseq = value)
	pac_bar <- pac1 %>% stringdist_inner_join(k1, by = c(pbio_bar = "barseq"), max_dist = 2)
	pac_bar <- unique(pac_bar)
	new_lines <- nrow(pac_bar)
	final_file <- rbind(final_file, pac_bar)
	cat("Adding", new_lines, "to", out_file,"....\n")
	opfile <- file(out_file)
	write.table(final_file, opfile, sep = "\t", row.names = FALSE, quote = FALSE)
	cat("Done writing file: ",file,"\n")
}
