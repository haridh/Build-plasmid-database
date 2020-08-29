suppressMessages(library(fuzzyjoin, quietly = TRUE))
suppressMessages(library(tibble, quietly = TRUE))
suppressMessages(library(dplyr, quietly = TRUE))

# Make a list of ETmapper output files
ET_file_list <- list.files("./hits")
for (ET_file in ET_file_list) {
	ET_file_sub <- paste("hits_subset/Subset", ET_file, sep = "_")
	ET_path <- paste("hits/", ET_file, sep = "")
	print (ET_path)

# For each hit file filter out lines that contain background plasmid ie pHLL and read the file..
	system(paste("grep -v \"pHLL\" ", ET_path, "| awk -F \"\t\" '{print $3\"\t\"$26}' | sort | uniq >", ET_file_sub))
	ET_bar <- read.table(ET_file_sub, sep = "\t", stringsAsFactor = FALSE, header = FALSE)
	head(ET_bar)

# Make a list of pacbio database files split by promoter contents..
	file_list <- list.files("./pcbio_db_partitions")
	for (file in file_list) {
		final_file <- data.frame("Promoter"="", "TPase"="", "pbio_bar"="", "ET_barcode"="", "ET_Genome"="")
		file_path <- paste("./pcbio_db_partitions", file, sep = "/")
		out_file <- paste(ET_file,file, sep = "_")

# Read each pacbio promoter database file and further split into five parts and process each of the five parts one by one. Why?? Steps below cant handle large files.
		pac <- read.table(file_path, sep = "\t", stringsAsFactor = FALSE)
		file_len <- nrow(pac)
		part_size <- round(file_len/5)

# Processing first 4 out of 5 parts

		for (part in c(1:4)) {
			part_stop <- part*part_size
			part_start <- part_stop - part_size
			cat("Merging with: ",file,part_start,"-",part_stop, "out of", file_len, "\n")
			pac1 <- pac[, c(1,2,4)]
			pac1 <- pac1[c(part_start:part_stop),]
			names(pac1)[1] <- paste("Promoter")
			names(pac1)[2] <- paste("TPase")
			names(pac1)[3] <- paste("pbio_bar")
			names(ET_bar)[1] <- paste("ET_Genome")
			names(ET_bar)[2] <- paste("ET_barcode")

# Merge the ETmapper hits file and pacbio part database based on match between pbio and ET barcodes. Allowing mismatches max of 2 dis away.
			k <- as_tibble(ET_bar)
			pac_bar <- pac1 %>% stringdist_inner_join(k, by = c(pbio_bar = "ET_barcode"), max_dist = 2)
			pac_bar <- unique(pac_bar)
			new_lines <- nrow(pac_bar)
			final_file <- rbind(final_file, pac_bar)
			cat("Adding", new_lines, "to", out_file,"....\n")
		}

# Processing 5th part as above. 5th part was expanded to include the lines till the end of file. This is because division by 5 doesnt lead to inclusion of entire file and sometimes exceeds max index.

		part_start <- part_stop
		cat("Merging with: ",file,part_start,"-",file_len, "out of", file_len, "\n")
		pac1 <- pac[, c(1,2,4)]
		pac1 <- pac1[c(part_start:file_len),]
		names(pac1)[1] <- paste("Promoter")
		names(pac1)[2] <- paste("TPase")
		names(pac1)[3] <- paste("pbio_bar")
		names(ET_bar)[1] <- paste("ET_Genome")
		names(ET_bar)[2] <- paste("ET_barcode")
		k <- as_tibble(ET_bar)
		pac_bar <- pac1 %>% stringdist_inner_join(k, by = c(pbio_bar = "ET_barcode"), max_dist = 2)
		pac_bar <- unique(pac_bar)
		new_lines <- nrow(pac_bar)
		final_file <- rbind(final_file, pac_bar)
		cat("Adding", new_lines, "to", out_file,"....\n")
		opfile <- file(out_file)
		write.table(final_file, opfile, sep = "\t", row.names = FALSE, quote = FALSE)
		cat("Done writing file: ",file,"\n")
	}
}
