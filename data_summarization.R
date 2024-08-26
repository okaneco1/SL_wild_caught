# Data Summarization

# Load libraries and data
library(tidyverse)

com_mat_sub1 <- read_csv(file = "sub1_community_matrix_organized.csv")
com_mat_sub2 <- read_csv(file = "sub2_community_matrix_organized.csv")

options(scipen = 999)
#---------- Data Organization

# can remove "plate" column from sub2 matrix
com_mat_sub2 <- select(com_mat_sub2, -plate_number)

# trim sub 1 community matrix to wild-caught samples
wild_caught_samples <- c("H22", "C22", "S22", "H23", "C23", "S23")

com_mat_sub1 <- com_mat_sub1 %>%
  filter(str_detect(sample_name, paste(wild_caught_samples, collapse = "|")))


# combine both community matrices
com_mat_full <- bind_rows(com_mat_sub1, com_mat_sub2)

# change any NA values to 0
com_mat_full[is.na(com_mat_full)] <- 0


# remove human reads
com_mat_full <- select(com_mat_full, -Homo_sapiens)

# add total reads column
com_mat_full <- com_mat_full %>%
  mutate(total_reads = rowSums(select(com_mat_full, -sample_id, -replicate, -sample_name)))

# add columns for lake, year, and phase
com_mat_full <- com_mat_full %>%
  mutate(phase = as.character(NA), .after = sample_name) %>%
  mutate(year = as.numeric(NA), .after = sample_name) %>%
  mutate(lake = as.character(NA), .after = sample_name)


for (i in 1:nrow(com_mat_full)) {
  # first add lake names
  if (grepl("C", com_mat_full$sample_name[i])) {
    com_mat_full$lake[i] <- "Champlain"
  } else if (grepl("H", com_mat_full$sample_name[i])) {
    com_mat_full$lake[i] <- "Huron"
  } else if (grepl("S", com_mat_full$sample_name[i])) {
    com_mat_full$lake[i] <- "Superior"
  } else {
    next
  }
  # then add years
  if (grepl("22", com_mat_full$sample_name[i])) {
    com_mat_full$year[i] <- 2022
  } else if (grepl("23", com_mat_full$sample_name[i])) {
    com_mat_full$year[i] <- 2023
  } else {
    next
  }
  # last, add phase
  if (grepl("A", com_mat_full$sample_name[i])) {
    com_mat_full$phase[i] <- "Adult"
  } else if (grepl("P", com_mat_full$sample_name[i])) {
    com_mat_full$phase[i] <- "Parasitic"
  } else{
    next
  }
}
  


# OTU Barcharts for Sequence Reads -------------------

# function for making a bar chart
OTU_barchart <- function(lake_id, year_id, phase_id) {
  # filter to specific samples
  specified_samples <- com_mat_full %>%
    filter(lake == lake_id,
           year == year_id,
           phase == phase_id) %>%
    select(-total_reads)
  # get sample size
  n <- nrow(specified_samples) # 63 total OTUs
  # select OTUs 
  otus <- specified_samples %>%
    select(-sample_name, -sample_id, -replicate, -lake, -phase, -year)
  # sum otus for specified samples
  otus_sum <- data.frame(reads = colSums(otus)) %>%
    rownames_to_column(var = "otu")
  # reverse order
  otus_sum$otu <- factor(otus_sum$otu, levels = rev(otus_sum$otu))
  # create barplot
  #plot <-
    ggplot(otus_sum, aes(x = reads, y = otu)) +
    geom_bar(stat = "identity") +
    labs(title = paste0(lake_id," ", year_id, ", ", phase_id, " (n = ", n, ")"),
         x = "Sequence Reads",
         y = "OTU (Operational Taxonomic Unit")
  # save plot
  #ggsave(filename = paste0("summary_plots/", lake_id, "_", year_id, "_", phase_id, "_summary.png"),
         #plot = plot, width = 8, height = 6)
}


# save all to a PDF
save_all_plots_to_pdf <- function() {
  # set up PDF device
  pdf("summary_plots/OTU_barcharts_summary.pdf")
  
  # generate and save each plot
  OTU_barchart("Champlain", 2022, "Adult")
  OTU_barchart("Champlain", 2023, "Adult")
  OTU_barchart("Huron", 2022, "Adult")
  OTU_barchart("Huron", 2023, "Adult")
  OTU_barchart("Superior", 2022, "Adult")
  OTU_barchart("Superior", 2023, "Adult")
  
  OTU_barchart("Champlain", 2022, "Parasitic")
  OTU_barchart("Champlain", 2023, "Parasitic")
  OTU_barchart("Huron", 2022, "Parasitic")
  OTU_barchart("Huron", 2023, "Parasitic")
  OTU_barchart("Superior", 2022, "Parasitic")
  OTU_barchart("Superior", 2023, "Parasitic")
  
  # close PDF device
  dev.off()
}

# run the function
save_all_plots_to_pdf()
