# Initial Community Matrix Analyses of Submission 2 Sequence Reads

# Libraries
library(tidyverse)
library(readxl)
library(vegan)

# Load in Data
com_mat <- read_excel("submission2_full_community_matrix_3.xlsx")


# Data Organization
colnames(com_mat)[1] <- "sample_id"
com_mat <- com_mat[,-61]

# add column to distinguish replicate 1 and 2
com_mat <- com_mat %>%
  mutate(replicate = NA, .after = sample_id)

for (i in 1:nrow(com_mat)) {
  if (grepl("_R.12S", com_mat$sample_id[i])) {
    com_mat$replicate[i] <- "replicate 2"
  } else {
    com_mat$replicate[i] <- "replicate 1"
  }
}

# add column for sample names (not the id)
com_mat <- com_mat %>%
  mutate(sample_name = NA, .after = replicate)

# need a function for identifying sample names
extract_sample_names2 <- function(sample) {
  # check to see if it is a replicate 2
  if (grepl("_R\\.12S$", sample)) {
    return(str_extract(sample, ".*(?=_R)"))
  # check to see if replicate 1 (ie does not have "_R")
  } else if (!grepl("_R", sample)) {
    return(str_extract(sample, ".*(?=\\.12S)"))
  } else {
    return()
  }
}

# assign sample names to sample name column
for (i in 1:nrow(com_mat)) {
  com_mat$sample_name[i] <- extract_sample_names2(com_mat$sample_id[i])
}


# ------------------------------
# Assigning Plate Numbers to Samples

# initialize column
com_mat <- com_mat %>%
  mutate(plate_number = NA, .after = sample_name)

# initialize an empty list to store the stacked data frames
stacked_plate_list <- list()

# loop over each sheet, saving each cell name to a stacked list
for (i in 1:7) {
  # read the sheet
  plate <- read_excel("./submission 2 PCR plate setup.xlsx",
                      sheet = i, col_names = FALSE,
                      range = "B5:M12") # full plate, adjustments for partial plates done above
  
  # unlists and stores all names for each plate as single column
  stacked_plate_list[[paste0("plate", i, "_stacked")]] <- data.frame(sample_id = unlist(plate, use.names = FALSE))
}

# view list here
stacked_plate_list
# save stacked plates to global environment
list2env(stacked_plate_list, envir = .GlobalEnv)
# filter partial plate
plate7_stacked <- plate7_stacked[!is.na(plate7_stacked$sample_id), ]
plate7_stacked <- data.frame(sample_id = plate7_stacked)

# do the same for the replicate sheets
# set up function
duplicate_plates <- function(plate) {
  #duplicate the plate 
  plate_duplicate <- plate
  # add "_R" to sample id
  plate_duplicate <- paste0(plate_duplicate$sample_id, "_R")
  # return duplicate
  data.frame(sample_id = plate_duplicate)
}
# duplicate all sheets
plate8_stacked <- duplicate_plates(plate1_stacked)
plate9_stacked <- duplicate_plates(plate2_stacked)
plate10_stacked <- duplicate_plates(plate3_stacked)
plate11_stacked <- duplicate_plates(plate4_stacked)
plate12_stacked <- duplicate_plates(plate5_stacked)
plate13_stacked <- duplicate_plates(plate6_stacked)
plate14_stacked <- duplicate_plates(plate7_stacked)

# assign plate numbers to community matrix data frame
for (i in 1:nrow(com_mat)) {
  if (com_mat$sample_name[i] %in% plate1_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 1
    } else {
      com_mat$plate_number[i] <- 8
    }
  }
  if (com_mat$sample_name[i] %in% plate2_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 2
    } else {
      com_mat$plate_number[i] <- 9
    }
  }
  if (com_mat$sample_name[i] %in% plate3_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 3
    } else {
      com_mat$plate_number[i] <- 10
    }
  }
  if (com_mat$sample_name[i] %in% plate4_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 4
    } else {
      com_mat$plate_number[i] <- 11
    }
  }
  if (com_mat$sample_name[i] %in% plate5_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 5
    } else {
      com_mat$plate_number[i] <- 12
    }
  }
  if (com_mat$sample_name[i] %in% plate6_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 6
    } else {
      com_mat$plate_number[i] <- 13
    }
  }
  if (com_mat$sample_name[i] %in% plate7_stacked$sample_id) {
    if (com_mat$replicate[i] == "replicate 1") {
      com_mat$plate_number[i] <- 7
    } else {
      com_mat$plate_number[i] <- 14
    }
  }
}


#-----------------------------------------------
# Bray-Curtis Dissimilarities Between Replicates

bray_curtis_dist <- vegdist(com_mat[,c(5:ncol(com_mat))], method = "bray")
bray_curtis_matrix <- as.matrix(bray_curtis_dist)[com_mat$plate_number == 1, com_mat$plate_number == 8]

# initialize dissimilarity data frame
dissimilarity_df <- data.frame(sample_name = character(), 
                               dissimilarity = numeric(), 
                               plate_number = numeric())

# loop through each sample, calculating dissimilarity between its replicate
for (i in unique(com_mat$sample_name)) {
  # subset data
  subset <- com_mat[com_mat$sample_name == i, ]
  # ensure subset has exactly 2 rows (2 replicates)
  if (nrow(subset) == 2) {
    # calculate Bray-Curtis dissimilarity
    bray_curtis_dist <- vegdist(subset[, c(5:ncol(subset))]/rowSums(subset[, c(5:ncol(subset))]), method = "bray")
    # extract dissimilarity value
    dissimilarity_value <- as.numeric(bray_curtis_dist)
    # store result in data frame
    dissimilarity_df <- rbind(dissimilarity_df, data.frame(
      sample_name = i,
      plate_number = subset$plate_number[1],
      dissimilarity = dissimilarity_value
    ))
  } else {
    # print message if the sample does not have exactly 2 replicates
    message(paste("Skipping sample", i, "because it does not have 2 replicates"))
  }
}

View(dissimilarity_df)

# can now calculate average dissimilarity between replicates for each plate

dissimilarity_comparisons <- dissimilarity_df %>%
  group_by(plate_number) %>%
  summarize(average_dissimilarity = mean(dissimilarity, na.rm = TRUE),
            sd_dissimilarity = sd(dissimilarity, na.rm = TRUE))

boxplot(dissimilarity_df$dissimilarity~dissimilarity_df$plate_number)



