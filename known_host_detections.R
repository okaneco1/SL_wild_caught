# Known Host Detections

# Load libraries and data
library(tidyverse)
library(openxlsx)

# load data
detection_data <- read_csv("full_detection_data_sub2.csv")

detection_data <- detection_data %>%
  mutate(total_detections = rowSums(select(., 9:(ncol(.) - 1)), na.rm = TRUE))


#------------------------------------------------------------------------
counts <- detection_data %>%
  count(known_host)

# combining "trout" and lake trout
detection_data$known_host[detection_data$known_host == "trout"] <- "lake_trout"

# combining chinook and king salmon
detection_data$known_host[detection_data$known_host == "king_salmon"] <- "chinook_salmon"



#---------------------------------------------------------------
# Build a table with the data
known_host_data <- detection_data %>%
  count(known_host) %>%
  rename(samples = n) %>%
  arrange(desc(samples)) %>%
  filter(!is.na(known_host)) %>%
  mutate(known_detected = NA,
         other_detected = NA,
         zero_detections = NA,
         top5_species = NA)

# set up some mapping 
host_variable_map <- list(
  whitefish = "Coregoninae_unclassified",
  lake_trout = "Salvelinus_namaycush",
  pink_salmon = "Oncorhynchus_nerka",
  rainbow_trout = "Oncorhynchus_mykiss",
  burbot = "Lota_lota",
  herring = "Clupeidae_unclassified",
  pickerel = "Esox_unclassified",
  salmon = "Salmonidae_unclassified",
  atlantic_salmon = "Salmo_salar",
  chinook_salmon = "Oncorhynchus_tshawytscha",
  sucker = "Catostomus_catostomus",
  walleye = "Sander_vitreus"
)

# How many times was the known host detected
for (i in 1:nrow(known_host_data)) {
  host <- known_host_data$known_host[i]  
  
  # get the corresponding column from the map
  variable <- host_variable_map[[host]]
  
  # check if the column exists in detection_data
  if (!is.null(variable) && variable %in% colnames(detection_data)) {
    # filter rows for the current host
    species <- detection_data %>% 
      filter(known_host == host)
    
    # sum relevant column and update known_detected
    known_host_data$known_detected[i] <- sum(species[[variable]], na.rm = TRUE)
  }
}


# How many had detections from other species-specific OTUs
for (i in 1:nrow(known_host_data)) {
  host <- known_host_data$known_host[i] 
  
  primary_species <- host_variable_map[[host]]
  
  if (!is.null(primary_species) && primary_species %in% colnames(detection_data)) {
    species_data <- detection_data %>%
      filter(known_host == host)
    
    # include Coregoninae_unclassified and other columns without "_unclassified" in their name
    other_species_count <- species_data %>%
      select(Coregoninae_unclassified:Micropterus_unclassified) %>%  
      select(Coregoninae_unclassified, matches("^(?!.*_unclassified).*$", perl = TRUE)) %>%  
      { if (primary_species %in% colnames(.)) select(., -primary_species) else . } %>%  
      mutate(total_detections = rowSums(.)) %>%  
      nrow()
    
    known_host_data$other_detected[i] <- other_species_count
  }
}


# What were the top 5 OTUs for each sample with known hosts
for (i in 1:nrow(known_host_data)) {
  host <- known_host_data$known_host[i]  
  
  # filter detection_data 
  species_data <- detection_data %>%
    filter(known_host == host)
  
  # summing detections for all OTUs 
  otu_sums <- species_data %>%
    select(Coregoninae_unclassified:Micropterus_unclassified) %>%  
    summarise(across(everything(), sum, na.rm = TRUE))  
  
  # convert to a tidy format for ranking
  otu_sums_tidy <- otu_sums %>%
    pivot_longer(cols = everything(), names_to = "OTU", values_to = "detection_sum") %>%
    filter(detection_sum > 0) %>%  
    arrange(desc(detection_sum))  
  
  # get top 5 OTUs by detections, excluding "Salmonidae_unclassified"
  top5 <- otu_sums_tidy %>%
    filter(OTU != "Salmonidae_unclassified") %>%  
    slice_head(n = 5)  
  
  # create the top 5 OTUs with counts in parentheses
  top5_with_counts <- top5 %>%
    mutate(label = paste(OTU, " (", detection_sum, ")", sep = "")) %>%  
    pull(label)  
  
  known_host_data$top5_species[i] <- paste(top5_with_counts, collapse = ", ")
}


# write out the table
write_csv(known_host_data, file = "known_host_data.csv")





#----------------------------------------------------------
# also just checking for total non-detections
adult_non_detections <- detection_data %>%
  filter(phase == "Adult") %>%
  filter(total_detections == 0) %>%
  nrow()

parasitic_non_detections <- detection_data %>%
  filter(phase == "Parasitic") %>%
  filter(total_detections == 0) %>%
  nrow()













