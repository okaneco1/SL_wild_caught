# Multidimensional Distance

# Load libraries and data
library(tidyverse)
library(openxlsx)
library(vegan)
library(patchwork)

# load data
full_data <- read.csv("full_community_matrix_wild_caught.csv")
full_detection_data <- read.csv("full_detection_data_sub2.csv")

#---------- prepping data for distance analyses
full_data <- full_data %>%
  mutate(phase = ifelse(sample_id == "H22_120_2.12S", "Parasitic", phase))

# check for empty rows and remove them
empty_rows <- full_data %>%
  select(Coregoninae_unclassified:Micropterus_unclassified) %>%
  rowSums() == 0
full_data <- full_data[!empty_rows, ]

# exclude metadata to just get sequences
sequence_read_data <- full_data %>%
  select(Coregoninae_unclassified:Micropterus_unclassified) %>%
  select(-Petromyzontidae_unclassified)



#----------------------------------------------------------------------
# NMDS Comparisons (Function)
#----------------------------------------------------------------------

#---------- function to subset, calculate distance and nmds scores, and get centroids
nmds_calculation <- function(subset_phase=NULL, subset_year=NULL, outlier=NULL){
  # copy full data to work with
  filtered_data <- full_data
  
  # remove outliers if specified
  if (!is.null(outlier)) {
    filtered_data <- filtered_data %>%
      filter(!sample_name %in% outlier)
  }
  
  # filter by phase if specified
  if (!is.null(subset_phase)) {
    filtered_data <- filtered_data %>%
      filter(phase == subset_phase)
  }
  
  # filter by year if specified
  if (!is.null(subset_year)) {
    filtered_data <- filtered_data %>%
      filter(year == subset_year)
  }
  
  # subset data
  subset <- filtered_data %>%
    column_to_rownames(var = "sample_id") %>%
    select(Coregoninae_unclassified:Micropterus_unclassified)

  # calculate BC distance
  bray_curtis_matrix <- vegdist(subset, method = "bray")
  # perform NMDS
  nmds_result <- metaMDS(bray_curtis_matrix,
                         distance = "bray",
                         maxit = 999,
                         k = 2, 
                         trymax = 50)
  
  # merge data
  nmds_scores <- as.data.frame(scores(nmds_result))
  nmds_df <- nmds_scores %>%
    rownames_to_column(var = "sample_id") %>%
    inner_join(., full_data, by = "sample_id")
  
  #output
  return(nmds_df)
}

# a rule of thumbs is that a stress value of 1.0-2.0 (as seen above) can 
# be considered "fair" and that some distances can be misleading for interpretation

# as some NMDS analyses had stress values > 2.0, these are considered less
# interpretable, and PCoA will be used for multidimentional analyses




#------------------------------------------------------------------------
# PCoA Analyses (Total Reads)
#------------------------------------------------------------------------

# function for pcoa analyses
pcoa_calculation <- function(subset_phase=NULL, subset_year=NULL, outlier=NULL){
  # copy full data to work with
  filtered_data <- full_data
  
  # remove outliers if specified
  if (!is.null(outlier)) {
    filtered_data <- filtered_data %>%
      filter(!sample_name %in% outlier)
  }
  
  # filter by phase if specified
  if (!is.null(subset_phase)) {
    filtered_data <- filtered_data %>%
      filter(phase == subset_phase)
  }
  
  # filter by year if specified
  if (!is.null(subset_year)) {
    filtered_data <- filtered_data %>%
      filter(year == subset_year)
  }
  
  # subset data
  subset <- filtered_data %>%
    column_to_rownames(var = "sample_id") %>%
    select(Coregoninae_unclassified:Micropterus_unclassified)
  
  # calculate BC distance
  bray_curtis_matrix <- vegdist(subset, method = "bray")
  
  #perform PCoA
  pcoa <- cmdscale(bray_curtis_matrix, k = 2)
  
  # extract PCoA scores
  colnames(pcoa) <- c("PCoA1", "PCoA2")
  
  #merge data
  pcoa_df <- pcoa %>% 
    as.data.frame() %>%
    rownames_to_column(var = "sample_id") %>%
    inner_join(full_data, by = "sample_id")
  
  #output
  return(pcoa_df)
}

# set up PCoA analyses
adult_22_pcoa <- pcoa_calculation(subset_phase="Adult", subset_year=2022)
adult_23_pcoa <- pcoa_calculation(subset_phase="Adult", subset_year=2023)
parasitic_22_pcoa <- pcoa_calculation(subset_phase="Parasitic", subset_year=2022)
parasitic_23_pcoa <- pcoa_calculation(subset_phase="Parasitic", subset_year=2023)

# year-specific
pcoa_2022 <- pcoa_calculation(subset_year=2022)
pcoa_2023 <- pcoa_calculation(subset_year=2023)


# Visualizations
adult_22_pcoa_plot <- ggplot(adult_22_pcoa, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Adult 2022 PCoA (n = ",nrow(adult_22_pcoa),")"),
       x = "PCoA1", y = "PCoA2")

adult_23_pcoa_plot <- ggplot(adult_23_pcoa, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Adult 2023 PCoA (n = ",nrow(adult_23_pcoa),")"),
       x = "PCoA1", y = "PCoA2")

parasitic_22_pcoa_plot <- ggplot(parasitic_22_pcoa, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Parasitic 2022 PCoA (n = ",nrow(parasitic_22_pcoa),")"),
       x = "PCoA1", y = "PCoA2")

parasitic_23_pcoa_plot <- ggplot(parasitic_23_pcoa, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Parasitic 2023 PCoA (n = ",nrow(parasitic_23_pcoa),")"),
       x = "PCoA1", y = "PCoA2")

# plot together
adult_22_pcoa_plot + adult_23_pcoa_plot + parasitic_22_pcoa_plot + parasitic_23_pcoa_plot


# year-specific
pcoa_2022_plot <- ggplot(pcoa_2022, aes(x = PCoA1, y = PCoA2, color = lake, shape=phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("All 2022 PCoA (n = ",nrow(pcoa_2022),")"),
       x = "PCoA1", y = "PCoA2")

pcoa_2023_plot <- ggplot(pcoa_2023, aes(x = PCoA1, y = PCoA2, color = lake, shape=phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("All 2023 PCoA (n = ",nrow(pcoa_2023),")"),
       x = "PCoA1", y = "PCoA2")

pcoa_2022_plot + pcoa_2023_plot



#--------------- Including Multiple Ellispses
pcoa_2023_plot <- ggplot(pcoa_2023, aes(x = PCoA1, y = PCoA2, color = lake, shape = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  theme_minimal() +
  labs(title = paste0("All 2023 PCoA (n = ", nrow(pcoa_2023), ")"),
       x = "PCoA1", y = "PCoA2")

# Add ellipses for each phase with different linetypes
pcoa_2023_plot <- pcoa_2023_plot +
  stat_ellipse(data = filter(pcoa_2023, phase == "Parasitic"), 
               aes(x = PCoA1, y = PCoA2, color = lake), 
               linetype = "dashed", 
               type = "t") +
  stat_ellipse(data = filter(pcoa_2023, phase == "Adult"), 
               aes(x = PCoA1, y = PCoA2, color = lake), 
               linetype = "solid", 
               type = "t")

# Print the plot
pcoa_2023_plot






#------------------------------------------------------------------
# Comparing Relative Abundance (PCoA)
#------------------------------------------------------------------

full_data_relative_reads <- full_data %>%
  select(Coregoninae_unclassified:Micropterus_unclassified)

full_data_relative_reads <- full_data_relative_reads/rowSums(full_data_relative_reads)
full_data_relative <- full_data %>%
  select(sample_id:phase) %>%
  cbind(., full_data_relative_reads)
# sum(rowSums(full_data_relative))==nrow(full_data_relative)

# NMDS Calculations with relative data
# Function
pcoa_relative_calculation <- function(subset_phase=NULL, 
                                      subset_year=NULL, 
                                      subset_lake=NULL,
                                      outlier=NULL){
  # copy full data to work with
  filtered_data <- full_data_relative
  
  # remove outliers if specified
  if (!is.null(outlier)) {
    filtered_data <- filtered_data %>%
      filter(!sample_name %in% outlier)
  }
  
  # filter by lake if specified
  if (!is.null(subset_lake)) {
    filtered_data <- filtered_data %>%
      filter(lake == subset_lake)
  }
  
  # filter by phase if specified
  if (!is.null(subset_phase)) {
    filtered_data <- filtered_data %>%
      filter(phase == subset_phase)
  }
  
  # filter by year if specified
  if (!is.null(subset_year)) {
    filtered_data <- filtered_data %>%
      filter(year == subset_year)
  }
  
  # calculate BC distance
  sequences <- filtered_data %>%
    column_to_rownames(var = "sample_id") %>%
    select(Coregoninae_unclassified:Micropterus_unclassified)
  sequences[is.na(sequences)] <- 0 # replace any missing values with 0
  sequences <- sequences %>%
    filter(rowSums(.) != 0)
  bray_curtis_matrix <- vegdist(sequences, method = "bray")
  
  #perform PCoA
  pcoa <- cmdscale(bray_curtis_matrix, k = 2)
  
  # extract PCoA scores
  colnames(pcoa) <- c("PCoA1", "PCoA2")
  
  #merge data
  pcoa_df <- pcoa %>% 
    as.data.frame() %>%
    rownames_to_column(var = "sample_id") %>%
    inner_join(full_data, by = "sample_id")
  
  #output
  return(pcoa_df)
}

#------------ LAKE DIFFERENCES
# run nmds calculations for each phase/year combination
adult_22_pcoa_rel <- pcoa_relative_calculation(subset_phase="Adult", subset_year=2022)
adult_23_pcoa_rel <- pcoa_relative_calculation(subset_phase="Adult", subset_year=2023,
                                               outlier = "S23_A11")
parasitic_22_pcoa_rel <- pcoa_relative_calculation(subset_phase="Parasitic", subset_year=2022)
parasitic_23_pcoa_rel <- pcoa_relative_calculation(subset_phase="Parasitic", subset_year=2023, outlier="S23_P4")


# VISUALIZATIONS
# maintain a common axis range
common_xlim <- c(-1, 1)
common_ylim <- c(-1, 1)

adult_22_pcoa_rel_plot <- ggplot(adult_22_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Adult 2022 PCoA (n = ",nrow(adult_22_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)

adult_23_pcoa_rel_plot <- ggplot(adult_23_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Adult 2023 PCoA (n = ",nrow(adult_23_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)

parasitic_22_pcoa_rel_plot <- ggplot(parasitic_22_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Parasitic 2022 PCoA (n = ",nrow(parasitic_22_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)

parasitic_23_pcoa_rel_plot <- ggplot(parasitic_23_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = lake)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Parasitic 2023 PCoA (n = ",nrow(parasitic_23_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)

# plot together
adult_22_pcoa_rel_plot + adult_23_pcoa_rel_plot + parasitic_22_pcoa_rel_plot + parasitic_23_pcoa_rel_plot



#------------ PHASE DIFFERENCES
huron_22_pcoa_rel <- pcoa_relative_calculation(subset_lake="Huron", subset_year=2022)
huron_23_pcoa_rel <- pcoa_relative_calculation(subset_lake="Huron", subset_year=2023)
champlain_22_pcoa_rel <- pcoa_relative_calculation(subset_lake="Champlain", subset_year=2022)
champlain_23_pcoa_rel <- pcoa_relative_calculation(subset_lake="Champlain", subset_year=2023)
superior_22_pcoa_rel <- pcoa_relative_calculation(subset_lake="Superior", subset_year=2022)
superior_23_pcoa_rel <- pcoa_relative_calculation(subset_lake="Superior", subset_year=2023)

# VISUALIZATIONS

# Huron 2022 (stress = )
h22_plot_rel <- ggplot(huron_22_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Huron 2022 (n = ",nrow(huron_22_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Huron 2023 (stress = )
h23_plot_rel <- ggplot(huron_23_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Huron 2023 (n = ",nrow(huron_23_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Champlain 2022 (stress = )
c22_plot_rel <- ggplot(champlain_22_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Champlain 2022 (n = ",nrow(champlain_22_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Champlain 2023 (stress = )
c23_plot_rel <- ggplot(champlain_23_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Champlain 2023 (n = ",nrow(champlain_23_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Superior 2022 (stress = )
s22_plot_rel <- ggplot(superior_22_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Superior 2022 (n = ",nrow(superior_22_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)

# Superior 2023 (stress = )
s23_plot_rel <- ggplot(superior_23_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = phase)) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Superior 2023 (n = ",nrow(superior_23_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)



h22_plot_rel + h23_plot_rel + c22_plot_rel + c23_plot_rel + s22_plot_rel + s23_plot_rel



#------------ YEAR DIFFERENCES
huron_adult_pcoa_rel <- pcoa_relative_calculation(subset_lake="Huron", subset_phase="Adult")
huron_parasitic_pcoa_rel <- pcoa_relative_calculation(subset_lake="Huron", subset_phase="Parasitic")
champlain_adult_pcoa_rel <- pcoa_relative_calculation(subset_lake="Champlain", subset_phase="Adult")
champlain_parasitic_pcoa_rel <- pcoa_relative_calculation(subset_lake="Champlain", subset_phase="Parasitic")
superior_adult_pcoa_rel <- pcoa_relative_calculation(subset_lake="Superior", subset_phase="Adult")
superior_parasitic_pcoa_rel <- pcoa_relative_calculation(subset_lake="Superior", subset_phase="Parasitic")

# VISUALIZATIONS

# Huron 2022 (stress = )
h_adult_plot_rel <- ggplot(huron_adult_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = as.factor(year))) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Huron Adult (n = ",nrow(huron_adult_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Huron 2023 (stress = )
h_parasitic_plot_rel <- ggplot(huron_parasitic_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = as.factor(year))) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Huron Parasitic (n = ",nrow(huron_parasitic_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Champlain 2022 (stress = )
c_adult_plot_rel <- ggplot(champlain_adult_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = as.factor(year))) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Champlain Adult (n = ",nrow(champlain_adult_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Champlain 2023 (stress = )
c_parasitic_plot_rel <- ggplot(champlain_parasitic_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = as.factor(year))) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Champlain Parasitic (n = ",nrow(champlain_parasitic_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2") +
  xlim(common_xlim) +
  ylim(common_ylim)

# Superior 2022 (stress = )
s_adult_plot_rel <- ggplot(superior_adult_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = as.factor(year))) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Superior Adult (n = ",nrow(superior_adult_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)

# Superior 2023 (stress = )
s_parasitic_plot_rel <- ggplot(superior_parasitic_pcoa_rel, aes(x = PCoA1, y = PCoA2, color = as.factor(year))) +
  geom_point(size = 1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Superior Parasitic (n = ",nrow(superior_parasitic_pcoa_rel),")"),
       x = "PCoA1", y = "PCoA2")+
  xlim(common_xlim) +
  ylim(common_ylim)



h_adult_plot_rel + h_parasitic_plot_rel + c_adult_plot_rel + c_parasitic_plot_rel + s_adult_plot_rel + s_parasitic_plot_rel









#------------------------------------------------------------------
# DETECTION DATA (PCoA)
#------------------------------------------------------------------
empty_rows <- full_detection_data %>%
  select(Coregoninae_unclassified:Micropterus_unclassified) %>%
  rowSums() == 0
full_detection_data <- full_detection_data[!empty_rows, ] # removed 101 rows

# exclude metadata to just get sequences
sequence_read_data <- full_detection_data %>%
  select(Coregoninae_unclassified:Micropterus_unclassified) %>%
  select(-Petromyzontidae_unclassified)

# nmds function
nmds_calculation_det <- function(subset_phase=NULL, subset_year=NULL, outlier=NULL){
  # copy full data to work with
  filtered_data <- full_detection_data
  
  # remove outliers if specified
  if (!is.null(outlier)) {
    filtered_data <- filtered_data %>%
      filter(!sample_name %in% outlier)
  }
  
  # filter by phase if specified
  if (!is.null(subset_phase)) {
    filtered_data <- filtered_data %>%
      filter(phase == subset_phase)
  }
  
  # filter by year if specified
  if (!is.null(subset_year)) {
    filtered_data <- filtered_data %>%
      filter(year == subset_year)
  }
  
  # subset data
  subset <- filtered_data %>%
    column_to_rownames(var = "sample_id") %>%
    select(Coregoninae_unclassified:Micropterus_unclassified)
  
  # calculate BC distance
  bray_curtis_matrix <- vegdist(subset, method = "bray")
  # perform NMDS
  nmds_result <- metaMDS(bray_curtis_matrix,
                         distance = "bray",
                         maxit = 999,
                         k = 2, 
                         trymax = 30)
  
  # merge data
  nmds_scores <- as.data.frame(scores(nmds_result))
  nmds_df <- nmds_scores %>%
    rownames_to_column(var = "sample_id") %>%
    inner_join(., full_detection_data, by = "sample_id")
  
  #output
  return(nmds_df)
}

# describe comparisons
adult_22_nmds_det <- nmds_calculation_det(subset_phase="Adult", subset_year=2022, outlier="C22_A83")
adult_23_nmds_det <- nmds_calculation_det(subset_phase="Adult", subset_year=2023)
parasitic_22_nmds_det <- nmds_calculation_det(subset_phase="Parasitic", subset_year=2022,
                                              outlier = c("H22_P165","H22_P14","H22_P54"))
parasitic_23_nmds_det <- nmds_calculation_det(subset_phase="Parasitic", subset_year=2023)

# VISUALIZATIONS
# Adult 2022 (stress =)
a22_plot_det <- ggplot(adult_22_nmds_det, aes(x = NMDS1, y = NMDS2, color = lake)) +
  geom_jitter(size = 1, width = 0.1, height = 0.1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Adult 2022 Detections (n = ",nrow(adult_22_nmds_det),")"),
       x = "NMDS1", y = "NMDS2")

# Adult 2023 (stress = 0.1852066)
a23_plot_det <- ggplot(adult_23_nmds_det, aes(x = NMDS1, y = NMDS2, color = lake)) +
  geom_jitter(size = 1, width = 0.1, height = 0.1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Adult 2023 Detections (n = ",nrow(adult_23_nmds_det),")"),
       x = "NMDS1", y = "NMDS2")

# Parasitic 2022 (stress = 0.1510647)
p22_plot_det <- ggplot(parasitic_22_nmds_det, aes(x = NMDS1, y = NMDS2, color = lake)) +
  geom_jitter(size = 1, width = 0.1, height = 0.1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 3) +
  theme_minimal() +
  labs(title = paste0("Parasitic 2022 Detections (n = ",nrow(parasitic_22_nmds_det),")"),
       x = "NMDS1", y = "NMDS2")

# Parasitic 2023 (stress = 0.1593906)
p23_plot_det <- ggplot(parasitic_23_nmds_det, aes(x = NMDS1, y = NMDS2, color = lake)) +
  geom_jitter(size = 1, width = 0.1, height = 0.1, alpha = 0.5) +
  stat_ellipse(type = "t") +  
  #geom_text(aes(label = sample_id), vjust = -1, size = 1) +
  theme_minimal() +
  labs(title = paste0("Parasitic 2023 Detections (n = ",nrow(parasitic_23_nmds_det),")"),
       x = "NMDS1", y = "NMDS2")

# print all together
a22_plot_det + a23_plot_det + p22_plot_det + p23_plot_det






