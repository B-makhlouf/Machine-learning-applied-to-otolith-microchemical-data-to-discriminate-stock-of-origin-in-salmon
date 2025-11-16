# CWAK Analysis - Random Forest Model Performance
# Analyzes RF GAM model performance for CWAK vs non-CWAK fish populations
#
# Prerequisites:
#   1. NatalToMarine_Processed_GAM.csv (preprocessed GAM data)
#   2. Fish_ID_Splits.csv (train/test split assignments)
#   3. GAM_RF_TOTAL_predictions.csv (RF model predictions)
#   4. All_Yukon_Genetics.csv (genetic assignments for Yukon fish)

library(tidyverse)
library(here)

# =============================================================================
# LOAD DATA
# =============================================================================

# Load GAM test data
gam_data <- read.csv(here("data", "LA_Data", "NatalToMarine_Processed_GAM.csv"))

# Load test fish IDs
fish_splits <- read.csv(here("data", "TrainingTesting", "Fish_ID_Splits.csv"))
test_fish_ids <- fish_splits$Fish_id[fish_splits$Split == "Test"]

# Load RF predictions
rf_predictions <- read.csv(here("Output", "ModelResultsPreCal", "Total", "GAM_RF_TOTAL_predictions.csv"))

# Load genetic data
all_genetic_data <- read.csv(here("data", "All_Yukon_Genetics.csv"))

# =============================================================================
# PREPARE TEST DATASET
# =============================================================================

# Filter GAM data to test fish only
gam_test <- gam_data %>%
  filter(Fish_id %in% test_fish_ids)

# Merge RF predictions with metadata
if ("Watershed" %in% colnames(rf_predictions)) {
  rf_with_metadata <- rf_predictions %>%
    mutate(
      Fish_id = gam_test$Fish_id,
      predicted_watershed = .pred_class,
      correct_prediction = (Watershed == .pred_class)
    )
} else {
  rf_with_metadata <- rf_predictions %>%
    mutate(
      Fish_id = gam_test$Fish_id,
      Watershed = gam_test$Watershed,
      predicted_watershed = .pred_class,
      correct_prediction = (gam_test$Watershed == .pred_class)
    )
}

# =============================================================================
# MERGE WITH GENETIC DATA
# =============================================================================

# Remove duplicate Watershed column if present in genetic data
if ("Watershed" %in% colnames(all_genetic_data)) {
  genetic_data_clean <- all_genetic_data %>%
    select(-Watershed)
} else {
  genetic_data_clean <- all_genetic_data
}

# Merge with genetic assignments
rf_with_genetics <- rf_with_metadata %>%
  left_join(genetic_data_clean, by = "Fish_id")

# =============================================================================
# CREATE CWAK GROUPINGS
# =============================================================================

# CWAK Definition:
# - All fish from Nushagak (Nush) and Kuskokwim (Kusko) watersheds
# - Yukon fish genetically assigned to Lower Yukon
# 
# non-CWAK:
# - Yukon fish genetically assigned to Middle or Upper Yukon

rf_with_cwak <- rf_with_genetics %>%
  mutate(
    cwak_group = case_when(
      Watershed %in% c("Nush", "Kusko") ~ "CWAK",
      Watershed == "Yukon" & genetic_group == "Lower" ~ "CWAK", 
      Watershed == "Yukon" & genetic_group %in% c("Middle", "Upper") ~ "non-CWAK",
      TRUE ~ NA_character_
    )
  )

# Filter to fish with clear CWAK classification
rf_final <- rf_with_cwak %>%
  filter(!is.na(cwak_group))

# =============================================================================
# ANALYZE PERFORMANCE
# =============================================================================

# Overall performance by CWAK group
cwak_performance <- rf_final %>%
  group_by(cwak_group) %>%
  summarise(
    n_fish = n(),
    correct_predictions = sum(correct_prediction),
    accuracy = correct_predictions / n_fish,
    .groups = "drop"
  )

cat("=== RF GAM MODEL: CWAK vs NON-CWAK PERFORMANCE ===\n")
print(cwak_performance)

# Detailed breakdown by watershed
detailed_performance <- rf_final %>%
  group_by(cwak_group, Watershed) %>%
  summarise(
    n_fish = n(),
    correct_predictions = sum(correct_prediction),
    accuracy = correct_predictions / n_fish,
    .groups = "drop"
  )

cat("\n=== DETAILED BREAKDOWN BY WATERSHED ===\n")
print(detailed_performance)

# Yukon fish breakdown by genetic group
yukon_genetic_breakdown <- rf_final %>%
  filter(Watershed == "Yukon") %>%
  group_by(genetic_group, cwak_group) %>%
  summarise(
    n_fish = n(),
    correct_predictions = sum(correct_prediction),
    accuracy = correct_predictions / n_fish,
    .groups = "drop"
  )

cat("\n=== YUKON FISH BY GENETIC GROUP ===\n")
print(yukon_genetic_breakdown)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Create output directory
output_dir <- here("Output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Save performance summaries
write.csv(cwak_performance, 
          here(output_dir, "CWAK_RF_Performance_Summary.csv"), 
          row.names = FALSE)

write.csv(detailed_performance, 
          here(output_dir, "CWAK_RF_Detailed_Performance.csv"), 
          row.names = FALSE)

write.csv(yukon_genetic_breakdown,
          here(output_dir, "CWAK_RF_Yukon_Genetic_Breakdown.csv"),
          row.names = FALSE)

write.csv(rf_final,
          here(output_dir, "RF_Predictions_with_CWAK_Groups.csv"), 
          row.names = FALSE)

cat("\n✓ Analysis complete. Results saved to:", output_dir, "\n")