# Streamlined Modeling Script - STARTS AFTER SPLITTING
# This script loads pre-existing train/test splits and runs the full analysis
library(here)
library(dplyr)
library(tidymodels)
library(ranger)
library(kernlab)
library(kknn)
library(ggplot2)
library(viridis)
library(scales)
library(forcats)
library(tidyr)

################################################################################
#### CONFIGURATION
################################################################################

NATAL_ISO_THRESHOLD <- 0.713
set.seed(123)

data_types <- c("RAW", "GAM", "MA")

# Input directories (where train/test splits already exist)
train_test_dir_total <- "/Users/benjaminmakhlouf/Research_repos/Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Data/TrainingTesting"
train_test_dir_overlap <- "/Users/benjaminmakhlouf/Research_repos/Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Data/TrainingTesting/Filtered"

# Output directories
models_dir_total <- "/Users/benjaminmakhlouf/Research_repos/Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Output/Models/Total"
models_dir_overlap <- "/Users/benjaminmakhlouf/Research_repos/Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Output/Models/Filtered"
results_dir_total <- "/Users/benjaminmakhlouf/Research_repos/Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Output/ModelResultsPreCal/Total"
results_dir_overlap <- "/Users/benjaminmakhlouf/Research_repos/0Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Output/ModelResultsPreCal/Filtered"
figures_dir <- "/Users/benjaminmakhlouf/Research_repos/Machine-learning-applied-to-otolith-microchemical-data-to-discriminate-stock-of-origin-in-salmon/Figures/ModelPerformance"

for(dir in c(models_dir_total, models_dir_overlap, results_dir_total, results_dir_overlap, figures_dir)) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}

################################################################################
#### STEP 0: Verify required files exist
################################################################################

cat("\n=== Checking for required train/test split files ===\n")

# Check Total analysis files
total_files_exist <- TRUE
for (data_type in data_types) {
  train_file <- file.path(train_test_dir_total, paste0("Train_", data_type, ".csv"))
  test_file <- file.path(train_test_dir_total, paste0("Test_", data_type, ".csv"))
  
  if (!file.exists(train_file)) {
    cat("✗ MISSING:", train_file, "\n")
    total_files_exist <- FALSE
  } else {
    cat("✓ Found:", basename(train_file), "\n")
  }
  
  if (!file.exists(test_file)) {
    cat("✗ MISSING:", test_file, "\n")
    total_files_exist <- FALSE
  } else {
    cat("✓ Found:", basename(test_file), "\n")
  }
}

# Check Overlap/Filtered analysis files
overlap_files_exist <- TRUE
for (data_type in data_types) {
  train_file <- file.path(train_test_dir_overlap, paste0("Train_", data_type, ".csv"))
  test_file <- file.path(train_test_dir_overlap, paste0("Test_", data_type, ".csv"))
  
  if (!file.exists(train_file)) {
    cat("✗ MISSING:", train_file, "\n")
    overlap_files_exist <- FALSE
  } else {
    cat("✓ Found (Filtered):", basename(train_file), "\n")
  }
  
  if (!file.exists(test_file)) {
    cat("✗ MISSING:", test_file, "\n")
    overlap_files_exist <- FALSE
  } else {
    cat("✓ Found (Filtered):", basename(test_file), "\n")
  }
}

# Check for Fish_ID_Splits.csv (optional but helpful)
fish_id_file <- file.path(train_test_dir_total, "Fish_ID_Splits.csv")
if (file.exists(fish_id_file)) {
  cat("✓ Found: Fish_ID_Splits.csv\n")
  fish_splits <- read.csv(fish_id_file)
  cat("  Train Fish IDs:", sum(fish_splits$Split == "Train"), "\n")
  cat("  Test Fish IDs:", sum(fish_splits$Split == "Test"), "\n")
} else {
  cat("⚠ Fish_ID_Splits.csv not found (optional)\n")
}

if (!total_files_exist || !overlap_files_exist) {
  stop("\n❌ ERROR: Missing required train/test split files. Please run the splitting script first.\n")
}

cat("\n✓ All required files found! Proceeding with analysis...\n")

################################################################################
#### STEP 1: Run analysis
################################################################################

run_analysis <- function(train_test_dir, models_dir, analysis_name, results_dir) {
  results <- data.frame()
  
  for (data_type in data_types) {
    train_file <- file.path(train_test_dir, paste0("Train_", data_type, ".csv"))
    test_file <- file.path(train_test_dir, paste0("Test_", data_type, ".csv"))
    
    if (!file.exists(train_file) || !file.exists(test_file)) next
    
    train_data <- read.csv(train_file) %>% mutate(Watershed = as.factor(Watershed))
    test_data <- read.csv(test_file) %>% mutate(Watershed = as.factor(Watershed))
    if (nrow(test_data) == 0) next
    
    cat("Processing", data_type, "- Train:", nrow(train_data), "Test:", nrow(test_data), "\n")
    
    base_recipe <- recipe(Watershed ~ ., data = train_data)
    n_predictors <- ncol(train_data) - 1
    
    models <- list(
      RF = rand_forest(trees = 500, mtry = floor(sqrt(n_predictors))) %>% set_engine("ranger") %>% set_mode("classification"),
      SVM = svm_rbf() %>% set_engine("kernlab") %>% set_mode("classification"),
      KNN = nearest_neighbor(neighbors = 5) %>% set_engine("kknn") %>% set_mode("classification")
    )
    
    for (model_name in names(models)) {
      set.seed(123)
      
      workflow_obj <- workflow() %>% add_recipe(base_recipe) %>% add_model(models[[model_name]]) %>% fit(train_data)
      saveRDS(workflow_obj, file.path(models_dir, paste0(data_type, "_", model_name, "_model.rds")))
      
      predictions <- workflow_obj %>% predict(test_data) %>% bind_cols(test_data %>% select(Watershed))
      pred_probs <- workflow_obj %>% predict(test_data, type = "prob")
      predictions_with_probs <- predictions %>% bind_cols(pred_probs) %>% mutate(Dataset = data_type, Model = model_name, Correct = Watershed == .pred_class)
      
      write.csv(predictions_with_probs, file.path(results_dir, paste0(data_type, "_", model_name, "_", analysis_name, "_predictions.csv")), row.names = FALSE)
      
      accuracy <- mean(predictions$Watershed == predictions$.pred_class)
      f1_score <- predictions %>% f_meas(truth = Watershed, estimate = .pred_class) %>% pull(.estimate)
      
      results <- rbind(results, data.frame(Dataset = data_type, Model = model_name, Accuracy = round(accuracy, 3), F1_Score = round(f1_score, 3)))
    }
  }
  
  results <- results[order(-results$Accuracy), ]
  cat("\n=== Results for", analysis_name, "===\n")
  print(results)
  return(results)
}

cat("\n=== Running TOTAL analysis ===\n")
results_total <- run_analysis(train_test_dir_total, models_dir_total, "TOTAL", results_dir_total)

cat("\n=== Running OVERLAP analysis ===\n")
results_overlap <- run_analysis(train_test_dir_overlap, models_dir_overlap, "OVERLAP", results_dir_overlap)

################################################################################
#### STEP 2: Combined heatmaps
################################################################################

create_combined_heatmaps <- function(results_total, results_overlap) {
  if (nrow(results_total) == 0 || nrow(results_overlap) == 0) return(NULL)
  
  combined_results <- bind_rows(
    results_total %>% mutate(Analysis = "Total"),
    results_overlap %>% mutate(Analysis = "Overlapping")
  ) %>%
    filter(Dataset %in% c("RAW", "GAM", "MA")) %>%
    mutate(
      Dataset_Combined = paste(Dataset, Analysis, sep = " - "),
      Model = factor(Model, levels = c("RF", "SVM", "KNN"), labels = c("Random Forest", "SVM", "KNN")),
      Dataset_Combined = factor(Dataset_Combined, levels = c(paste(c("RAW", "GAM", "MA"), "- Total"), paste(c("RAW", "GAM", "MA"), "- Overlapping")))
    )
  
  # Accuracy heatmap
  accuracy_plot <- ggplot(combined_results, aes(x = Model, y = Dataset_Combined, fill = Accuracy)) +
    geom_tile(color = NA, width = 1, height = 1) +
    geom_text(aes(label = sprintf("%.3f", Accuracy)), color = "black", size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#9AB87A", limits = c(min(combined_results$Accuracy) * 0.99, max(combined_results$Accuracy) * 1.01), labels = scales::percent_format(accuracy = 0.1)) +
    labs(title = "Accuracy", x = "Model Type", y = "Data Source", fill = "Accuracy") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)),
          axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 15)),
          axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 15), angle = 90),
          axis.text.x = element_text(size = 14, margin = margin(t = 8)),
          axis.text.y = element_text(size = 14, margin = margin(r = 8)),
          legend.position = "right", legend.key.height = unit(1.5, "cm"),
          plot.margin = margin(20, 25, 20, 25)) +
    annotate("segment", x = 0.5, xend = 3.5, y = 3.5, yend = 3.5, color = "white", size = 2)
  
  # F1-Score heatmap
  f1_plot <- ggplot(combined_results, aes(x = Model, y = Dataset_Combined, fill = F1_Score)) +
    geom_tile(color = NA, width = 1, height = 1) +
    geom_text(aes(label = sprintf("%.3f", F1_Score)), color = "black", size = 6, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "#9AB87A", limits = c(min(combined_results$F1_Score) * 0.99, max(combined_results$F1_Score) * 1.01), labels = scales::percent_format(accuracy = 0.1)) +
    labs(title = "F1-Score", x = "Model Type", y = "Data Source", fill = "F1-Score") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)),
          axis.title.x = element_text(face = "bold", size = 16, margin = margin(t = 15)),
          axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 15), angle = 90),
          axis.text.x = element_text(size = 14, margin = margin(t = 8)),
          axis.text.y = element_text(size = 14, margin = margin(r = 8)),
          legend.position = "right", legend.key.height = unit(1.5, "cm"),
          plot.margin = margin(20, 25, 20, 25)) +
    annotate("segment", x = 0.5, xend = 3.5, y = 3.5, yend = 3.5, color = "white", size = 2)
  
  ggsave(file.path(figures_dir, "Combined_Accuracy_Heatmap.pdf"), accuracy_plot, width = 10, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(figures_dir, "Combined_F1Score_Heatmap.pdf"), f1_plot, width = 10, height = 6, dpi = 300, bg = "white")
  
  cat("✓ Combined heatmaps saved\n")
  return(list(accuracy_plot = accuracy_plot, f1_plot = f1_plot))
}

combined_plots <- create_combined_heatmaps(results_total, results_overlap)

################################################################################
#### STEP 3: GAM RF specific figures
################################################################################

create_gam_rf_figures <- function(results_dir_total, results_dir_overlap) {
  gam_rf_total_file <- file.path(results_dir_total, "GAM_RF_TOTAL_predictions.csv")
  gam_rf_overlap_file <- file.path(results_dir_overlap, "GAM_RF_OVERLAP_predictions.csv")
  
  if (!file.exists(gam_rf_total_file) || !file.exists(gam_rf_overlap_file)) return(NULL)
  
  combined_gam_rf <- bind_rows(
    read.csv(gam_rf_total_file) %>% mutate(Analysis = "Total"),
    read.csv(gam_rf_overlap_file) %>% mutate(Analysis = "Overlapping")
  )
  
  # Confusion matrices
  create_confusion_data <- function(predictions, analysis_name) {
    # Get all unique watershed levels from the predictions
    all_watersheds <- sort(unique(c(predictions$Watershed, predictions$.pred_class)))
    
    conf_data <- predictions %>%
      count(Watershed, .pred_class, .drop = FALSE) %>%
      group_by(Watershed) %>%
      mutate(
        Total = sum(n),
        Percentage = ifelse(Total > 0, n / Total, 0)
      ) %>%
      ungroup() %>%
      complete(
        Watershed = all_watersheds, 
        .pred_class = all_watersheds, 
        fill = list(n = 0, Total = 0, Percentage = 0)
      ) %>%
      mutate(
        Analysis = analysis_name,
        Label = ifelse(n > 0, as.character(n), "")
      )
    
    return(conf_data)
  }
  
  gam_rf_total <- combined_gam_rf %>% filter(Analysis == "Total")
  gam_rf_overlap <- combined_gam_rf %>% filter(Analysis == "Overlapping")
  
  combined_conf <- bind_rows(
    create_confusion_data(gam_rf_total, "Total"),
    create_confusion_data(gam_rf_overlap, "Overlapping")
  ) %>%
    mutate(Analysis = factor(Analysis, levels = c("Total", "Overlapping")))
  
  confusion_plot <- ggplot(combined_conf, aes(x = .pred_class, y = fct_rev(Watershed))) +
    geom_tile(aes(fill = ifelse(Watershed == .pred_class, "Correct", "Incorrect")), color = "white", size = 1.5) +
    geom_text(aes(label = paste0(sprintf("%.1f%%", Percentage * 100), "\n(", Label, ")")), 
              color = "white", size = 4, fontface = "bold", lineheight = 0.9) +
    scale_fill_manual(
      values = c("Correct" = "#e74c3c", "Incorrect" = "#95a5a6"),
      labels = c("Correct", "Incorrect")
    ) +
    facet_wrap(~Analysis, labeller = labeller(Analysis = c("Total" = "Total Analysis", "Overlapping" = "Restricted Analysis"))) +
    labs(
      title = "Classification Confusion Matrices", 
      subtitle = "GAM Random Forest Model", 
      x = "Predicted Watershed", 
      y = "Actual Watershed", 
      fill = "Classification"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18, margin = margin(b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray50", margin = margin(b = 20)),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 11),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 13, face = "bold"),
      strip.background = element_rect(fill = "gray90", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_equal()
  
  ggsave(file.path(figures_dir, "GAM_RF_Confusion_Matrix.pdf"), confusion_plot, width = 12, height = 6, dpi = 300, bg = "white")
  
  cat("✓ GAM RF confusion matrix saved\n")
  return(confusion_plot)
}

gam_rf_plot <- create_gam_rf_figures(results_dir_total, results_dir_overlap)

cat("\n=== Analysis Complete ===\n")