
## 01_ModelTesting.R

**Purpose:** Tests machine learning models (Random Forest, SVM, KNN) on watershed classification tasks using preprocessed fish data.

**What it does:**
- Loads pre-split training and testing datasets (RAW, GAM, and MA data types)
- Trains three model types (Random Forest, SVM, K-Nearest Neighbors) on each dataset
- Evaluates model performance on both "Total" and "Overlapping" datasets
- Saves trained models and prediction results

**Output files:**
- Trained models (`.rds` files) in `Output/Models/`
- Prediction CSV files in `Output/ModelResultsPreCal/`
- **Combined_Accuracy_Heatmap.pdf** - Visual comparison of model accuracy
- **Combined_F1Score_Heatmap.pdf** - Visual comparison of F1-scores
- **GAM_RF_Confusion_Matrix.pdf** - Confusion matrices for GAM Random Forest model

---

## 02_PCA_Figures.R

**Purpose:** Creates Principal Component Analysis (PCA) visualizations to explore variation in fish otolith isotope data.

**What it does:**
- Performs PCA on GAM-smoothed and raw Sr isotope time series data
- Creates summary PCA plots showing separation between watersheds
- Generates individual fish time series plots with PC loadings overlay
- Produces 2D PCA comparison plots (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3)

**Output files:**
- **Four_Panel_PCA_Loadings_Comparison.pdf** - Multi-panel figure showing PC loadings for selected fish
- **2D PCA plots** (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3) showing watershed separation
- Various PCA summary plots for different data subsets

---

## 03_CWAK_Analysis.R

**Purpose:** Analyzes Random Forest model performance specifically for CWAK (Coastal Western Alaska) vs non-CWAK fish populations.

**What it does:**
- Loads test fish predictions from the GAM Random Forest model
- Merges predictions with genetic data to classify Yukon fish as CWAK or non-CWAK
- Groups fish into CWAK (Nush, Kusko, Lower Yukon) and non-CWAK (Middle/Upper Yukon)
- Calculates classification accuracy for each group

**Output files:**
- **CWAK_RF_Performance_Summary.csv** - Overall accuracy by CWAK group
- **CWAK_RF_Detailed_Performance.csv** - Breakdown by watershed
- **CWAK_RF_Yukon_Genetic_Breakdown.csv** - Performance by Yukon genetic groupings
- **RF_Predictions_with_CWAK_Groups.csv** - Full dataset with CWAK classifications

---

## 04_CombinedCalibration.R

**Purpose:** Calibrates Random Forest model prediction probabilities and evaluates performance across different confidence thresholds.

**What it does:**
- Applies isotonic regression calibration to RF model probabilities
- Compares before/after calibration using log-loss and Brier scores
- Tests classification performance at various probability thresholds (60-90%)
- Creates line plots showing accuracy vs. threshold for each watershed

**Output files:**
- **RF_Calibration_[DataType]_[Analysis].png** - Before/after calibration plots for each model
- **RF_Calibration_Summary.csv** - Summary of calibration improvements
- **GAM_RF_TOTAL_Performance_Line_Plot.png** - Threshold performance for total dataset
- **GAM_RF_Restricted_Performance_Line_Plot.png** - Threshold performance for restricted dataset
- **GAM_RF_Combined_Performance_Line_Plot.png** - Side-by-side comparison of both datasets

# Files 

NatalToMarine_Processed_GAM.csv: trimmed LA-ICPMS data gam smoothed 
NatalToMarine_Processed_RAW.csv: timmed LA-ICPMS data raw data (no ts smoothing) 
All_Yukon_Genetics.csv: Genetic posterior assignments for all yukon individuals following methods in Makhlouf et al., 2025 
Fish_ID_Splits.csv: Fish ID splits for training and testing sets, used to make sure the same individuals are in all testing/train sets. 
Test_GAM.csv: Testing dataset of GAM smoothed ts 
Test_MA.csv: Testing dataset of Moving average smoothed ts 
Test_RAW.csv: Testing dataset of raw (non smoothed) ts 
Train_GAM.csv: Training dataset of GAM smoothed ts 
Train_MA.csv: Training dataset of MA smoothed ts 
Train_RAW.csv: Training dataset of raw (non_smoothed) ts 


