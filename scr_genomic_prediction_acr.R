# Genomic Prediction for Winter Wheat - Across Environments Analysis
# Prediction Type: ACR (Across environments)
# Models: Additive+Dominance and Additive+Dominance+Epistatic
# Cross-validation: 80/20 random split

# Load required packages
library(BGLR)
library(qs)
library(dplyr)
library(AGHmatrix)
library(feather)

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}

PREDICTION_TYPE <- "acr"  # across environments
TEST_PROPORTION <- 0.2
RANDOM_SEED <- 123
THIN <- 5

# Subset parameters for debugging
N_GENO_SUBSET <- 1000 # Can be removed, but then script would need to be modified. 
N_ITER <- 150  # can raise to 15000
BURN_IN <- 20  # can raise to 2000

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================

# Create directories relative to project_path
dirs_to_create <- c("results", "tmp", "logs")
for (dir in dirs_to_create) {
  dir_path <- file.path(project_path, dir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
}

# Set up logging
log_file <- file.path(project_path, "logs", paste0("genomic_prediction_", PREDICTION_TYPE, ".log"))

write_log <- function(message, log_file_path, append = TRUE) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste(timestamp, "-", message)
  cat(log_entry, "\n", file = log_file_path, append = append)
  cat(log_entry, "\n")
}

write_log("Starting G×E genomic prediction analysis for ACR", log_file, append = FALSE)

# =============================================================================
# DATA LOADING
# =============================================================================

write_log("Loading data files...", log_file)

tryCatch({
  pheno_data <- qread(file.path(project_path, "data", paste0("pheno_", PREDICTION_TYPE, ".qs")))
  geno_add <- qread(file.path(project_path, "data", "g_data_add.qs"))
  geno_dom <- qread(file.path(project_path, "data", "g_data_dom.qs"))
  
  write_log(paste("Loaded phenotype data:", nrow(pheno_data), "observations"), log_file)
  write_log(paste("Loaded genomic data:", nrow(geno_add), "genotypes x", ncol(geno_add), "markers"), log_file)
}, error = function(e) {
  write_log(paste("ERROR loading data:", e$message), log_file)
  stop("Data loading failed")
})

# =============================================================================
# DATA PREPARATION
# =============================================================================

write_log("Preparing phenotype data...", log_file)

# Subset for debugging
geno_ori <- unique(pheno_data$Geno_n)
set.seed(RANDOM_SEED)
sampled_geno <- sample(geno_ori, min(N_GENO_SUBSET, length(geno_ori)))

pheno_clean <- pheno_data %>%
  distinct(Geno_n, BLUEs_acr, Type) %>%
  filter(!is.na(BLUEs_acr)) %>%
  filter(Geno_n %in% sampled_geno) %>%
  mutate(idx_cv = row_number())

geno_add_final <- geno_add[pheno_clean$Geno_n, ]
geno_dom_final <- geno_dom[pheno_clean$Geno_n, ]

write_log(paste("Final dataset:", nrow(pheno_clean), "genotypes"), log_file)

# =============================================================================
# CROSS-VALIDATION SETUP
# =============================================================================

write_log("Setting up cross-validation...", log_file)

set.seed(RANDOM_SEED)
n_total <- nrow(pheno_clean)
test_size <- round(n_total * TEST_PROPORTION)
test_indices <- sample(pheno_clean$idx_cv, test_size)
train_indices <- setdiff(pheno_clean$idx_cv, test_indices)

# Add to pheno data
pheno_final <- pheno_clean %>%
  mutate(set = ifelse(idx_cv %in% train_indices, "Train", "Test"),
         observed = ifelse(idx_cv %in% train_indices, BLUEs_acr, NA))

write_log(paste("Training set:", length(train_indices), "observations"), log_file)
write_log(paste("Test set:", length(test_indices), "observations"), log_file)

# =============================================================================
# EXPORT FILES FOR PYTHON INTEGRATION
# =============================================================================

write_log("Exporting files for Python integration...", log_file)

# Export genomic data
geno_df <- data.frame(
  idx = rownames(geno_add_final),
  geno_add_final,
  row.names = NULL
)

write_feather(geno_df, 
              file.path(project_path, "results", "g_data_add_acr.feather"))

# Export phenotype data
write_feather(pheno_final, file.path(project_path, "results", "pheno_final_acr.feather"))

write_log("Files exported successfully for Python integration", log_file)

# =============================================================================
# KINSHIP MATRICES
# =============================================================================

write_log("Creating kinship matrices...", log_file)

# Additive kinship matrix using VanRaden method
G_add <- Gmatrix(geno_add_final, method = "VanRaden", integer = FALSE)

# Dominance kinship matrix
YY_t <- tcrossprod(geno_dom_final)
G_dom <- YY_t / mean(diag(YY_t))

# Epistatic kinship matrix (additive x additive)
G_epi <- G_add * G_add

# Clean up - keep config, functions, final data, and kinship matrices
rm(list = setdiff(ls(), c("project_path", "PREDICTION_TYPE", "RANDOM_SEED", 
                          "N_ITER", "BURN_IN", "THIN", 
                          "log_file", "write_log",
                          "pheno_final", "G_add", "G_dom", "G_epi")))

# =============================================================================
# MODEL FITTING FUNCTION
# =============================================================================

# Model fitting function
fit_model <- function(ETA_list, model_name, pheno_data_subset, 
                      n_iter, burn_in, thin, prediction_type, log_file_path, project_path) {
  
  write_log(paste("Fitting", model_name), log_file_path)
  
  t0 <- Sys.time()
  model <- BGLR(
    y = pheno_data_subset$observed,
    ETA = ETA_list,
    nIter = n_iter,
    burnIn = burn_in,
    thin = thin,
    saveAt = file.path(project_path, "tmp", paste0(prediction_type, "_", model_name, "_")),
    verbose = FALSE
  )
  t1 <- Sys.time()
  
  # Extract results
  results <- pheno_data_subset %>%
    mutate(predicted = model$yHat,
           model_name = model_name)
  
  # Calculate accuracy by type
  accuracy <- results %>% 
    filter(set == "Test") %>%
    group_by(Type) %>% 
    summarize(accuracy = cor(BLUEs_acr, predicted, use = "complete.obs"), .groups = "drop") %>%
    as.data.frame()
  
  runtime <- round(as.numeric(t1 - t0, units = "mins"), 2)
  write_log(paste(model_name, "completed in", runtime, "minutes"), log_file_path)
  
  for (i in 1:nrow(accuracy)) {
    write_log(paste(model_name, accuracy$Type[i], "accuracy:", round(accuracy$accuracy[i], 3)), log_file_path)
  }
  
  return(list(results = results, accuracy = accuracy, model = model))
}

# =============================================================================
# RUN MODELS
# =============================================================================

write_log("Running models...", log_file)

# Model 1: Additive + Dominance
ETA_model1 <- list(
  additive = list(K = G_add, model = "RKHS"),
  dominance = list(K = G_dom, model = "RKHS")
)

model1_output <- fit_model(
  model_name = "M1",
  ETA = ETA_model1,
  pheno_data_subset = pheno_final,
  n_iter = N_ITER,
  burn_in = BURN_IN,
  thin = THIN,
  prediction_type = PREDICTION_TYPE,
  log_file_path = log_file,
  project_path = project_path
)

# Model 2: Additive + Dominance + Epistatic
ETA_model2 <- list(
  additive = list(K = G_add, model = "RKHS"),
  dominance = list(K = G_dom, model = "RKHS"),
  epistatic = list(K = G_epi, model = "RKHS")
)

model2_output <- fit_model(
  model_name = "M2",
  ETA = ETA_model2,
  pheno_data_subset = pheno_final,
  n_iter = N_ITER,
  burn_in = BURN_IN,
  thin = THIN,
  prediction_type = PREDICTION_TYPE,
  log_file_path = log_file,
  project_path = project_path
)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Combine all results
all_results <- rbind(model1_output$results, model2_output$results)
all_accuracy <- rbind(
  cbind(model1_output$accuracy, model = "M1_basic", prediction_type = PREDICTION_TYPE),
  cbind(model2_output$accuracy, model = "M2_complex", prediction_type = PREDICTION_TYPE)
)

results_file <- file.path(project_path, "results", paste0("prediction_results_", PREDICTION_TYPE, ".qs"))
accuracy_file <- file.path(project_path, "results", paste0("accuracy_summary_", PREDICTION_TYPE, ".qs"))

qsave(all_results, results_file)
qsave(all_accuracy, accuracy_file)

write_log("Analysis completed", log_file)

# Clean up temporary files
temp_files <- list.files(file.path(project_path, "tmp"), pattern = paste0(PREDICTION_TYPE, "_"), full.names = TRUE)
if (length(temp_files) > 0) {
  file.remove(temp_files)
}

# Final cleanup - keep only essential results for potential further analysis
rm(list = setdiff(ls(), c("all_results", "all_accuracy", 
                          "results_file", "accuracy_file", "log_file")))
