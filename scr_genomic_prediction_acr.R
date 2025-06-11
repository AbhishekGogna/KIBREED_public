# Genomic Prediction for Winter Wheat - Across Environments Analysis
# Prediction Type: ACR (Across environments)
# Models: Additive+Dominance and Additive+Dominance+Epistatic
# Cross-validation: 80/20 random split

# Load required packages
library(BGLR)
library(qs)
library(dplyr)
library(AGHmatrix)

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================

PREDICTION_TYPE <- "acr"  # across environments
TEST_PROPORTION <- 0.2
RANDOM_SEED <- 123
N_ITER <- 150  # can raise to 15000 for production
BURN_IN <- 20  # can raise to 2000 for production
THIN <- 5

# Subset parameters for debugging
N_GENO_SUBSET <- 1000

# =============================================================================
# SETUP AND INITIALIZATION
# =============================================================================

# Create directories
dirs_to_create <- c("results", "tmp", "logs")
for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

# Set up logging
log_file <- paste0("logs/genomic_prediction_", PREDICTION_TYPE, ".log")

write_log <- function(message, log_file_path) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- paste(timestamp, "-", message)
  cat(log_entry, "\n", file = log_file_path, append = TRUE)
  cat(log_entry, "\n")
}

write_log("Starting genomic prediction analysis for ACR", log_file)

# =============================================================================
# DATA LOADING
# =============================================================================

write_log("Loading data files...", log_file)

tryCatch({
  pheno_data <- qread(paste0("data/pheno_", PREDICTION_TYPE, ".qs"))
  geno_add <- qread("data/g_data_add.qs")
  geno_dom <- qread("data/g_data_dom.qs")
  
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
  mutate(Type = ifelse(Type == "Hybrid", "Hybrid", "Line"))

# Match genotypes between phenotype and genomic data
common_genos <- intersect(pheno_clean$Geno_n, rownames(geno_add))

if (length(common_genos) == 0) {
  stop("No common genotypes found between phenotype and genomic data")
}

# Subset data to common genotypes
pheno_final <- pheno_clean %>%
  filter(Geno_n %in% common_genos) %>%
  arrange(Geno_n)

geno_add_final <- geno_add[pheno_final$Geno_n, ]
geno_dom_final <- geno_dom[pheno_final$Geno_n, ]

write_log(paste("Final dataset:", nrow(pheno_final), "genotypes"), log_file)

# =============================================================================
# CROSS-VALIDATION SETUP
# =============================================================================

write_log("Setting up cross-validation...", log_file)

set.seed(RANDOM_SEED)
n_total <- nrow(pheno_final)
test_size <- round(n_total * TEST_PROPORTION)
test_indices <- sample(1:n_total, test_size)
train_indices <- setdiff(1:n_total, test_indices)

# Prepare response variable
y <- pheno_final$BLUEs_acr
y_train <- y
y_train[test_indices] <- NA  # Set test observations to NA

write_log(paste("Training set:", length(train_indices), "observations"), log_file)
write_log(paste("Test set:", length(test_indices), "observations"), log_file)

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

# =============================================================================
# MODEL FITTING FUNCTION
# =============================================================================

fit_and_extract <- function(ETA, model_name, model_label, y_response, pheno_data_final, 
                            test_idx, n_total_obs, n_iter, burn_in, thin, 
                            prediction_type, log_file_path) {
  
  write_log(paste("Fitting", model_name), log_file_path)
  
  t0 <- Sys.time()
  model <- BGLR(
    y = y_response,
    ETA = ETA,
    nIter = n_iter,
    burnIn = burn_in,
    thin = thin,
    saveAt = paste0("tmp/", prediction_type, "_", model_label, "_"),
    verbose = FALSE
  )
  t1 <- Sys.time()
  
  # Extract results
  results <- data.frame(
    genotype = pheno_data_final$Geno_n,
    type = pheno_data_final$Type,
    observed = pheno_data_final$BLUEs_acr,
    predicted = model$yHat,
    dataset = ifelse(1:n_total_obs %in% test_idx, "test", "train"),
    model = model_label,
    prediction_type = prediction_type
  )
  
  # Calculate accuracy for test set only
  test_accuracy <- results %>% 
    filter(dataset == "test") %>%
    group_by(type) %>% 
    summarize(accuracy = cor(observed, predicted, use = "complete.obs"), .groups = "drop") %>%
    as.data.frame()
  
  runtime <- round(as.numeric(t1 - t0, units = "mins"), 2)
  write_log(paste(model_name, "completed in", runtime, "minutes"), log_file_path)
  
  for (i in 1:nrow(test_accuracy)) {
    write_log(paste(model_name, test_accuracy$type[i], "accuracy:", round(test_accuracy$accuracy[i], 3)), log_file_path)
  }
  
  return(list(results = results, accuracy = test_accuracy, model = model))
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

model1_output <- fit_and_extract(
  ETA = ETA_model1,
  model_name = "Model 1: Additive + Dominance",
  model_label = "add_dom",
  y_response = y_train,
  pheno_data_final = pheno_final,
  test_idx = test_indices,
  n_total_obs = n_total,
  n_iter = N_ITER,
  burn_in = BURN_IN,
  thin = THIN,
  prediction_type = PREDICTION_TYPE,
  log_file_path = log_file
)

# Model 2: Additive + Dominance + Epistatic
ETA_model2 <- list(
  additive = list(K = G_add, model = "RKHS"),
  dominance = list(K = G_dom, model = "RKHS"),
  epistatic = list(K = G_epi, model = "RKHS")
)

model2_output <- fit_and_extract(
  ETA = ETA_model2,
  model_name = "Model 2: Additive + Dominance + Epistatic",
  model_label = "add_dom_epi",
  y_response = y_train,
  pheno_data_final = pheno_final,
  test_idx = test_indices,
  n_total_obs = n_total,
  n_iter = N_ITER,
  burn_in = BURN_IN,
  thin = THIN,
  prediction_type = PREDICTION_TYPE,
  log_file_path = log_file
)

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Combine all results
all_results <- rbind(model1_output$results, model2_output$results)
all_accuracy <- rbind(
  cbind(model1_output$accuracy, model = "add_dom", prediction_type = PREDICTION_TYPE),
  cbind(model2_output$accuracy, model = "add_dom_epi", prediction_type = PREDICTION_TYPE)
)

# Save results
results_file <- paste0("results/prediction_results_", PREDICTION_TYPE, ".qs")
accuracy_file <- paste0("results/accuracy_summary_", PREDICTION_TYPE, ".qs")

qsave(all_results, results_file)
qsave(all_accuracy, accuracy_file)

write_log("Analysis completed", log_file)

# Clean up temporary files
temp_files <- list.files("tmp", pattern = paste0(PREDICTION_TYPE, "_"), full.names = TRUE)
if (length(temp_files) > 0) {
  file.remove(temp_files)
}