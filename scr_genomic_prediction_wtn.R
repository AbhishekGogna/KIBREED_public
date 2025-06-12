# Genomic Prediction for Winter Wheat - Within Environment G×E Analysis
# Prediction Type: WTN (Within environments)
# Models: 6 different G×E models with modular ETA components
# Cross-validation: CV3 scenario (genotype-environment combinations)

# Load required packages
library(BGLR)
library(qs)
library(dplyr)
library(AGHmatrix)

# =============================================================================
# CONFIGURATION PARAMETERS
# =============================================================================
if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}

PREDICTION_TYPE <- "wtn"
TEST_PROPORTION <- 0.33
RANDOM_SEED <- 123
THIN <- 5

# Subset parameters for debugging
N_GENO_SUBSET <- 1000 # Can be removed, but then script would need to be modified. 
N_ENV_SUBSET <- 10 # Can be removed, but then script would need to be modified. 
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

write_log("Starting G×E genomic prediction analysis for WTN", log_file, append = FALSE)

# =============================================================================
# DATA LOADING
# =============================================================================

# Load data
write_log("Loading data files...", log_file)

tryCatch({
  pheno_data <- qread(file.path(project_path, "data", paste0("pheno_", PREDICTION_TYPE, ".qs")))
  geno_add <- qread(file.path(project_path, "data", "g_data_add.qs"))
  geno_dom <- qread(file.path(project_path, "data", "g_data_dom.qs"))
  ev_data <- qread(file.path(project_path, "data", "ev_data.qs"))
  
  write_log(paste("Loaded phenotype data:", nrow(pheno_data), "observations"), log_file)
  write_log(paste("Loaded genomic data:", nrow(geno_add), "genotypes x", ncol(geno_add), "markers"), log_file)
}, error = function(e) {
  write_log(paste("ERROR loading data:", e$message), log_file)
  stop("Data loading failed")
})

# =============================================================================
# CROSS-VALIDATION SETUP
# =============================================================================

# CV3 cross-validation setup
write_log("Setting up CV3 cross-validation...", log_file)

# Get unique genotypes and environments
geno_ori <- unique(pheno_data$Geno_n)
env_ori <- unique(pheno_data$Env_n)

# Subset for debugging
set.seed(RANDOM_SEED)
sampled_geno <- sample(geno_ori, min(N_GENO_SUBSET, length(geno_ori)))
sampled_envs <- sample(env_ori, min(N_ENV_SUBSET, length(env_ori)))

# Create indexed phenotype data
pheno_data_idx <- pheno_data %>%
  filter(Env_n %in% sampled_envs, Geno_n %in% sampled_geno) %>%
  mutate(idx_cv = row_number()) 

# Genotype and environment lists
geno <- unique(pheno_data_idx$Geno_n)
ngeno <- length(geno)
env <- unique(pheno_data_idx$Env_n)
nenv <- length(env)

# Sample test genotypes
test_geno <- as.vector(geno[sample(1:ngeno, round(TEST_PROPORTION * ngeno))])

# Get test environments based on test genotype distribution
test_geno_dist <- pheno_data_idx %>% 
  filter(Geno_n %in% test_geno) %>% 
  count(Env_n) %>%
  arrange(desc(n)) %>%
  filter(n > (0.05 * length(test_geno)))

target_test_env <- test_geno_dist %>% pull(Env_n)
n_target_test_env <- length(target_test_env)
test_env <- as.vector(target_test_env[sample(1:n_target_test_env, min(round(TEST_PROPORTION * nenv), n_target_test_env))])

# Define training sets
train_geno <- setdiff(geno, test_geno)
train_env <- setdiff(env, test_env)

# Create CV3 split indices
train_indices <- pheno_data_idx %>%
  filter(Geno_n %in% train_geno & Env_n %in% train_env) %>%
  pull(idx_cv)

test_indices <- pheno_data_idx %>%
  filter(Geno_n %in% test_geno & Env_n %in% train_env) %>%
  pull(idx_cv)

# Prepare final dataset with CV splits
pheno_data_idx <- pheno_data_idx %>%
  mutate(
    set = case_when(
      idx_cv %in% train_indices ~ "train",
      idx_cv %in% test_indices ~ "test",
      TRUE ~ "unused"
    )
  ) %>%
  filter(set != "unused") %>%
  mutate(observed = ifelse(set == "train", BLUEs_wtn, NA)) %>%
  arrange(Env_n, Geno_n)

# Update genotype and environment lists
geno <- unique(pheno_data_idx$Geno_n)
env <- unique(pheno_data_idx$Env_n)

write_log(paste("Training set:", length(train_indices), "observations"), log_file)
write_log(paste("Test set:", length(test_indices), "observations"), log_file)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Helper functions
gaussian_kernel <- function(x, h = 1) {
  d <- as.matrix(dist(x, upper = TRUE, diag = TRUE))^2
  q <- median(d)
  return(exp(-h * d/q))
}

create_eigen_Z <- function(Z_base, K_matrix, threshold = 1e-8) {
  eigen_result <- eigen(K_matrix)
  keep_idx <- eigen_result$values > threshold
  eigen_vectors <- eigen_result$vectors[, keep_idx]
  eigen_values <- eigen_result$values[keep_idx]
  Z_eigen <- Z_base %*% sweep(eigen_vectors, 2, sqrt(eigen_values), FUN = '*')
  return(Z_eigen)
}

# =============================================================================
# CREATE ETA COMPONENTS
# =============================================================================

# Create all ETA components
write_log("Creating ETA components...", log_file)

# Basic incidence matrices
Z_e <- model.matrix(~ -1 + Env_n, data = pheno_data_idx)
Z_g <- model.matrix(~ -1 + Geno_n, data = pheno_data_idx)

# Genomic relationship matrices
geno_add_final <- geno_add[geno, ]
geno_dom_final <- geno_dom[geno, ]

# Additive kinship matrix
G_add <- Gmatrix(geno_add_final, method = "VanRaden", integer = FALSE)

# Dominance kinship matrix
YY_t <- tcrossprod(geno_dom_final)
G_dom <- YY_t / mean(diag(YY_t))

# Epistatic kinship matrix
G_epi <- G_add * G_add

# Environmental relationship matrices
env_matrix_erm <- tcrossprod(ev_data[env, ])
erm_l <- env_matrix_erm/(sum(diag(env_matrix_erm))/nrow(env_matrix_erm))

# Non-linear ERM using Gaussian kernel
erm_nl <- gaussian_kernel(x = ev_data[env, ])

# Create eigendecomposed Z matrices
Z_add <- create_eigen_Z(Z_g, G_add)
Z_dom <- create_eigen_Z(Z_g, G_dom)
Z_epi <- create_eigen_Z(Z_g, G_epi)
Z_erm_l <- create_eigen_Z(Z_e, erm_l)
Z_erm_nl <- create_eigen_Z(Z_e, erm_nl)

# G×E interaction kernels
GE_kernel_add_linear <- tcrossprod(Z_add) * tcrossprod(Z_erm_l)
GE_kernel_add_nonlinear <- tcrossprod(Z_add) * tcrossprod(Z_erm_nl)

# Define ETA components
ETA_components <- list(
  # Basic effects
  E_i = list(X = Z_e, model = "BRR"),
  G_i = list(X = Z_g, model = "BRR"),
  
  # Genomic effects (BRR)
  G_a = list(X = Z_add, model = "BRR"),
  G_d = list(X = Z_dom, model = "BRR"),
  G_aa = list(X = Z_epi, model = "BRR"),
  
  # Environmental effects (BRR)
  ERM_l = list(X = Z_erm_l, model = "BRR"),
  ERM_nl = list(X = Z_erm_nl, model = "BRR"),
  
  # G×E interactions (RKHS)
  G_a_ERM_l = list(K = GE_kernel_add_linear, model = "RKHS"),
  G_a_ERM_nl = list(K = GE_kernel_add_nonlinear, model = "RKHS")
)

# =============================================================================
# MODEL SPECIFICATIONS
# =============================================================================

# Model specifications
model_specifications <- list(
  M1 = list(
    eta = list(
      environment = ETA_components$E_i,
      genotype = ETA_components$G_i
    )
  ),
  
  M2 = list(
    eta = list(
      environment = ETA_components$E_i,
      additive = ETA_components$G_a,
      dominance = ETA_components$G_d,
      epistatic = ETA_components$G_aa
    )
  ),
  
  M3 = list(
    eta = list(
      environment_linear = ETA_components$ERM_l,
      additive = ETA_components$G_a,
      dominance = ETA_components$G_d,
      epistatic = ETA_components$G_aa
    )
  ),
  
  M4 = list(
    eta = list(
      environment_linear = ETA_components$ERM_l,
      additive = ETA_components$G_a,
      dominance = ETA_components$G_d,
      epistatic = ETA_components$G_aa,
      GE_add_linear = ETA_components$G_a_ERM_l
    )
  ),
  
  M5 = list(
    eta = list(
      environment_nonlinear = ETA_components$ERM_nl,
      additive = ETA_components$G_a,
      dominance = ETA_components$G_d,
      epistatic = ETA_components$G_aa
    )
  ),
  
  M6 = list(
    eta = list(
      environment_nonlinear = ETA_components$ERM_nl,
      additive = ETA_components$G_a,
      dominance = ETA_components$G_d,
      epistatic = ETA_components$G_aa,
      GE_add_nonlinear = ETA_components$G_a_ERM_nl
    )
  )
)

# Clean up - keep config, functions, final data, and model specifications
rm(list = setdiff(ls(), c("project_path", "PREDICTION_TYPE", "TEST_PROPORTION", "RANDOM_SEED", 
                          "N_ITER", "BURN_IN", "THIN", "N_GENO_SUBSET", 
                          "N_ENV_SUBSET", "log_file", "write_log",
                          "pheno_data_idx", "model_specifications")))

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
    filter(set == "test") %>%
    group_by(Series, Type) %>% 
    summarize(accuracy = cor(BLUEs_wtn, predicted, use = "complete.obs"), .groups = "drop") %>%
    group_by(Type) %>%
    summarize(accuracy_mean = mean(accuracy), .groups = "drop") %>%
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

# Run example models
write_log("Running models...", log_file)

# Run Model 1 (simplest)
model1_output <- fit_model(
  model_name = "M1",
  ETA_list = model_specifications$M1$eta,
  pheno_data_subset = pheno_data_idx,
  n_iter = N_ITER,
  burn_in = BURN_IN,
  thin = THIN,
  prediction_type = PREDICTION_TYPE,
  log_file_path = log_file,
  project_path = project_path
)

# Run Model 6 (most complex)
model6_output <- fit_model(
  model_name = "M6",
  ETA_list = model_specifications$M6$eta,
  pheno_data_subset = pheno_data_idx,
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
all_results <- rbind(model1_output$results, model6_output$results)
all_accuracy <- rbind(
  cbind(model1_output$accuracy, model = "M1_basic", prediction_type = PREDICTION_TYPE),
  cbind(model6_output$accuracy, model = "M6_complex", prediction_type = PREDICTION_TYPE)
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