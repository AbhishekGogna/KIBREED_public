#!/usr/bin/env Rscript

# R version 4.0.1
args = commandArgs(trailingOnly=TRUE)

# load arguments
cv_type = args[[1]]
cv_id = args[[2]] # index of the run info file
model_name = args[[3]]
model_spec = args[[4]]
in_cc = as.logical(args[[5]])
options(scipen = 999)

source("/proj/renv/activate.R")

# Functions
"%!in%" <- Negate("%in%")

load_mat <- function(type){
  #type can be one of "GRM_a", "GRM_d" or "GRM_d"
  if (!exists(type)) {
    mat_path <- get(paste0(type, "_path"))
    if(!is.null(mat_path) & file.exists(mat_path)) {
      mat <- qread(mat_path)
    } else {
      print("Either the path to GRM_a was not provided or no file exists at the provided path")
    }
  }
  return(mat)
}

# define variables and file paths
if (in_cc){
  ext_parse <- "results/R"
  write_at <- sprintf("/proj/results/R/generate_prediction_data/%s/run_data/%s", cv_type, cv_id)
} else { 
  ext_parse <- "ext_dir"
  write_at <- "/proj"
}

if (grepl("cv_acr_5f", cv_type)){
  idx_col <- "unique_idx"
  deselect <- "idx_with_series"
} else {
  idx_col <- "idx_with_series"
  deselect <- "unique_idx"
}

log_at <- sprintf("%s/logs/%s.log", write_at, model_name)
result_at <- sprintf("%s/preds/%s.qs", write_at, model_name)
dump_dir <- sprintf("%s/tmp_data/%s_dump", write_at, model_name)
tmp_data_dir <- sprintf("%s/tmp_data/%s_tmp", write_at, model_name)
for (dir in c(dump_dir, tmp_data_dir)) if(!dir.exists(dir)){dir.create(dir, recursive = T)}
dump_at <- paste0(dump_dir, "/", model_name, "_")

nIter <- 15000
burnIn <- 2000

# define needed packages
from_cran <- c("dplyr",
               "tidyr",
               "stringr",
               "BGLR",
               "logger",
               "qs",
               "jsonlite",
               "AGHmatrix") # done installation from rserver

if(sum(unlist(lapply(from_cran, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE))))) == length(from_cran)) {
  cat("All required packages are present and are loaded. Version check was not done.", file = log_at, sep = "\n")
} else {
  cat("Some packages were not loaded or are not installed. Please install and load packages manually", file = log_at, sep = "\n")
}

# register the log file for further logs
log_appender(appender_file(log_at))

# load data
run_data <- read_json(sprintf("/proj/%s/generate_prediction_data/%s/%s.json", ext_parse, cv_type, cv_type))
pheno_data <- qread(sprintf("/proj/%s/KIBREED_data_generation/BLUES_acr_env.qs", ext_parse))
G_a_path <- sprintf("/proj/%s/generate_prediction_data/%s/eigen_data/evd_a_mat.qs", ext_parse, cv_type)
G_aa_path <- sprintf("/proj/%s/generate_prediction_data/%s/eigen_data/evd_aa_mat.qs", ext_parse, cv_type)
G_d_path <- sprintf("/proj/%s/generate_prediction_data/%s/eigen_data/evd_d_mat.qs", ext_parse, cv_type)
log_info("All data loaded")

# generate run data
specific_run <- run_data[[cv_id]]
run_name <- paste0(cv_type, "_", cv_id , "_", model_name)
test_index <- do.call(c, specific_run$test)
train_index <- do.call(c, c(specific_run$train, specific_run$val))
all_data <- sort(c(test_index, train_index))
log_info("Run data defined")

# generate prediction specific data
p_data_0 <- pheno_data %>% mutate(obs = BLUEs) %>% rename(idx_col = idx_col) %>%
  select(-all_of(deselect)) %>%
  mutate(present_in = ifelse(idx_col %in% test_index, "inTest",
                             ifelse(idx_col %!in% all_data, "Nowhere", "inTrain"))) %>%
  distinct(idx_col, .keep_all = TRUE)

p_data <- p_data_0 %>%
  mutate(BLUEs = ifelse(present_in == "inTest" | present_in == "Nowhere", NA, BLUEs), # This is needed since i account for any genotype occuring in test set as a result of sampling scheme
         Type = as.character(Type)) %>%    # Convert 'Type' column to character type
  mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
  mutate(Type = as.factor(Type), # Convert 'Type' column to a factor
         run_name = run_name) %>%
  rename(series = Series,
         geno = Geno_new,
         type = Type,
         blues = BLUEs) %>%
  select(run_name, series, type, 
         geno, idx_col, 
         connect_geno_data,
         blues, obs, present_in)

geno_order <- p_data %>% pull(connect_geno_data) %>% as.character()
matrices <- unlist(lapply(str_split(str_split(model_spec, "&")[[1]], "@"), function(x) x[[1]]))

ETA <- list()
for(j in strsplit(model_spec, "&")[[1]]){
  mat_name = strsplit(j, "@")[[1]][1]
  model_type = strsplit(j, "@")[[1]][2]
  if(model_type == "BRR"){
    mat <- load_mat(mat_name)[["eigen"]]
    PC <- mat$vectors
    if(any(mat$values < 0)) log_info(sprintf("Eigen decomposition has some negative values for %s", mat_name))
    PC <- PC[, mat$values>1e-5] # has some negative values. check it.
    for(i in 1:ncol(PC)){PC[,i]=PC[,i]*sqrt(mat$values[i]) }
    ETA[[mat_name]] <- list(X = PC, model = model_type, saveEffects = TRUE)
    ETA[[mat_name]]$g_order <- identical(mat[["geno"]], geno_order)
  } else if (model_type == "RKHS"){
    mat <- load_mat(mat_name)[["eigen"]]
    ETA[[mat_name]] <- list(V=mat$vectors,d=mat$values,model=model_type, saveEffects = TRUE)
    ETA[[mat_name]]$g_order <- identical(mat[["geno"]], geno_order)
  }
}
log_info("ETA defined")

g_order_check <- sum(unlist(lapply(ETA, function(x) x[["g_order"]])))

# Do predictions
if(!file.exists(result_at)){
  if(g_order_check == length(matrices)){
    fmA <- BGLR(y=p_data$blues,
                ETA=ETA,
                nIter=nIter,
                burnIn=burnIn,
                saveAt=sprintf("%s%s_", dump_at, run_name),
                verbose=FALSE)
    log_info("Prediction model fitted")
    
    ## Extract vars and sd
    out_var_sd <- NULL
    out_var_sd_names <- NULL
    for (i in matrices){
      if (i %in% names(fmA[["ETA"]])){
        fit <- fmA[["ETA"]][[i]]
        interest_names <- grep("var\\S", names(fit), value = T, perl = T)
        sd_id <- grep("SD.", interest_names)
        sd <- fit[[interest_names[sd_id]]]
        var <- fit[[interest_names[-sd_id]]]
        out_var_sd <- c(out_var_sd, var, sd)
        out_var_sd_names <- c(out_var_sd_names, paste0(i, c("_var", "_sd")))
      } else {
        var <- NA
        sd <- NA
        out_var_sd <- c(out_var_sd, var, sd)
        out_var_sd_names <- c(out_var_sd_names, paste0(i, c("_var", "_sd")))
      }
    }
    out_var_sd <- c(out_var_sd, fmA$varE)
    combined_vars_sd <- as.data.frame(rbind(out_var_sd), row.names = run_name)
    colnames(combined_vars_sd) <-  c(out_var_sd_names, "Var.E")
    combined_vars_sd$run_name <- run_name
    p_data$pred <- as.numeric(fmA$yHat)

    # save results
    output <- list()
    output[["vars"]] <- combined_vars_sd
    output[["preds"]] <- p_data
    
    qsave(output, result_at)
    
    log_info("Results written")
  } else {
    log_info("Results writing failed since orde of phenotypes and genotypic data was different in atleast one of the input matrices. please check")
  }
} else {
  log_info("Result file was present, execution skipped")
}