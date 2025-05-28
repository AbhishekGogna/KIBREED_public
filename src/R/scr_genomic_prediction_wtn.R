#!/usr/bin/env Rscript

source("/proj/renv/activate.R")

# R version 4.0.5
args = commandArgs(trailingOnly=TRUE)

# load arguments
cv_type = args[[1]]
cv_id = args[[2]] # index of the run info file
model_name = args[[3]]
model_spec = args[[4]]
in_cc = as.logical(args[[5]])
options(scipen = 999)

# to dubug
#cv_type = "cv_wtn_tra" 
#cv_id = "run_1_cv1"
#model_name = "M_6"
#model_spec = "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS"
#in_cc = TRUE

# Functions
get_product <- function(incidence_mat, kin_mat, mat_names){
  kin_mat_ordered <- kin_mat[colnames(incidence_mat), colnames(incidence_mat)]
  t0 <- Sys.time()
  mat_eigen <- eigen(kin_mat_ordered)
  mat_eigen$vectors<-mat_eigen$vectors[,mat_eigen$values>1e-8]
  mat_eigen$values<-mat_eigen$values[mat_eigen$values>1e-8]
  mat_PC <- sweep(mat_eigen$vectors,MARGIN=2,STATS=sqrt(mat_eigen$values),FUN='*')
  # mat_PC<-t(apply(mat_eigen$vectors, 1, function(x) x*sqrt(mat_eigen$values)))
  out_mat <- incidence_mat %*% mat_PC
  if (!is.null(mat_names)){
    rownames(out_mat) <- mat_names
  }
  colnames(out_mat) <- paste0("var_", 1:ncol(out_mat))
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- out_mat
  out[["time_min"]] <- period
  return(out)
}

get_had_prod <- function(mat_1, mat_2){
  t0 <- Sys.time()
  mat_1_sym <- mat_1 %*% t(mat_1)
  mat_2_sym <- mat_2 %*% t(mat_2)
  mat_had <- mat_1_sym * mat_2_sym
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- mat_had
  out[["time_min"]] <- period
  
  return(out)
}

get_incidence_matrix <- function(data, column_name) {
  
  if (length(column_name) > 1){
    subset_0 <- paste(data[, column_name[1]], data[, column_name[2]], sep = ":")
    subset <- unique(subset_0)
    column_name <- paste0(column_name, collapse = ":")
  } else {
    subset <- NULL
  }
  
  formula <- sprintf("~ -1 + %s", column_name)
  t0 <- Sys.time()
  matrix <- model.matrix(as.formula(formula), data)
  
  
  if (is.null(subset)){
    colnames(matrix) <- gsub(column_name, '', colnames(matrix), perl = TRUE)  
  } else {
    part_1 <- gsub(strsplit(column_name, ":")[[1]][1], '', colnames(matrix), perl = TRUE)
    colnames(matrix) <- gsub(strsplit(column_name, ":")[[1]][2], '', part_1, perl = TRUE)
    matrix <- matrix[, subset]
  }
  
  t1 <- Sys.time()
  period <- round(as.numeric(t1 - t0, units = "mins"), 2)
  ## generate output
  out <- list()
  out[["mat"]] <- matrix
  out[["time_min"]] <- period
  return(out)
}

create_and_get_product <- function(col_name, pheno_data, kin_mat_path) {
  incidence_mat <- get_incidence_matrix(data = pheno_data,
                                        column_name = col_name)
  kin_mat <- qread(kin_mat_path)
  mat_out <- get_product(incidence_mat[["mat"]], kin_mat, NULL)
  return(mat_out)
}

wrapper_functions <- list(
  wrap_inc = function(pheno_data, matrix_params) {
    col_name <- matrix_params[["col_name"]]
    out_0 <- get_incidence_matrix(pheno_data, col_name)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR_kin = function(pheno_data, matrix_params) {
    col_name <- matrix_params[["col_name"]]
    kin_path <- paths[[matrix_params[["path"]]]]
    incidence_mat <- get_incidence_matrix(data = pheno_data,
                                          column_name = col_name)
    kin_mat <- qread(kin_path)
    out_0 <- get_product(incidence_mat[["mat"]], kin_mat, NULL)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS = function(pheno_data, matrix_params) {
    
    matrix_params <- matrix_params[1:2]
    
    matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
    matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
    
    matrix1 <- create_and_get_product(matrix_params_1[["col_name"]], pheno_data, paths[[matrix_params_1[["path"]]]])
    matrix2 <- create_and_get_product(matrix_params_2[["col_name"]], pheno_data, paths[[matrix_params_2[["path"]]]])
    
    out_0 <- get_had_prod(matrix1[["mat"]], matrix2[["mat"]])
    out <- list(K = out_0[["mat"]])
    return(out)
  }
  )

return_elements <- function(element, pheno_data, paths, 
                            tmp_data_dir, dump_dir, cv_id) {
  matrix_mapping <- list(
    "E_i@BRR"        = list(col_name = "connect_climate", 
                            type = "BRR", wrapper = "wrap_inc"),
    "G_i@BRR"        = list(col_name = "connect_geno", 
                            type = "BRR", wrapper = "wrap_inc"),
    "G_a@BRR"        = list(col_name = "connect_geno", path = "G_a_RM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_d@BRR"        = list(col_name = "connect_geno", path = "G_d_RM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_aa@BRR"       = list(col_name = "connect_geno", path = "G_aa_RM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "ERM_l@BRR"      = list(col_name = "connect_climate", path = "ERM_l",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "ERM_nl@BRR"     = list(col_name = "connect_climate", path = "ERM_nl",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "G_a_ERM_l@RKHS" = list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_a_RM", "ERM_l"), 
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_d_ERM_l@RKHS" = list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_d_RM", "ERM_l"),
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_a_ERM_nl@RKHS"= list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_a_RM", "ERM_nl"),
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_d_ERM_nl@RKHS"= list(col_name = c("connect_geno", "connect_climate"),
                            path = c("G_d_RM", "ERM_nl"),
                            type = "RKHS", wrapper = "wrap_RKHS")
  )
  
  matrix_params <- matrix_mapping[[element]] # get needed params
  
  wrapper_function <- wrapper_functions[[matrix_params[["wrapper"]]]]

  mat_out <- wrapper_function(pheno_data, matrix_params) # define paths explicitly. at the moment they are read in because of scope jump
  
  output_0 <- list(model = matrix_params[["type"]],
                   saveEffects = TRUE,
                   name = as.character(element))
  
  output <- c(output_0, mat_out)
  
  log_info(paste("ETA element for column =", element, "defined"))
  
  return(output)
}

# to_debug
if (in_cc){
  ext_parse <- "results/R"
  write_at <- sprintf("/proj/results/R/generate_prediction_data/%s/run_data/%s", cv_type, cv_id)
} else {
  stop()
}

# generate paths for files
log_at <- sprintf("%s/logs/%s.log", write_at, model_name)
result_at <- sprintf("%s/preds/%s.qs", write_at, model_name)
dump_dir <- sprintf("%s/tmp_data/%s_dump", write_at, model_name)
tmp_data_dir <- sprintf("%s/tmp_data/%s_tmp", write_at, model_name)

for (dir in c(dump_dir, tmp_data_dir)) if(!dir.exists(dir)){dir.create(dir, recursive = T)}
dump_at <- paste0(dump_dir, "/", model_name, "_")

# generate environment
from_cran <- c("dplyr",
               "tidyr",
               "hablar",
               "stringr",
               "BGLR",
               "logger",
               "qs",
               "jsonlite",
               "AGHmatrix",
               "feather") # done installation from rserver

if(sum(unlist(lapply(from_cran, function(x) suppressPackageStartupMessages(require(x, character.only = TRUE))))) == length(from_cran)) {
  cat("All required packages are present and are loaded. Version check was not done.", file = log_at, sep = "\n")
} else {
  cat("Some packages were not loaded or are not installed. Please install and load packages manually", file = log_at, sep = "\n")
}

# register the log file for further logs
log_appender(appender_file(log_at))

# load necessary data
log_info(sprintf("ext parse = %s", ext_parse))

kibreed_data_blues_within <- qread("/proj/data/pheno_data_wtn.qs")
run_data <- read_json(sprintf("/proj/%s/generate_prediction_data/%s/%s.json", ext_parse, cv_type, cv_type))

paths <- list(
  "G_a_RM" = "/proj/data/grm_a.qs", # based on genetic data
  "G_d_RM" = "/proj/data/grm_d.qs", # based on genetic data
  "G_aa_RM" = "/proj/data/grm_aa.qs", # based on genetic data
  "ERM_l" = "/proj/data/erm_l.qs", # based on environment data
  "ERM_nl" = "/proj/data/erm_nl.qs" # based on environment data
)

log_info("All data loaded")

# generate run data
specific_run <- run_data[[cv_id]]
run_name <- paste0(cv_type, "_", model_name, "_", cv_id)
test_index <- do.call(c, specific_run$test)
train_index <- do.call(c, c(specific_run$train, specific_run$val))
all_data <- sort(c(test_index, train_index))
log_info("Run data defined")

## set missing  
phenoGE_with_NA_subset <- kibreed_data_blues_within %>%
  filter(idx_cv %in% all_data) %>%
  transmute(
    run_name = run_name,
    env = Env,
    type = if_else(as.character(Type) == "Hybrid", "Hybrid", "Non_hybrid"),
    geno = Geno_new,
    connect_geno = connect_geno,
    connect_climate = Env,
    blues = if_else(idx_cv %in% test_index, NA_real_, BLUES_dt),
    obs = BLUES_dt
  ) %>%
  select(run_name, env, type, geno, connect_geno, connect_climate, blues, obs) %>%
  convert(fct(env),
          chr(type, run_name, starts_with("connect")),
          num(blues, obs))

## for debugging
##phenoGE_with_NA_subset <- phenoGE_with_NA_subset[sample(1:nrow(phenoGE_with_NA_subset), 1000), ]

#execute = TRUE # switch for quick debugging
#if (execute){
t0 = Sys.time()

if(!file.exists(result_at)){
  # define matrices
  elements <- str_split(model_spec, "&")[[1]]
  
  #mat_names <- paste0(phenoGE_with_NA_subset$env, ":", phenoGE_with_NA_subset$geno)
  #matrices <- unlist(lapply(str_split(elements, "@"), function(x) x[[1]]))
  
  ETA <- lapply(elements, return_elements,
                pheno_data = phenoGE_with_NA_subset, 
                paths = paths,
                tmp_data_dir = tmp_data_dir,
                dump_dir = dump_dir,
                cv_id = cv_id)
  names(ETA) <- sapply(ETA, function(x) x[["name"]])
  
  t1 = Sys.time()
  log_info(sprintf("ETAs defined for model %s i.e. %s . The whole process starting with defining of E_i Took %s min",
                   model_name, 
                   model_spec,
                   round(as.numeric(t1 - t0, units = "mins"), 2)))
  
  # check if any element has missing values 
  check_na <- sapply(ETA, function(x) ifelse(sum(is.na(x[[grep("X|K", names(x), value = T)]])) == 0, TRUE, FALSE))
  
  if(any(check_na)){
    # run model
    #RhpcBLASctl::blas_set_num_threads(1)
    fm_base <- BGLR(y = phenoGE_with_NA_subset$blues,
                    ETA = ETA,
                    nIter = 15000,
                    burnIn = 2000,
                    thin = 5,
                    saveAt = dump_at,
                    verbose = FALSE)
    t2 <- Sys.time()
    log_info(sprintf("%s took %s minutes.\nPredictions finished.",
                     model_name, 
                     round(as.numeric(t2 - t1, units = "mins"), 2)))
    
    # extract outputs
    
    ## Yield
    phenoGE_with_NA_subset$pred <- as.numeric(fm_base$yHat)
    
    ## vars and sd
    out_var_sd <- NULL
    out_var_sd_names <- NULL
    for (i in elements){
      if (i %in% names(fm_base[["ETA"]])){
        fit <- fm_base[["ETA"]][[i]]
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
    out_var_sd <- c(out_var_sd, fm_base$varE)
    out_var_sd_names <- c(out_var_sd_names, "err_var")
    
    combined_vars_sd <- as.data.frame(rbind(out_var_sd), row.names = run_name)
    colnames(combined_vars_sd) <- out_var_sd_names
    
    # save results
    output <- list()
    output[["vars"]] <- combined_vars_sd
    output[["preds"]] <- phenoGE_with_NA_subset
    
    qsave(output, result_at)
    log_info("Results written")
  } else {
    log_info("Predictors have missing values, execution skipped")
    log_info(paste(names(check_na), check_na, sep = ": ", collapse = ", "))
  }
} else {
  log_info("Result file was present, execution skipped")
}
