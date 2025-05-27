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

# functions
ifrm <- function(obj, env = globalenv()) {
  obj <- deparse(substitute(obj))
  if(exists(obj, envir = env)) {
    rm(list = obj, envir = env)
  }
}

prep_as_df_for_feather <- function(data){
  data_df <- cbind("idx" = rownames(data), data)
  data_df <- as.data.frame(data_df)
  return(data_df)
}

get_product <- function(incidence_mat, kin_mat, mat_names, save_at){
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
  
  if (!is.null(save_at)){
    if(!file.exists(save_at)) {
      out_df <- prep_as_df_for_feather(out[["mat"]])
      write_feather(out_df, save_at)
    }
  }
  
  return(out)
}

get_had_prod <- function(mat_1, mat_2, save_at){
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
  
  if (!is.null(save_at)){
    if(!file.exists(save_at)) {
      out_df <- prep_as_df_for_feather(out[["mat"]])
      write_feather(out_df, save_at)
    }
  }
  
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

get_BRR_pred <- function(pheno_data, col_name, kin_path, var_path, 
                         tmp_data_dir, dump_dir, 
                         cv_id, nIter, burnIn){
  
  # create tmp dirs
  dump_at_temp <- paste0(dump_dir, "/par_pred")
  if(!dir.exists(dump_at_temp)) dir.create(dump_at_temp, recursive = T)
  
  # create tmp_dir
  tmp_data <- paste0(tmp_data_dir, "/par_pred")
  if(!dir.exists(tmp_data)) dir.create(tmp_data)
  
  pheno_data_0 <- pheno_data %>%
    mutate(idx = row_number())
  
  g_data <- qread(kin_path)
  param_data <- qread(var_path)
  param_data_df <- as.data.frame(param_data, row.names = NA)
  param_data_df[, col_name] <- rownames(param_data)
  
  pheno_data_mod <- pheno_data_0 %>% mutate(obs_set = ifelse(is.na(blues), "test", "train")) %>%
    arrange(desc(obs_set)) %>% left_join(param_data_df, by = col_name)
  
  envs <- as.character(unique(pheno_data_mod$env))
  meta_info <- NULL
  param_data_out <- NULL
  
  for (env in envs){
    ifrm(pheno_data_0_site, env = environment())
    pheno_data_0_site <- pheno_data_0 %>% filter(env == !!as.character(env))
    param_data_site <- pheno_data_mod %>% filter(env == !!as.character(env)) %>%
      mutate(across(all_of(colnames(param_data)), ~ifelse(!is.na(blues), ., NA), .names = "to_pred_{.col}"))
    geno_site <- param_data_site %>% pull(connect_geno) %>% as.character()
    geno_data_site <- g_data[geno_site, geno_site]
    for (param in colnames(param_data)){
      cols_to_select <- c("idx", colnames(pheno_data), "obs_set", param, paste0("to_pred_", param))
      ifrm(param_data_site_param, env = environment())
      param_data_site_param <- param_data_site %>% select(all_of(cols_to_select)) %>%
        rename("param_raw" = all_of(param),
               "param" = all_of(paste0("to_pred_", param)))
      model_fit <- BGLR(y = param_data_site_param[, "param"],
                        ETA = list(G=list(K=geno_data_site,model='RKHS')),
                        nIter = nIter,
                        burnIn = burnIn,
                        thin = 5,
                        saveAt = sprintf("%s/%s_%s_%s_", dump_at_temp, env, param, cv_id),
                        verbose = FALSE)
      param_data_site_param$pred_param <- model_fit$yHat
      
      meta <- data.frame("env" = env,
                         "param" = param,
                         "cv_id" = cv_id,
                         "var_G" = model_fit$ETA$G$varU,
                         "var_E" = model_fit$varE,
                         "rep" = model_fit$ETA$G$varU/(model_fit$ETA$G$varU+model_fit$varE),
                         "train_set_cor_avg" = param_data_site_param %>% 
                           filter(obs_set == "train") %>% group_by(env) %>% 
                           summarize(corr_v = cor(param_raw, pred_param)) %>% ungroup() %>% 
                           pull(corr_v) %>% as.numeric() %>% mean(),
                         "test_set_cor_avg" =  param_data_site_param %>% 
                           filter(obs_set == "test") %>% group_by(env) %>% 
                           summarize(corr_v = cor(param_raw, pred_param)) %>% ungroup() %>% 
                           pull(corr_v) %>% as.numeric() %>% mean())
      
      meta_info <- rbind(meta_info, meta)
      
      param_data_site_param <- param_data_site_param %>% mutate(param = ifelse(is.na(param), pred_param, param)) %>%
        select(-param_raw, -pred_param, -obs_set) %>%
        relocate(param, .after = last_col())
      colnames(param_data_site_param)[ncol(param_data_site_param)] <- param
      pheno_data_0_site <- pheno_data_0_site %>% left_join(param_data_site_param, by = c("idx", colnames(pheno_data))) # binds columns
    }
    if(is.null(param_data_out)) {
      param_data_out <- pheno_data_0_site
    } else {
      param_data_out <- param_data_out %>% bind_rows(pheno_data_0_site) # binds rows
    }
  }
  
  # remove temp_data
  system(sprintf("rm -rf %s", dump_at_temp))
  
  # derive_output
  param_data_out_mat_0 <- param_data_out %>% select(-all_of(c(colnames(pheno_data), "idx")))
  param_data_out_mat <- as.matrix(param_data_out_mat_0)
  rownames(param_data_out_mat) <- param_data_out[, col_name]
  out <- list("mat" = param_data_out_mat)
  
  # write_data
  write_feather(param_data_out, sprintf("%s/%s.feather", tmp_data, cv_id))
  write.table(meta_info, sprintf("%s/meta_%s.csv", tmp_data, cv_id))
  
  return(out)
}

get_mat <- function(pheno_data, col_name, BRR_mat){
  
  # calulate kinship
  param_data_df <- as.data.frame(BRR_mat, row.names = NA)
  param_data_df[, col_name] <- rownames(BRR_mat)
  
  pheno_data_mod <- pheno_data %>% 
    mutate(obs_set = ifelse(is.na(blues), "test", "train")) %>%
    arrange(desc(obs_set)) %>% left_join(param_data_df, by = col_name) %>% 
    filter(obs_set == "train")
  param_data_kin <- pheno_data_mod %>% 
    select(all_of(colnames(BRR_mat)), "latlong", "geno") %>% 
    group_by(latlong) %>% 
    summarize(across(all_of(colnames(BRR_mat)), ~ mean(as.numeric(.x), na.rm = TRUE)))
  
  param_data_kin_scaled <- scale(param_data_kin[, -1], scale = T, center = T)
  rownames(param_data_kin_scaled) <- param_data_kin %>% pull(latlong) %>% as.vector()
  param_data_kin_scaled_miss <- as.matrix(param_data_kin_scaled)
  
  RM <- tcrossprod(param_data_kin_scaled_miss)
  RM_scaled <- RM/(sum(diag(RM))/nrow(RM))
  
  # calculate incidence
  inc_mat <- get_incidence_matrix(pheno_data, "latlong")[["mat"]]
  
  # generate output
  if(any(colnames(RM_scaled) %in% colnames(inc_mat))) {
    mat_out <- get_product(inc_mat, RM_scaled, NULL, NULL)
  } else {
    mat_out <- NULL
  }
  return(mat_out)
}

create_and_get_product <- function(col_name, pheno_data, kin_mat_path) {
  incidence_mat <- get_incidence_matrix(data = pheno_data,
                                        column_name = col_name)
  kin_mat <- qread(kin_mat_path)
  mat_out <- get_product(incidence_mat[["mat"]], kin_mat, NULL, NULL)
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
    out_0 <- get_product(incidence_mat[["mat"]], kin_mat, NULL, NULL)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS = function(pheno_data, matrix_params) {
    
    matrix_params <- matrix_params[1:2]
    
    matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
    matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
    
    matrix1 <- create_and_get_product(matrix_params_1[["col_name"]], pheno_data, paths[[matrix_params_1[["path"]]]])
    matrix2 <- create_and_get_product(matrix_params_2[["col_name"]], pheno_data, paths[[matrix_params_2[["path"]]]])
    
    out_0 <- get_had_prod(matrix1[["mat"]], matrix2[["mat"]], NULL)
    out <- list(K = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR = function(pheno_data, matrix_params) {
    col_name <- matrix_params[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params[["path"]]]])
    out <- list(X = BRR_mat[pheno_data[, col_name], ])
    return(out)
  },
  
  wrap_BRR_2 = function(pheno_data, matrix_params) {
    col_name <- matrix_params[["col_name"]]
    kin_path <- paths[[matrix_params[["path"]][1]]]
    var_path <- paths[[matrix_params[["path"]][2]]]
    tmp_data_dir <- matrix_params[["tmp_data_dir"]]
    dump_dir <- matrix_params[["dump_dir"]]
    cv_id  <- matrix_params[["cv_id"]]
    out_0 <- get_BRR_pred(pheno_data, col_name, kin_path, var_path, 
                          tmp_data_dir, dump_dir, 
                          cv_id, nIter = 15000, burnIn = 2000)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR_3 = function(pheno_data, matrix_params) {
    col_name <- matrix_params[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params[["path"]]]])
    out_0 <- get_mat(pheno_data, col_name, BRR_mat)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS_1 = function(pheno_data, matrix_params) {
    
    matrix_params <- matrix_params[1:2]
    
    matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
    matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
    
    matrix1 <- create_and_get_product(matrix_params_1[["col_name"]], pheno_data, paths[[matrix_params_1[["path"]]]])
    
    col_name <- matrix_params_2[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params_2[["path"]]]])
    matrix2 <- get_mat(pheno_data, col_name, BRR_mat) 
    
    out_0 <- get_had_prod(matrix1[["mat"]], matrix2[["mat"]], NULL)
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
    "S_i@BRR"        = list(col_name = "site", 
                            type = "BRR", wrapper = "wrap_inc"),
    "Y_i@BRR"        = list(col_name = "year", 
                            type = "BRR", wrapper = "wrap_inc"),
    "G_a_S_i@BRR"    = list(col_name = c("connect_geno", "latlong"), 
                            type = "BRR", wrapper = "wrap_inc"),
    "S@BRR"          = list(col_name = "site", path = "SRM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
    "Y@BRR"          = list(col_name = "year", path = "YRM",
                            type = "BRR", wrapper = "wrap_BRR_kin"),
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
                            type = "RKHS", wrapper = "wrap_RKHS"),
    "G_a_Y@RKHS"      = list(col_name = c("connect_geno", "year"), 
                             path = c("G_a_RM", "YRM"),  
                             type = "RKHS", wrapper = "wrap_RKHS"),
    "G_a_S@BRR"      = list(col_name = "connect_param", path = "G_a_S",
                            type = "BRR", wrapper = "wrap_BRR"),
    "G_a_S_p@BRR"    = list(col_name = "connect_param", 
                            path = c("G_a_RM", "G_a_S"),
                            type = "BRR", wrapper = "wrap_BRR_2"),
    "S_al@BRR"       = list(col_name = "connect_param", path = "G_a_S", 
                            type = "BRR", wrapper = "wrap_BRR_3"),
    "G_a_S_al@RKHS"  = list(col_name = c("connect_geno", "connect_param"), 
                            path = c("G_a_RM", "G_a_S"),
                            type = "RKHS", wrapper = "wrap_RKHS_1")
  )
  
  matrix_params <- matrix_mapping[[element]] # get needed params
  
  wrapper_function <- wrapper_functions[[matrix_params[["wrapper"]]]]
  if (element == "G_a_S_p@BRR"){
    matrix_params[["tmp_data_dir"]] <- tmp_data_dir 
    matrix_params[["dump_dir"]] <- dump_dir 
    matrix_params[["cv_id"]] <- cv_id
  }
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
  ext_parse <- "ext_dir"
  write_at <- "/proj"
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

kibreed_data_blues_within <- qread(sprintf("/proj/%s/process_cgm_data/BLUEs_within_env_cgm.qs", ext_parse))
run_data <- read_json(sprintf("/proj/%s/generate_prediction_data/%s/%s.json", ext_parse, cv_type, cv_type))

paths <- list(
  "G_a_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_a.qs", ext_parse),  # based on genetic data
  "G_d_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_d.qs", ext_parse),  # based on genetic data
  "G_aa_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_aa.qs", ext_parse),  # based on genetic data
  "ERM_l" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_linear.qs", ext_parse), # based on environment data
  "ERM_nl" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_non_linear.qs", ext_parse), # based on environment data
  "SRM" = sprintf("/proj/%s/KIBREED_data_generation/SRM.qs", ext_parse), # based on environment data
  "YRM" = sprintf("/proj/%s/KIBREED_data_generation/YRM.qs", ext_parse), # based on environment data
  "G_a_S" = sprintf("/proj/%s/process_cgm_data/g_s_mat.qs", ext_parse) # based on CGM output
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
  mutate(obs = BLUES_dt) %>%               # Create a new column 'obs' with 'BLUES_dt' values
  mutate(BLUES_dt = ifelse(idx_cv %in% test_index, NA, BLUES_dt)) %>%  # Replace 'BLUES_dt' values with NA where 'test_index' is TRUE
  filter(idx_cv %in% all_data)  %>%   # Subset rows based on 'all_data' condition
  mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
  mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
  mutate(Type = as.factor(Type), # Convert 'Type' column to a factor
         run_name = run_name) %>%
  rename(env = Env,
         geno = Geno_new,
         type = Type,
         site = Site,
         year = Year,
         blues = BLUES_dt) %>%
  select(run_name, env, type, 
         site, latlong, year, geno, 
         connect_climate, connect_param, connect_geno,
         blues, obs) %>%
  convert(fct(env, site, latlong, year, geno),
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
