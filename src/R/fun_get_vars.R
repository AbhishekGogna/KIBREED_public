# Core functions
"%!in%" <- Negate("%in%")

create_dir_file <- function(paths, file = TRUE) {
  for (path in paths) {
    if (file) {
      if (!file.exists(path)) {
        tryCatch({
          file.create(path, recursive = TRUE)
        }, error = function(e) {
          message("Error creating file '", path, "': ", e$message)
        })
      }
    } else {
      if (!dir.exists(path)) {
        tryCatch({
          dir.create(path, recursive = TRUE)
        }, error = function(e) {
          message("Error creating directory '", path, "': ", e$message)
        })
      }
    }
  }
  return(paths)
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
  
  pheno_data_mod <- pheno_data_0 %>% mutate(obs_set = ifelse(is.na(blues_var), "test", "train")) %>%
    arrange(desc(obs_set)) %>% left_join(param_data_df, by = col_name)
  
  envs <- as.character(unique(pheno_data_mod$env))
  meta_info <- NULL
  param_data_out <- NULL
  
  for (env in envs){
    ifrm(pheno_data_0_site, env = environment())
    pheno_data_0_site <- pheno_data_0 %>% filter(env == !!as.character(env))
    param_data_site <- pheno_data_mod %>% filter(env == !!as.character(env)) %>%
      mutate(across(all_of(colnames(param_data)), ~ifelse(!is.na(blues_var), ., NA), .names = "to_pred_{.col}"))
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
    mutate(obs_set = ifelse(is.na(blues_var), "test", "train")) %>%
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
  
  wrap_BRR_kin = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    kin_path <- paths[[matrix_params[["path"]]]]
    incidence_mat <- get_incidence_matrix(data = pheno_data,
                                          column_name = col_name)
    kin_mat <- qread(kin_path)
    out_0 <- get_product(incidence_mat[["mat"]], kin_mat, NULL, NULL)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS_kin = function(pheno_data, matrix_params, paths){
    col_name <- matrix_params[["col_name"]]
    kin_path <- paths[[matrix_params[["path"]]]]

    kin_mat <- qread(kin_path)
    geno_order <- pheno_data %>% pull(col_name) %>% as.vector()
    
    mat <- kin_mat[geno_order, geno_order]
    out <- list(K = mat)
  },
  
  wrap_RKHS = function(pheno_data, matrix_params, paths) {
    
    matrix_params <- matrix_params[1:2]
    
    matrix_params_1 <- sapply(matrix_params, function(x) x[[1]], USE.NAMES = TRUE, simplify = FALSE)
    matrix_params_2 <- sapply(matrix_params, function(x) x[[2]], USE.NAMES = TRUE, simplify = FALSE)
    
    matrix1 <- create_and_get_product(matrix_params_1[["col_name"]], pheno_data, paths[[matrix_params_1[["path"]]]])
    matrix2 <- create_and_get_product(matrix_params_2[["col_name"]], pheno_data, paths[[matrix_params_2[["path"]]]])
    
    out_0 <- get_had_prod(matrix1[["mat"]], matrix2[["mat"]], NULL)
    out <- list(K = out_0[["mat"]])
    return(out)
  },
  
  wrap_BRR = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params[["path"]]]])
    out <- list(X = BRR_mat[pheno_data[, col_name], ])
    return(out)
  },
  
  wrap_BRR_2 = function(pheno_data, matrix_params, paths) {
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
  
  wrap_BRR_3 = function(pheno_data, matrix_params, paths) {
    col_name <- matrix_params[["col_name"]]
    BRR_mat <- qread(paths[[matrix_params[["path"]]]])
    out_0 <- get_mat(pheno_data, col_name, BRR_mat)
    out <- list(X = out_0[["mat"]])
    return(out)
  },
  
  wrap_RKHS_1 = function(pheno_data, matrix_params, paths) {
    
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
    "G_a@RKHS"        = list(col_name = "connect_geno", path = "G_a_RM",
                            type = "RKHS", wrapper = "wrap_RKHS_kin"),
    "G_d@RKHS"        = list(col_name = "connect_geno", path = "G_d_RM",
                            type = "RKHS", wrapper = "wrap_RKHS_kin"),
    "G_aa@RKHS"       = list(col_name = "connect_geno", path = "G_aa_RM",
                            type = "RKHS", wrapper = "wrap_RKHS_kin"),
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
  
  wrapper <- matrix_params[["wrapper"]]
  wrapper_function <- wrapper_functions[[wrapper]]
  if (element == "G_a_S_p@BRR"){
    matrix_params[["tmp_data_dir"]] <- tmp_data_dir 
    matrix_params[["dump_dir"]] <- dump_dir 
    matrix_params[["cv_id"]] <- cv_id
  }
  
  if (wrapper == "wrap_inc"){
    mat_out <- wrapper_function(pheno_data, matrix_params)
  } else {
    mat_out <- wrapper_function(pheno_data, matrix_params, paths)
  }
  
  output_0 <- list(model = matrix_params[["type"]],
                   saveEffects = TRUE,
                   name = as.character(element))
  
  output <- c(output_0, mat_out)
  
  log_info(paste("ETA element for column =", element, "defined"))
  
  return(output)
}

# Tasks functions
calculate_vars <- function(paths, pheno_data, index, log_at, tmp_data){
  # define paths
  create_dir_file(paste0(tmp_data, "/tmp"), file = FALSE)
  create_dir_file(paste0(tmp_data, "/dump"), file = FALSE)
  
  tmp_data_dir <- paste0(tmp_data, "/tmp/", index, "_")
  dump_dir <- paste0(tmp_data, "/dump/", index, "_")
  
  # define pheno data
  pheno_data <- pheno_data %>% distinct(geno, .keep_all = TRUE)
  
  # define matrices
  elements <- str_split(index, "&")[[1]]
  
  ETA <- lapply(elements, return_elements,
                pheno_data = pheno_data, 
                paths = paths,
                tmp_data_dir = tmp_data_dir,
                dump_dir = dump_dir,
                cv_id = "vars_pred")
  names(ETA) <- sapply(ETA, function(x) x[["name"]])
  
  fm_base <- BGLR(y = pheno_data$blues_var,
                  ETA = ETA,
                  nIter = 15000,
                  burnIn = 2000,
                  thin = 5,
                  saveAt = dump_dir,
                  verbose = FALSE)
  
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
  
  # produce output
  return(combined_vars_sd)
}
get_vars <- function(existing_data,
                     write_at,
                     log_at,
                     tmp_at,
                     model_info,
                     key){
  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/%s.log", log_at, key)
  create_dir_file(log_at, file = FALSE)
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_at, run_instance)
  
  create_dir_file(tmp_data, file = FALSE)
  
  cat(sprintf("Logs for task: %s ------------------", key),
      file = log_file,
      sep = "\n")
  
  # Put a results dir
  create_dir_file(write_at, file = FALSE)
  
  # Produce the output

  if (grepl("acr", key)){
    
    in_paths <- existing_data[c("G_a_RM", "G_d_RM", "G_aa_RM")]
    names(in_paths)[grepl("acr", names(in_paths))] <- "pheno_data"
    
    pheno_data <- qread(existing_data[["pheno_acr"]]) %>%
      mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
      mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
      mutate(Type = as.factor(Type), # Convert 'Type' column to a factor
             run_name = run_name) %>%
      rename(series = Series,
             geno = Geno_new,
             type = Type,
             blues_var = BLUEs,
             connect_geno = connect_geno_data) %>%
      select(run_name, series, type, 
             geno, 
             connect_geno,
             blues_var) %>%
      convert(fct(series, geno),
              chr(type, run_name, starts_with("connect")),
              num(blues_var)) 
  } else if (grepl("wtn", key)){
    
    in_paths <- existing_data[c("G_a_RM", "G_d_RM", "G_aa_RM", "ERM_l", "ERM_nl", "SRM", "YRM", "G_a_S")]
    names(in_paths)[grepl("wtn", names(in_paths))] <- "pheno_data"
    
    pheno_data <- qread(existing_data[["pheno_wtn"]]) %>%
      mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
      mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
      mutate(Type = as.factor(Type), # Convert 'Type' column to a factor
             run_name = run_name) %>%
      rename(env = Env,
             geno = Geno_new,
             type = Type,
             site = Site,
             year = Year,
             blues_var = BLUES_dt) %>%
      select(run_name, env, type, 
             site, latlong, year, geno, 
             connect_climate, connect_param, connect_geno,
             blues_var) %>%
      convert(fct(env, site, latlong, year, geno),
              chr(type, run_name, starts_with("connect")),
              num(blues_var)) 
  }
  
  # Set blas threads
  RhpcBLASctl::blas_set_num_threads(5)
  
  cl <- makeCluster(nrow(model_info), type = "FORK")
  
  registerDoParallel(cl)
  
  data_var <- foreach(idx_str = model_info$model_specs,
                      .packages = c("tidyverse", 
                                    "hablar", 
                                    "BGLR",
                                    "AGHmatrix",
                                    "qs")) %dopar%
    calculate_vars(paths = in_paths,
                   pheno_data = pheno_data,
                   index = idx_str,
                   log_at = log_at,
                   tmp_data = tmp_data)
  
  stopImplicitCluster()
  
  # Reset blas thread count
  RhpcBLASctl::blas_set_num_threads(1)
  
  data_var_df_0 <- as.data.frame(do.call(bind_rows, data_var), row.names = NULL)
  data_var_df <- cbind(data_var_df_0, "model" = model_info$models)
  
  if(!is.null(data_var_df)){
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot_df <- data_var_df %>% 
      pivot_longer(cols = all_of(colnames(data_var_df_0)), names_to = "type", values_to = "value") %>% 
      filter(grepl("var", type, ignore.case = T), !is.na(value)) %>% 
      rowwise() %>% 
      mutate(comp = str_split(type, "@")[[1]][1], prior = str_split(type, "@")[[1]][2]) %>%
      ungroup() %>%
      group_by(model) %>%
      mutate(p_cent = value/sum(value)) %>%
      ungroup()
    
    color_df_0 <- plot_df %>% distinct(comp)
    color_df <- color_df_0 %>% bind_cols("color" = col_vector[1:nrow(color_df_0)])
    
    plot_df_col <- plot_df %>%
      left_join(color_df, by = "comp")
    
    if(grepl("wtn", key)){
      plot_df_col$comp <- factor(plot_df_col$comp,levels = c("E_i", "G_i", "G_a", "G_d", "G_aa",
                                                             "ERM_l", "G_a_ERM_l", "ERM_nl",
                                                             "G_a_ERM_nl", "S", "Y", "G_a_S_i",
                                                             "G_a_S", "G_a_Y", "err_var"))
      
    } else {
      plot_df_col$comp <-factor(plot_df_col$comp,levels = c("G_a", "G_d", "G_aa",
                                                            "err_var"))
    }

    var_plot <- plot_df_col %>%
      ggplot(aes(x = model, y = p_cent, fill = comp)) +
      geom_bar(stat = "identity", position = "fill", colour = "black") +
      theme_classic() +
      scale_fill_manual(values = color_df %>% deframe())
    
      # geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5)) + # if labels are desired
      # facet_grid(run_type ~ cv_idx, scales = "free_x")
    
    if(grepl("acr", key)){
      wd <- 8.4
      ht <- 8.4
    } else if (grepl("wtn", key)){
      wd <- 16.8
      ht <- 12.6
    }

    ggsave(plot = var_plot, filename = sprintf("%s/%s_var_plot.png", write_at, key), 
           width = wd, height = ht, units = "cm", dpi = 600)
  } else {
    var_plot <- NULL
  }
  
  
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- write_at
  out[[sprintf("vars_%s", key)]] <- data_var_df
  out[[sprintf("plot_vars_%s", key)]] <- var_plot
  
  return(out)
  
}
make_acr_plots <- function(linked_proj,
                           data) {
  # Sequester data
  log_file <- data[["log_file"]]
  write_at <- data[["write_at"]]
  var_plot <- data[["plot_vars_pheno_data_acr"]][["data"]] %>%
    filter(model %in% c("M_2", "M_3")) %>%
    convert(chr(comp)) %>%
    mutate(comp = ifelse(comp == "err_var", "err", comp)) %>%
    convert(fct(comp, .args = list(levels = c("G_a", "G_d", "G_aa", "err")))) %>%
    ggplot(aes(x = model, y = p_cent, fill = comp)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    theme_classic() +
    scale_x_discrete(labels=c("M_2" = "GBLUP_D", "M_3" = "E-GBLUP_D")) + 
    ylab("Proportion explained") +
    xlab("Model") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 11.25, vjust = 0.5), 
          legend.position = "bottom") +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  pred_data <- tar_read(pred_data, store = sprintf("/proj/store/%s", linked_proj))
  acr_pred_data <- pred_data[grep("acr", names(pred_data))]
  
  # Produce plots
  p1_data <- acr_pred_data$plot_acr_cv$data %>%
    filter(name %in% c("GBLUP_D", "E-GBLUP_D", "acr_CNN")) %>%
    mutate(type = as.factor(ifelse(type == "Non_hybrid", "Lines", "Hybrids")))
  p1_data_2 <- acr_pred_data$plot_acr_cv$layers[[3]]$data %>%
    filter(name %in% c("GBLUP_D", "E-GBLUP_D", "acr_CNN")) %>%
    mutate(type = as.factor(ifelse(type == "Non_hybrid", "Lines", "Hybrids")))
  p1 <- ggplot(p1_data,
               aes(x = name, y = value, color = type)) +
    geom_jitter(size = 0.75,  position = position_jitterdodge(), alpha = 0.5) +
    coord_cartesian(ylim = c(0, 1))+
    geom_hline(yintercept  = 0.75, linetype ="dashed", color = "red") + 
    geom_text(aes(label = value, y = 0.1), data = p1_data_2, 
              angle = 90, size = 3, position = position_dodge(1), show.legend = FALSE)+
    labs(x = "Model", y = "Prediction ability", color = "Type")+
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 11.25, vjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = c("grey", "white"))
  
  p2_data <- acr_pred_data$plot_acr_sce_1$data %>%
    filter(name %in% c("E-GBLUP_D", "acr_CNN")) %>%
    mutate(type = as.factor(ifelse(type == "Non_hybrid", "Lines", "Hybrids")))
  
  p2 <- ggplot(p2_data,
               aes(x = train_val, y = value, color = type)) +
    geom_point(size = 0.1) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    facet_grid(name~type) +
    ylab("Prediction ability") +
    xlab("Training set size") +
    geom_smooth(method = "loess") +
    stat_cor(aes(label =paste(after_stat(r.label), cut(after_stat(p), 
                                                breaks = c(-Inf, 0.05, Inf),
                                                labels = c("'*'", "''")), 
                              sep = "~")),
             p.accuracy = 0.01, r.accuracy = 0.01, 
             label.x = 0, label.y = 0.9, size = 3) +
    theme_classic(base_size = 10) +
    scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme(legend.position = "none") +
    scale_fill_manual(values = c("grey", "white"))
  
  jp_1 <- ggarrange(var_plot, p1,
                    labels = c("", ""),
                    ncol = 1, 
                    common.legend = FALSE,
                    font.label = list(size = 10, color = "black", face = "plain", family = NULL))
  
  joint_plot_1 <- ggarrange(jp_1, p2, ncol = 2, labels = c("", ""), widths = c(1, 2),
                            font.label = list(size = 10, color = "black", face = "plain", family = NULL))
  
  ggsave(joint_plot_1, filename = sprintf("%s/%s", write_at, "pred_acr_res.png"),
         width = 17, height = 16, units = "cm", dpi = 600, bg = "white")
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- write_at
  out[["joint_plot_1"]] <- joint_plot_1
  
  return(out)
}
