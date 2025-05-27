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

cmd_RD_func<-function(SNP, print_status = T){
  ## calculate the matrix of Roger's distance
  ## SNP must as 0,1,2 
  ## Transfer element class from integer/character to numeric
  ## Data may have missing values, use the na.rm function
  nmarker <- ncol(SNP)
  No.Geno<-dim(SNP)[1]             # Auto input number of genotypes
  Geno<-rownames(SNP)
  RD <- matrix(0,No.Geno,No.Geno)
  if(print_status == T){
    print(paste0("Start at = ", Sys.time()))
  }
  system.time(
    for(i in 1:No.Geno)
    {
      if(is.element("2",names(table(SNP[,1:10])))=="FALSE"){
        print("SNP marker may not fit 0,1,2")
        break
      }
      if(length(SNP[i,])!=nmarker){
        print("input matrix is wrong")
        break}
      for (j in 1:No.Geno)
      {
        RD[i,j] <- mean(abs(SNP[i,]-SNP[j,])/2,na.rm=TRUE)
        if(length(SNP[i,])!=nmarker){
          print("input matrix is wrong")
          break}
      }
      if(print_status == T){
        if(i%%1000==0){print(paste("row",i,"finished",Sys.time()))}
      }
    }
  )
  if(print_status == T){
    print(paste0("End at = ", Sys.time()))
  }
  rownames(RD)<-Geno
  colnames(RD)<-Geno
  return(RD)
}

cmd_RD_func_mod <- function(SNP, print_status = TRUE) {
  # Calculate the matrix of Roger's distance
  # SNP must be 0, 1, 2 
  # Transfer element class from integer/character to numeric
  # Data may have missing values, use the na.rm function
  nmarker <- ncol(SNP)
  No.Geno <- nrow(SNP)  # Auto input number of genotypes
  Geno <- rownames(SNP)
  RD <- matrix(0, No.Geno, No.Geno)
  
  if (print_status) {
    print(paste0("Start at = ", Sys.time()))
  }
  
  system.time({
    for (i in 1:No.Geno) {
      if (is.element("2", names(table(SNP[, 1:10]))) == FALSE) {
        print("SNP marker may not fit 0, 1, 2")
        break
      }
      if (length(SNP[i, ]) != nmarker) {
        print("input matrix is wrong")
        break
      }
      for (j in i:No.Geno) {
        RD[i, j] <- mean(abs(SNP[i, ] - SNP[j, ])/2, na.rm = TRUE)
        if (length(SNP[i, ]) != nmarker) {
          print("input matrix is wrong")
          break
        }
      }
      if (print_status && i %% 1000 == 0) {
        print(paste("row", i, "finished", Sys.time()))
      }
    }
  })
  
  # Mirror the lower half to the upper half
  for (col in 1:No.Geno) {
    is_miss <- which(RD[, col] == 0)
    
    for (row in is_miss[-1]) {
      RD[row, col] <- RD[col, row]
    }
  }
  
  if (print_status) {
    print(paste0("End at = ", Sys.time()))
  }
  
  rownames(RD) <- Geno
  colnames(RD) <- Geno
  return(RD)
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
  
  print(paste("ETA element for column =", element, "defined"))
  
  return(output)
}

calculate_pred <- function(paths, p_data_idx, p_data, index, log_at, tmp_data){
  
  # Specify data
  pheno_data <- p_data[[p_data_idx]]
  tmp_data <- paste0(tmp_data, "/", p_data_idx)
  
  # define paths
  create_dir_file(paste0(tmp_data, "/tmp"), file = FALSE)
  create_dir_file(paste0(tmp_data, "/dump"), file = FALSE)
  
  tmp_data_dir <- paste0(tmp_data, "/tmp/", index, "_")
  dump_dir <- paste0(tmp_data, "/dump/", index, "_")
  
  # define matrices
  elements <- str_split(index, "&")[[1]]
  
  ETA <- lapply(elements, return_elements,
                pheno_data = pheno_data, 
                paths = paths,
                tmp_data_dir = tmp_data_dir,
                dump_dir = dump_dir,
                cv_id = "vars_pred")
  names(ETA) <- sapply(ETA, function(x) x[["name"]])
  
  fm_base <- BGLR(y = pheno_data$blues,
                  ETA = ETA,
                  nIter = 15000,
                  burnIn = 2000,
                  thin = 5,
                  saveAt = dump_dir,
                  verbose = FALSE)
  
  pheno_data$pred <- as.numeric(fm_base$yHat)
  
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
  
  # Write a section in the log file
  cat("Done with predicting data",
      file = log_at,
      sep = "\n",
      append = TRUE)
  
  # produce output
  out <- list()
  out[["vars"]] <- combined_vars_sd 
  out[["pred"]] <- pheno_data
  return(out)
}

create_cluster_plot <- function(mds.df, wss_values, label, write_at, n_clust, labs, col_vector) {
  kmclusters <- wss_values[[n_clust]][["clusters"]]
  mds.df$groups <- kmclusters
  mds.df$env <- rownames(mds.df)
  
  clust_plot <- ggscatter(mds.df,
                          x = "V1", 
                          y = "V2", 
                          xlab = sprintf("PC1 = %s%%", as.character(round(labs[1], 1))),
                          ylab = sprintf("PC2 = %s%%", as.character(round(labs[2], 1))),
                          label = label,
                          color = "groups", 
                          #palette = col_vector[1:n_clust], 
                          size = 1, 
                          ellipse = TRUE, 
                          ellipse.type = "convex", 
                          repel = TRUE,
                          ggtheme = theme_gray(),
                          title = sprintf("Clusters = %s", n_clust))
  
  ggsave(plot = clust_plot, filename = sprintf("%s/clust_plot_%s_%s.png", write_at, n_clust, ifelse(is.null(label), "no_label", label)), 
         width = 16.8, height = 16.8, units = "cm", dpi = 600)
  return(clust_plot)
}

record_raster_grob <- function(width, height, dpi = 300, expr){
  t <- tempfile(fileext = ".png")
  on.exit(unlink(t), add = TRUE)
  png(t, width = width, height = height, res = dpi, units = "cm")
  tryCatch( expr , finally = dev.off() )
  grid::rasterGrob(png::readPNG(t))
}

standardize_pair <- function(pair) {
  parts <- strsplit(pair, "&")[[1]]
  return(paste(sort(parts), collapse = "&"))
}

plot_series <- function(series, data_list) {
  # Open a new plotting device
  par(mfrow = c(length(series), length(series)), 
      mar = c(2,2,2,1), xpd = NA, cex = 0.5, font = 1) # b,l,t,r

  # Loop through each element in the list and plot histograms
  for (i in series) {
    for (j in series){
      #print(id)
      if(i == j | j > i){
        id <- sprintf("Exp_%s&Exp_%s", i, j)
        data <- data_list[[id]]
        data_mean <- round(mean(data), 2)
        if(i == j){
          plot_axes <- TRUE
        } else {
          plot_axes <- FALSE
        }
        #plot_axes <- TRUE
        h <- hist(data, plot = FALSE, breaks = 25)
        h$density = h$counts/sum(h$counts)*100
        max = as.integer(max(h$density))
        plot(h,freq=FALSE, main = "", 
             xlim = c(0, 0.5), ylim = c(0, 40), 
             xlab = "", ylab = "", axes = plot_axes,
             col = "lightblue", border = "black")
        title(main = sprintf("%s\nh%s m%s", id, max, data_mean),
              adj = 0.3, line = -1)
      } else {
        id_2 <- sprintf("Exp_%s&Exp_%s", i, j)
        data_2 <- data_list[[id_2]]
        data_2_mean <- round(mean(data_2), 2)
        wtn_series <- sprintf("Exp_%s&Exp_%s", j, j)
        data_series <- data_list[[wtn_series]]
        data_series_mean <- round(mean(data_series), 2)
        pcent_diff <- round(100*((data_2_mean - data_series_mean)/data_series_mean), 2)
        plot(x = 0:10, y = 0:10, ann = F,bty = "n",type = "n",
             xaxt = "n", yaxt = "n")
        text(x = 5,y = 5, 
             as.character(pcent_diff),
             cex = 3)
      }
    }
  }
  
  # reset the plotting device
  par(mfrow = c(1, 1))
}

# Task functions
get_RD_mat <- function(existing_data,
                       write_at,
                       log_at,
                       tmp_at,
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
  
  # Load data
  g_data <- qread(existing_data[["g_data"]])
  
  # Write a section in the log file
  cat("GN data loaded",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  save_at <- sprintf("%s/%s/RD_data.qs", core_paths[["results_R"]], run_name)
  if(!file.exists(save_at)){
    RD <- cmd_RD_func_mod(g_data)
    qsave(RD, save_at)
  } else {
    RD <- qread(save_at)
  }
  
  # Write a section in the log file
  cat("RD data generated",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Check availability
  check <- sum(rownames(RD) == row.names(g_data)) == nrow(g_data)
  
  if(check) {
    RD_out <- RD
  } else {
    RD_out <- NULL
  }

  # return output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["RD_out"]] <- RD_out
  return(out)
}

get_core_set <- function(existing_data, data) {
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  RD_out <- data[["RD_out"]]
  p_wtn <- qread(existing_data[["p_wtn"]]) %>%
    rename(connect_climate = Connect_at,
           connect_geno = connect_geno_data) %>%
    filter(Type != "Hybrid",
           BLUES_dt > 0)
  
  # Write a section in the log file
  cat("Starting with core sampling",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output: Subset RD_mat for genotypes of interest
  geno <- p_wtn %>% distinct(Env, connect_geno) %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = "connect_geno", names_from = "Env", values_from = "present") %>%
    rowwise() %>%
    mutate(freq = sum(c_across(where(is.numeric)), na.rm = TRUE))
  
  overview <- table(cut(geno$freq, breaks = c(0, 2, 4, 8, 16, 32, 64, 128)))
  
  RD_sub <- RD_out[geno$connect_geno, geno$connect_geno]
  
  core_obj <- list(
    objective("EN", "PD", weight = 1)
    #, objective("AN", "GD", weight = 1)
  ) # multiple objectives (custom weight)
  
  set.seed(12345)
  core <- sampleCore(data = distances(RD_sub),
                     obj = core_obj, # default. can include more than one, in which case normalization is done for combined results
                     size = 500,
                     #always.selected = integer(0), #geno vector to be always be selected in the core collection
                     #never.selected = integer(0), #geno vector to be never be selected in the core collection
                     mode = "default",
                     normalize = TRUE, # For single-objective configurations, this argument is ignored.
                     time = NA,
                     impr.time = NA, # Defaults to ten or two seconds, mode = default or fast
                     steps = NA,
                     impr.steps = NA,
                     indices = FALSE,
                     verbose = FALSE)
  
  # Write a section in the log file
  cat("Done with core sampling",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Extra table
  ## identify best and worst performers per environment cluster?
  core_dist <- geno %>% filter(connect_geno %in% core$sel) %>% 
    relocate(freq, .after= "connect_geno")
  
  interval_size <- 3
  
  breaks_in <- seq(0, max(core_dist$freq) + interval_size, interval_size)
  core_dist_series_0 <- core_dist %>% 
    pivot_longer(!c("connect_geno", "freq"), names_to = "env") %>% 
    filter(!is.na(value), env != "freq") %>% 
    mutate(freq_env = freq,
           freq_env_class_integer = cut(freq_env, breaks = breaks_in, include.lowest = TRUE, labels = FALSE),
           freq_env_class = paste0("(", (as.numeric(freq_env_class_integer) - 1) * interval_size + 1, ", ", as.numeric(freq_env_class_integer) * interval_size, "]"), #upperbound of lower class +1 to upper bound of this class
           series = gsub("^([[:alnum:]_]+_\\d_[[:alnum:]]+).*", "\\1", env),
           value = 1) 
  multi_series <- core_dist_series_0 %>%
    select(connect_geno, freq_env_class, series) %>%
    distinct() %>% count(connect_geno) %>% filter(n > 1)
  core_dist_series <- core_dist_series_0 %>%
    select(connect_geno, freq_env_class, series) %>% 
    mutate(series = ifelse(connect_geno %in% multi_series$connect_geno, "Exp_multi", series)) %>%
    distinct() %>%
    group_by(freq_env_class, series) %>% 
    summarize(freq_ser = n(), .groups = "drop") %>% 
    pivot_wider(id_cols = freq_env_class, names_from = series, values_from = freq_ser) %>%
    rowwise() %>%
    mutate(class_lower_bound = as.numeric(gsub("\\((\\d+)\\,.*", "\\1", freq_env_class)),
           total = sum(c_across(starts_with("Exp")), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(class_lower_bound) %>%
    select(freq_env_class, Exp_1_Hywheat , Exp_2_Zucht, Exp_3_Zucht, Exp_4_Zucht, Exp_5_Zucht, Exp_6_KWS, Exp_7_GABI, Exp_multi, total)
  
  # Extra figure
  ## how many of the core set are in the extremes
  core_ext <- p_wtn %>%
    mutate(present_in = ifelse(connect_geno %in% core$sel, "core", "general"))
  core_ext_plot <- ggplot(aes(x = BLUES_dt), data = core_ext %>% filter(present_in == "general")) +
    geom_histogram(bins = 50, fill = "#00AFBB", color = "black") +
    geom_rug(aes(x = BLUES_dt), data = core_ext %>% filter(present_in == "core"), color = "#E7B800") +
    labs(x= "yield (quintal per ha)", y = "frequency")+
    theme_classic()
    
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["core"]] <- core
  out[["geno_dist"]] <- geno
  out[["geno_dist_intervals"]] <- as.data.frame(overview)
  out[["core_dist_series"]] <- core_dist_series
  out[["core_ext_plot"]] <- core_ext_plot
  
  return(out)
}

get_training_data <- function(existing_data, data) {
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  core_set <- data[["core"]]
  p_wtn <- qread(existing_data[["p_wtn"]]) %>%
    mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
    mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
    mutate(Type = as.factor(Type)) %>% # Convert 'Type' column to a factor
    rename(connect_climate = Connect_at,
           connect_geno = connect_geno_data,
           env = Env,
           geno = Geno_new,
           type = Type,
           site = Site,
           year = Year,
           blues = BLUES_dt) %>%
    select(env, type, 
           site, latlong, year, geno, 
           connect_climate, connect_geno,
           blues) %>%
    convert(fct(env, site, latlong, year, geno),
            chr(type, starts_with("connect")),
            num(blues)) %>%
    filter(!is.na(blues),
           type != "Hybrid")
  
  # Write a section in the log file
  cat("Starting with defining a training data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  core <- p_wtn %>% 
    filter(connect_geno %in% core_set$sel) %>%
    mutate(set = "core")
  rest <- p_wtn %>% 
    filter(connect_geno %!in% core_set$sel) %>%
    mutate(set = "rest")
  
  geno_connection <- p_wtn %>%
    distinct(geno, connect_geno, type)

  core_reformatted <- core %>%
    select(-connect_geno) %>%
    group_by(env, site, latlong, year, connect_climate, geno) %>%
    summarise(blues = mean(blues), .groups = "drop") %>%
    pivot_wider(id_cols = c("env",
                            "site", "latlong", 
                            "year", "connect_climate"),
                names_from = "geno",
                values_from = "blues") %>%
    pivot_longer(!(c("env",
                     "site", "latlong", 
                     "year", "connect_climate")),
                 names_to = "geno",
                 values_to = "blues") %>%
    left_join(geno_connection, by = "geno") %>%
    select(all_of(colnames(p_wtn))) %>%
    mutate(set = "core")
  
  training_data <- rest %>% bind_rows(core_reformatted)
  #training_data <- core_reformatted
  
  miss_prop <- round(sum(is.na(training_data$blues))/nrow(training_data), 2)
  print(sprintf("Miss prop = %s", miss_prop))
  cat(sprintf("Miss prop = %s", miss_prop),
      file = log_file,
      sep = "\n",
      append = TRUE)
  # Write a section in the log file
  cat("Done with defining a training data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["training_data"]] <- training_data
  return(out)
}

get_training_data_2 <- function(existing_data, data) {
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  #core_set <- data[["core"]]
  p_wtn <- qread(existing_data[["p_wtn"]]) %>%
    mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
    mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
    mutate(Type = as.factor(Type)) %>% # Convert 'Type' column to a factor
    rename(connect_climate = Connect_at,
           connect_geno = connect_geno_data,
           env = Env,
           geno = Geno_new,
           type = Type,
           site = Site,
           year = Year,
           blues = BLUES_dt) %>%
    select(env, type, 
           site, latlong, year, geno, 
           connect_climate, connect_geno,
           blues) %>%
    convert(fct(env, site, latlong, year, geno),
            chr(type, starts_with("connect")),
            num(blues)) %>%
    filter(!is.na(blues),
           type != "Hybrid")
  
  # Write a section in the log file
  cat("Starting with defining a training data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Subset the extreme lines based on BLUEs
  BLUEs <- qread("/proj/results/R/KIBREED_data_generation/BLUES_acr_env.qs")
  select_from <- BLUEs %>% distinct(unique_idx, .keep_all = TRUE) %>%
    filter(Type == "Line") 
  needed <- as.integer(nrow(select_from) * 0.1)
  extreme_n <- as.integer(needed * 0.7) 
  random_n <- as.integer(needed * 0.3) 
  
  extremes <- select_from %>%
    arrange(desc(BLUEs)) %>%
    top_n(extreme_n, BLUEs) %>%
    slice(1:extreme_n) %>%
    #top_frac(0.07, BLUEs) %>%
    pull(connect_geno_data) %>%
    as.vector()
  
  set.seed(12345)
  randoms <- select_from %>%
    filter(connect_geno_data %!in% extremes) %>%
    sample_n(random_n) %>%
    pull(connect_geno_data) %>%
    as.vector()
  
  selected <- unique(c(extremes, randoms))
  
  # Produce output
  pred_set <- p_wtn %>% 
    filter(connect_geno %in% selected) %>%
    mutate(set = "pred_set")
  rest <- p_wtn %>% 
    filter(connect_geno %!in% selected) %>%
    mutate(set = "rest",
           set_type = "rest")
  
  geno_connection <- p_wtn %>%
    distinct(geno, connect_geno, type)
  
  pred_set_reformatted <- pred_set %>%
    select(-connect_geno) %>%
    group_by(env, site, latlong, year, connect_climate, geno) %>%
    summarise(blues = mean(blues), .groups = "drop") %>%
    pivot_wider(id_cols = c("env",
                            "site", "latlong", 
                            "year", "connect_climate"),
                names_from = "geno",
                values_from = "blues") %>%
    pivot_longer(!(c("env",
                     "site", "latlong", 
                     "year", "connect_climate")),
                 names_to = "geno",
                 values_to = "blues") %>%
    left_join(geno_connection, by = "geno") %>%
    select(all_of(colnames(p_wtn))) %>%
    mutate(set = "pred_set",
           set_type = ifelse(connect_geno %in% extremes, "extremes",
                             ifelse(connect_geno %in% randoms, "randoms", "rest")))
  
  training_data <- rest %>% bind_rows(pred_set_reformatted)
  #training_data <- core_reformatted
  
  miss_prop <- round(sum(is.na(training_data$blues))/nrow(training_data), 2)
  print(sprintf("Miss prop = %s", miss_prop))
  cat(sprintf("Miss prop = %s", miss_prop),
      file = log_file,
      sep = "\n",
      append = TRUE)
  # Write a section in the log file
  cat("Done with defining a training data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["extremes"]] <- extremes
  out[["randoms"]] <- randoms
  out[["training_data"]] <- training_data
  return(out)
}

predict_missing <- function(existing_data, data_1, data_2, model, debug){
  # Sequester data
  log_file <- data_1[["log_file"]]
  tmp_data <- data_1[["tmp_data"]]
  write_at <- data_1[["write_at"]]
  training_data <- list("core" = data_1[["training_data"]], 
                        "max_yield" = data_2[["training_data"]])
  
  if (debug) {
    training_data <- lapply(training_data, function(x) {
      all_geno <- sample(x$geno, 10)
      out <- x %>%
        filter(geno %in% all_geno)
    })
  }
  
  # Write a section in the log file
  cat("Start with predicting data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  RhpcBLASctl::blas_set_num_threads(10)
  n_cases <- length(training_data)
  c1 <- makeCluster(n_cases, type = "FORK", outfile = sprintf("%s/%s", tmp_data, "cluster.log"))
  registerDoParallel(c1)
  system.time({
    out_raw <- foreach(parts = c("core", "max_yield"),
                       .packages = "BGLR",
                       .export = c("%!in%")) %dopar% {
                         predicted_data <- calculate_pred(paths = existing_data,
                                                          p_data_idx = parts, 
                                                          p_data = training_data,
                                                          index = model,
                                                          log_at = sprintf("%s/predict_missing.log", tmp_data),
                                                          tmp_data = tmp_data)
                         return(predicted_data)
                       }
  })
  stopCluster(c1)
  
  core_data <- out_raw[[1]]
  max_yield <- out_raw[[2]]
  
  train_pred <- core_data$pred %>% filter(set == "core")

  # Write a section in the log file
  cat(sprintf("Done with predicting data. obs/pred = %s/%s. train acc = %s", 
              train_pred %>% filter(!is.na(blues)) %>% nrow(),
              train_pred %>% filter(is.na(blues)) %>% nrow(),
              train_pred %>% filter(!is.na(blues)) %>% 
                summarize(corr = cor(blues, pred)) %>%
                pull(corr) %>%
                as.numeric()),
      file = log_file,
      sep = "\n",
      append = TRUE)

  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["predicted_data_core"]] <- core_data
  out[["predicted_max_yield"]] <- max_yield

  return(out)
}

get_residuals <- function(existing_data, data){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  predicted_data <- data[["predicted_data_core"]]$pred %>%
    convert(fct(geno)) %>% 
    filter(set == "core") %>%
    select(-set) # modify this part in training set creation where set == "core" is lost
  
  # Write a section in the log file
  cat("Start with deriving residuals",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  asr_fit <- asreml(fixed = pred ~ 1,
                    random = ~ geno + env,
                    data = predicted_data)
  
  asr_fir_sum <- summary(asr_fit)$varcomp
  
  fit <- asr_fir_sum["geno", "component"]/(asr_fir_sum["geno", "component"]+asr_fir_sum["units!R", "component"])
  
  predicted_data$resid <- asr_fit$residuals
  
  # Write a section in the log file
  cat(sprintf("Done with deriving residuals. h2 = %s", fit),
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["residual_vals"]] <- predicted_data
  return(out)
}

get_env_clusters <- function(existing_data, data, ec_at){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  residual_vals <- data[["residual_vals"]]
  ec_mat <- tar_read(KIBREED_data_full, store = ec_at)[["ec_mat"]]
  ec_mat_scaled <- tar_read(overviews, store = ec_at)[["ec_mat_scaled"]]
  
  # Write a section in the log file
  cat("Data loaded",
      file = log_file,
      sep = "\n",
      append = TRUE)

  # Reformat to exg ------------------------------------------------------------
  residual_vals_wide <- residual_vals %>%
    distinct(env, connect_climate, geno, year, resid) %>%
    pivot_wider(id_cols = c("env", "connect_climate", "year"), names_from = "geno", values_from = "resid") %>%
    select(-env) %>%
    distinct(connect_climate, .keep_all = T) %>%
    mutate(Site = gsub("^\\S+\\_\\S+\\_(\\S+)\\_\\d+", "\\1", connect_climate),
           Year = gsub("^\\S+\\_\\S+\\_\\S+\\_(\\d+)", "\\1", connect_climate),
           abb_1 = gsub("(\\S{4}).*", "\\1", Site), 
           abb_2 = gsub("\\d{2}(\\d+)", "\\1", Year),
           abb = paste0(abb_1, abb_2)) %>%
    group_by(abb) %>%
    mutate(ids_mod = ifelse(row_number() != 1, paste0(abb, "_", row_number()), abb)) %>%
    ungroup() %>%
    select(-Site, -Year, -abb_1, -abb_2, -abb) %>%
    relocate(ids_mod)
  
  colors_years <- residual_vals_wide %>% distinct(year)
  n_years <- nrow(colors_years)
  #coul <- brewer.pal(n = 12, name = "Paired")[1:n_years]
  coul <- c("#6EB5FF", "#1E5DAB", "#7FFF7F", "#2E8B57", "#FF6347", "#CD5C5C", "#BA55D3", "#800080", "#FFD700", "#DAA520", "#FF6961")
  
  colors_years$color <- coul
  
  residual_vals_wide <- residual_vals_wide %>%
    left_join(colors_years, by = "year") %>%
    relocate(color, .after = "year")
  
  cat("Reformatted residuals to rectangular matrix",
      file = log_file,
      sep = "\n",
      append = TRUE)
  # Derive relationship matrix -------------------------------------------------
  mat <- as.matrix(residual_vals_wide[, 5:ncol(residual_vals_wide)])
  row.names(mat) <- residual_vals_wide$ids_mod
  #mat_scaled <- scale(mat, center = T, scale = T) # not done since all data points are residual values
  
  cat("Relationship matrix derived",
      file = log_file,
      sep = "\n",
      append = TRUE)
  # perform hierarchical clustering and produce plots
  #tree_data <- mat_scaled %>%
  #  dist(method = "euclidean") %>% 
  #  hclust() %>% as.dendrogram() %>% 
  #  set("branches_lwd", c(3,2,3)) %>% 
  #  set("branches_lty", c(1,1,3,1,1,2)) %>% 
  #  set("labels_cex", c(.9,1.2)) 
  #label_order = data.frame("ids_mod" = tree_data %>% labels) 
  #colors_to_use <- label_order %>% 
  #  left_join(residual_vals_wide %>% distinct(ids_mod, color), 
  #            by = "ids_mod") %>%
  #  pull(color) %>% as.vector()
  #labels_colors(tree_data) <- colors_to_use
  #
  #pdf(sprintf("%s/%s", write_at, "env_geno_resid_phylo_plot.pdf"), width = 9.5, height = 9.5)
  #par(mar = c(2,2,2,2), bg=NA)
  #(cic_plot <- circlize_dendrogram(tree_data))
  #dev.off()
  
  # Derive optimal clusters  ---------------------------------------------------
  options("ggrepel.max.overlaps" = Inf)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) # different from those used for hierarchical clustering
  
  dist.object <- dist(mat, method="euclidean", diag=T, upper=T)
  dist.matrix <- as.matrix(dist.object)
  mds <- cmdscale(dist.matrix,eig=TRUE, k=2, add = T)
  mds.df <- as.data.frame(mds$points)
  labs <- 100*(mds$eig/sum(mds$eig))[1:2]
  
  k.values <- 1:50
  wss <- function(df, k) {
    k_fit <- kmeans(df, k, nstart = 10 )
    
    # silhouette score
    if(k > 1){
      s <- silhouette(k_fit$cluster, dist(df))
      ss_mean <- mean(s[, 3])
      ss_median <- quantile(s[, 3], probs = 0.5, names = FALSE)
      ss_quantiles <- quantile(s[, 3], probs = c(0.25, 0.75), names = FALSE)
    } else {
      ss_mean <- NA
      ss_median <- NA
      ss_quantiles <- c(NA, NA)
    }
    
    out <- list()
    out[["wss"]] <- k_fit$tot.withinss
    out[["n_clust"]] <- k
    out[["min_size"]] <- min(k_fit$size) 
    out[["clusters"]] <- as.factor(k_fit$cluster)
    out[["ss"]] <- ss_median
    out[["ss_q25"]] <- ss_quantiles[1]
    out[["ss_q75"]] <- ss_quantiles[2]
    return(out)
  }
  set.seed(2345)
  wss_values <- lapply(k.values, wss, df = mds.df)
  names(wss_values) <- k.values
  
  s_vals <- sapply(wss_values, function(x) x[["ss"]])
  s_q25 <- sapply(wss_values, function(x) x[["ss_q25"]])
  s_q75 <- sapply(wss_values, function(x) x[["ss_q75"]])
  max_ss <- sort(s_vals[!is.na(s_vals)], decreasing = T)[1]
  n_ss <- as.numeric(names(max_ss))
    
  n_clust_qual <- sapply(wss_values, function(x) x[["min_size"]])
  n_clust_0 <- n_clust_qual[which(n_clust_qual >= 3)]
  n_clust <- as.numeric(names(n_clust_0[length(n_clust_0)]))
  wss_vals <- as.numeric(sapply(wss_values, function(x) x[["wss"]]))
  # Initialize a variable to store the kink
  kink <- 1
  
  # Loop through the vector
  for (i in 2:length(wss_vals)) {
    # Check if the current value is greater than the previous value
    if (wss_vals[i] > wss_vals[i - 1]) {
      kink <- i
      break  # Break the loop when the condition is met
    }
  }
  
  kink_fin <- kink - 1
  
  kmean_n_clust <- data.frame("clusters" = k.values,
                              "wss" = wss_vals) %>%
    ggplot(aes(x = clusters, y = wss)) +
    geom_line(size = 0.25) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = n_clust, linetype="dotted", size = 0.25, color = "blue") +
    geom_vline(xintercept = kink_fin, linetype="dotted", size = 0.25, color = "green") +
    theme_classic()
  
  ss_plot <- data.frame("clusters" = k.values,
                        "ss" = s_vals,
                        "q_25" = s_q25,
                        "q_75" = s_q75) %>%
    ggplot(aes(x = clusters, y = ss)) +
    geom_errorbar(aes(ymin = q_25, ymax = q_75), size = 0.25) +
    geom_line(linewidth = 0.25) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = n_clust, linetype="dotted", size = 0.25, color = "blue") +
    geom_vline(xintercept = kink_fin, linetype="dotted", size = 0.25, color = "green") +
    geom_vline(xintercept = n_ss, linetype="dotted", size = 0.25, color = "red") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic()
  
  joint_ss_elbow <- ggarrange(kmean_n_clust, ss_plot,
                              labels = c("a", "b"),
                              font.label = list(size = 10, color = "black", face = "plain", family = NULL))
  
  ggsave(plot = joint_ss_elbow, filename = sprintf("%s/kmean_n_clust.png", write_at), 
         width = 16.8, height = 8.4, units = "cm", dpi = 600)
  
  clust_plot_0 <- create_cluster_plot(mds.df = mds.df, label = "env",
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = n_clust, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  clust_plot_1 <- create_cluster_plot(mds.df = mds.df, 
                                      wss_values = wss_values,
                                      label = "env", 
                                      write_at = write_at, 
                                      n_clust = kink_fin, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  clust_plot_2 <- create_cluster_plot(mds.df = mds.df, label = "env",
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = n_ss, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  clust_plot_3 <- create_cluster_plot(mds.df = mds.df, label = "env",
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = 7, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  # Write a section in the log file
  cat("Done with deriving env_clusters",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Add clusters to original mat -----------------------------------------------
  mds.df_0 <- as.data.frame(clust_plot_0$data) %>% 
    tibble::rownames_to_column(var = "ids_mod") %>%
    left_join(residual_vals_wide %>% select(ids_mod, connect_climate, year, color), by = "ids_mod") %>%
    select(ids_mod, connect_climate, year, color, groups) %>%
    mutate(lat = as.numeric(gsub("(\\S+)_\\S+_\\S+_\\S+", "\\1", connect_climate)),
           long = as.numeric(gsub("\\S+_(\\S+)_\\S+_\\S+", "\\1", connect_climate)),
           loc = gsub("\\S+_\\S+_(\\S+)_\\S+", "\\1", connect_climate),
           year = as.numeric(gsub("\\S+_\\S+_\\S+_(\\S+)", "\\1", connect_climate)))
  
  cat("Added cluster information to original data",
      file = log_file,
      sep = "\n",
      append = TRUE)
  # Visualize the identified clusters on a physical map ------------------------
  #bw_map <- get_googlemap(center = c(9.9019, 49.8430), zoom = 4, scale = 2,
  #                        style = "feature:road|visibility:off&style=element:labels|visibility:off&style=feature#:administrative|visibility:on")
  
  #load("/proj/ext_dir/KIBREED/source_data/map.Rdata")
  
  bw_map <- qread("/proj/results/R/KIBREED_data_generation/stadia_map.qs")
  
  cluster_map <- ggmap(bw_map) +
    geom_point(aes(x = long, y = lat, color = as.factor(year)),
               data = mds.df_0, size = 1) +
    #scale_color_discrete(type = mds.df_0$color)+
    labs(x = "Longitude", y = "Latitude") +
    facet_wrap(~groups, nrow = 2, ncol = 5) +
    guides(color = guide_legend(title = "Year")) +
    scale_colour_manual(values=coul)
  
  ggsave(plot = cluster_map,
         filename = sprintf("%s/clust_map.png", write_at), 
         width = 16.4, height = 8.4, dpi = 600, units = "cm") 
  
  cat("Clusters visualized on physical map",
      file = log_file,
      sep = "\n",
      append = TRUE)
  # Check difference in clustering with different inputs -----------------------
  
  ## match the two matrices
  mat_resi <- mat # unscaled, rectangular
  mat_ec <- ec_mat_scaled[rownames(mat_resi), ] # scaled, rectangular
  
  ## get trees
  tree_data_1 <- mat_ec %>% 
    dist(method = "euclidean") %>% 
    hclust() %>% as.dendrogram() %>%
    set("branches_lwd", 1) %>%
    #set("branches_lty", c(1,1,3,1,1,2)) %>%
    set("labels_cex", 0.5)
  
  tree_data_2 <- mat_resi %>% 
    dist(method = "euclidean") %>% 
    hclust() %>% as.dendrogram() %>%
    set("branches_lwd", 1) %>%
    #set("branches_lty", c(1,1,3,1,1,2)) %>%
    set("labels_cex", .5)

  dl <- dendlist("Without GxE" = tree_data_1, 
                 "With GxE" = tree_data_2)
  
  plot_1 <- record_raster_grob(16.8, 16.8, 300, expr = {
    dl %>%
      tanglegram(#sort = TRUE,
        #common_subtrees_color_lines = FALSE,
        #common_subtrees_color_branches = TRUE,
        highlight_distinct_edges  = FALSE,
        lwd = 2,
        main = paste("Entanglement =", round(entanglement(dl), 2)),
        main_left = "Without GxE",
        main_right = "With GxE",
        cex_main = 1
      )
  })
  
  ## check entanglement 
  methods <- c("labels", "ladderize", "random", "step1side", "step2side", "DendSer")
  entg <- sapply(methods,
                 function(x, dend) {
                   untangled_tree <- dend %>%
                     untangle(method = x)
                   value <- round(entanglement(untangled_tree), 2)
                   return(value)
                 }, dend = dl)
  
  dl_2 <- dl %>% 
    untangle(method = names(which(entg == min(entg))))
  
  plot_2 <- record_raster_grob(16.8, 16.8, 300, expr = {
    dl_2 %>%
      tanglegram(#sort = TRUE,
        #common_subtrees_color_lines = FALSE,
        #common_subtrees_color_branches = TRUE,
        highlight_distinct_edges  = FALSE,
        lwd = 2,
        main = paste("Entanglement =", round(entanglement(dl_2), 2)),
        main_left = "Without GxE",
        main_right = "With GxE",
        cex_main = 1
        )
  })
  
  tree_cor <- data.frame("all.equal" = round(as.numeric(gsub(".*(\\d\\.\\d+).*", "\\1", all.equal(dl), perl = TRUE)), 2),
                         "cor.dendlist" = round(cor.dendlist(dl)[1, 2], 2),
                         "cor_bakers_gamma" = round(cor_bakers_gamma(dl)[1], 2))
  
  combined_plot <- wrap_elements(plot_1) + wrap_elements(plot_2)
  
  ggsave(plot = combined_plot, filename = sprintf("%s/tanglegram.png", write_at), 
         width = 16.8, height = 8.4, units = "cm", dpi = 600)
  
  ggsave(plot = plot_2, filename = sprintf("%s/tanglegram_resolved.png", write_at), 
         width = 16.8, height = 16.8, units = "cm", dpi = 600)
  
  cat("Tanglegram created",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Combine with EC data -------------------------------------------------------
  final_data <- ec_mat %>%
    left_join(mds.df_0, by = c("harvest_env" = "connect_climate")) %>%
    filter(!is.na(groups)) %>%
    relocate(all_of(c("ids_mod", "year", "color", "groups")), .after = "env")
  
  write.csv(final_data, file = sprintf("%s/clust_data.csv", write_at), row.names = F)
  
  cat("Final data written. Ready for random forests!",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output  -----------------------------------------------------------
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["input_mat"]] <- mat
  #out[["dend_plot"]] <- cic_plot
  out[["kmean_n_clust"]] <- kmean_n_clust
  out[["kmean_raw"]] <- wss_values
  out[["clust_plot"]] <- clust_plot_0
  out[["final_data"]] <- final_data
  out[["labs"]] <- labs
  out[["tanglegram_cor"]] <- tree_cor
  
  return(out)
}

get_feature_imp_plots <- function(existing_data, data, scores_at, ref_time){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  feature_file_time <- file.info(scores_at)$ctime
  
  if(difftime(feature_file_time, ref_time) < 0){
    msg <- "The feature importance score file is old. Maybe rerun the analysis"
    message(msg)
  }
  
  feature_imp_0 <- read.csv(scores_at, header = T)
  
  # Write a section in the log file
  cat("Start with feature importance scores",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce output
  period_to_month <- data.frame(period = c("0-30", "30-60", "60-90", 
                                           "90-120", "120-150", "150-180", 
                                           "180-210", "210-240", "240-270", 
                                           "270-300", "300-330"),
                                month = c(month.abb[10:12], month.abb[1:8]))
  
  full_names <- c(
    "air_temp_5cm_avg..C.", "air_temp_5cm_max..C.",
    "air_temp_5cm_min..C.", "air_temp_avg..C.",
    "air_temp_max..C.", "air_temp_min..C.",
    "dew_point_avg..C.", "dew_point_max..C.",
    "dew_point_min..C.", "long_wave_radiation_avg..W.m.2.",
    "long_wave_radiation_max..W.m.2.", "long_wave_radiation_min..W.m.2.",
    "pet_period..mm.", "precip_acc_period..mm.",
    "precip_acc_period_adjusted..mm.", "precip_acc_period_raw..mm.",
    "relative_humidity_avg....", "relative_humidity_max....",
    "relative_humidity_min....", "short_wave_radiation_avg..W.m.2.",
    "short_wave_radiation_max..W.m.2.", "wind_speed_2m_avg..kph.",
    "wind_speed_2m_max..kph.", "wind_speed_2m_min..kph.",
    "wind_speed_avg..kph.", "wind_speed_max..kph.",
    "wind_speed_min..kph."
  )
  short_names <- c(
    "at_5_avC", "at_5_maC",
    "at_5_miC", "at_avC",
    "at_maC", "at_miC",
    "dp_avC", "dp_maC",
    "dp_miC", "lwr_avW",
    "lwr_maW", "lwr_miW",
    "pre_MM", "pre_a_MM",
    "pre_ad_MM", "pre_r_MM",
    "rh_avP", "rh_maP",
    "rh_miP", "swr_avW",
    "swr_maW", "ws_2_avK",
    "ws_2_mxK", "ws_2_miK",
    "ws_avK", "ws_maK",
    "ws_miK"
  )
  var_names <- data.frame(variable = full_names, var_short = short_names)
  
  feature_imp <- feature_imp_0 %>%
    mutate(variable = gsub("(\\S+)\\_\\(\\S+", "\\1", old, perl = T),
           period = gsub("(\\S+)\\_\\((\\S+)\\]", "\\2", old, perl = T),
           period = gsub(",", "-", period)) %>%
    left_join(period_to_month, by = "period") %>%
    #filter(period != "BLUEs_raw") %>%
    left_join(var_names, by = "variable") %>%
    convert(fct(var_short)) %>%
    arrange(desc(gain_scaled))
  feature_imp$period <- factor(feature_imp$period, levels = period_to_month$period)
  feature_imp$month <- factor(feature_imp$month, levels = period_to_month$month)
  
  ## mean feature importances
  feature_imp_mean <- feature_imp %>%
    group_by(var_short) %>%
    summarize(mean_gain_scaled = mean(gain_scaled), .groups = "drop") %>%
    arrange(desc(mean_gain_scaled))
  feature_imp_mean$var_short <- factor(feature_imp_mean$var_short, levels = feature_imp_mean$var_short)
  
  mean_imp_plot <- ggplot(data=feature_imp_mean, aes(x=var_short, y=mean_gain_scaled)) +
    geom_bar(stat="identity") +
    labs(x = "") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  ggsave(mean_imp_plot, filename = sprintf("%s/%s", write_at, "mean_imp_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  ## feature importance per growth stage
  feature_imp$var_short <- factor(feature_imp$var_short, levels = feature_imp_mean$var_short)
  breaks_in <- c(seq(0, 0.01, 0.001), seq(0.02, 0.05, 0.01))
  imp_stage_plot <- ggplot(feature_imp %>% rename(gain = gain_scaled), 
                           aes(x = var_short, y = month, fill = gain)) +
    geom_tile() +
    labs(x = "", y = "") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    scale_fill_binned(limits=c(0, 0.04), breaks=breaks_in, type = "viridis")
  
  ggsave(imp_stage_plot, filename = sprintf("%s/%s", write_at, "imp_stage_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  # Write a section in the log file
  cat("Done with feature importance scores",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["mean_imp_plot"]] <- mean_imp_plot
  out[["imp_stage_plot"]] <- imp_stage_plot
  
  return(out)
}

get_yield_gain <- function(existing_data, data){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  write_at <- data[["write_at"]]
  predicted_data <- data[["predicted_max_yield"]]$pred
  env_cluster_information <- existing_data[["final_data"]] %>%
    select("harvest_env", "env", "ids_mod", "year", "color", "groups")
  connection <- predicted_data %>% distinct(env, connect_climate) %>%
    left_join(env_cluster_information %>% 
                select("harvest_env", "ids_mod"),
              by = c("connect_climate" = "harvest_env")) %>%
    rowwise() %>%
    mutate(site = gsub("(\\w{4})\\d{2}(\\_\\d)?", "\\1\\2", ids_mod)) %>%
    ungroup() %>%
    arrange(site) %>%
    mutate(env_no_2 = row_number()) %>%
    convert(fct(site))
  
  # Write a section in the log file
  cat("Start with getting yield gain",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Get the top performers from BLUEs across environments
  BLUEs <- qread("/proj/results/R/KIBREED_data_generation/BLUES_acr_env.qs")
  select_from <- BLUEs %>% distinct(unique_idx, .keep_all = TRUE) %>%
    filter(Type == "Line") 
  
  extremes <- predicted_data %>% filter(set_type == "extremes") %>% 
    distinct(connect_geno) %>% pull(connect_geno) %>% as.vector()
  
  randoms <- predicted_data %>% filter(set_type == "randoms") %>% 
    distinct(connect_geno) %>% pull(connect_geno) %>% as.vector()
  
  extreme_50 <- select_from %>%
    filter(connect_geno_data %in% extremes) %>%
    slice_max(order_by = BLUEs, n = 50) %>%
    pull(connect_geno_data) %>%
    as.vector() # best overall
    
  # Format data
  output <- NULL
  output_common <- NULL
  for(i in unique(predicted_data$env)){
    ## Subset data
    sub_data <- predicted_data %>% 
      filter(env == i, set == "pred_set") %>%
      arrange(desc(pred))
    
    ## Get env cluster information
    connect_env <- sub_data %>%
      distinct(connect_climate) %>%
      pull(connect_climate)
    group <- env_cluster_information %>% 
      filter(harvest_env == connect_env) %>%
      pull(groups) %>% as.numeric()
    
    ## Get needed stats
    ### For total data as baseline
    median_base <- sub_data %>% 
      filter(connect_geno %in% c(randoms, extremes)) %>% 
      pull(pred) %>%
      median() %>%
      as.numeric()
    mean_base <- sub_data %>% 
      filter(connect_geno %in% c(randoms, extremes)) %>% 
      pull(pred) %>%
      mean() %>%
      as.numeric()
    
    ### For top 10 percent in each environment
    top_per <- sub_data %>% 
      slice_max(order_by = pred, prop = 0.1)
    median_top <- top_per %>%
      pull(pred) %>%
      median() %>%
      as.numeric()
    mean_top <- top_per %>%
      pull(pred) %>%
      mean() %>%
      as.numeric()
    prop_random <- top_per %>% count(set_type) %>%
      filter(set_type == "randoms") %>%
      pull(n) %>%
      as.numeric()
    if(identical(prop_random, numeric(0))){
      prop_random <- 0
    }
    
    ### performance of extremes
    extreme_per <- sub_data %>% 
      #slice_max(order_by = pred, prop = 0.1) %>%
      filter(connect_geno %in% extreme_50)
    median_best_overall <- extreme_per %>%
      pull(pred) %>%
      median() %>%
      as.numeric()
    mean_best_overall <- extreme_per %>%
      pull(pred) %>%
      mean() %>%
      as.numeric()
    
    ## Generate output
    out <- data.frame("env" = i,
                      "group" = group,
                      "val_median_overall" = median_best_overall,
                      "val_mean_overall" = mean_best_overall,
                      "val_median_top" = median_top,
                      "val_mean_top" = mean_top,
                      "val_prop_random" = prop_random,
                      "val_median" = median_base,
                      "val_mean" = mean_base)
    if(is.null(output)){
      output <- out
    } else {
      output <- rbind(output, out)
    }
    if(is.null(output_common)){
      output_common <- top_per
    } else {
      output_common <- rbind(output_common, top_per)
    }
  }
  
  # derive additional information
  output$delta_median_base_top <- output$val_median_top - output$val_median
  output$delta_mean_base_top <- output$val_mean_top - output$val_mean
  output$delta_median_base_over <-  output$val_median_overall - output$val_median
  output$delta_mean_base_over <-  output$val_mean_overall - output$val_mean
  
  # add env some additional information
  output_mod <- output %>%
    mutate(series = gsub("(Exp_\\d).*", "\\1", env, perl = TRUE),
           year = gsub("\\S+(\\d{4}).*", "\\1", env, perl = TRUE)) %>%
    convert(fct(series)) %>%
    arrange(series, year) %>%
    mutate(env_no = row_number())
    
  #output$env_no <- 1:nrow(output)
  
  # get common geno_info
  geno_identity <- output_common %>%
    mutate(present = 1) %>%
    pivot_wider(id_cols = "env", names_from = "connect_geno", values_from = present) %>%
    left_join(output %>% select(env, group), by = "env") %>%
    relocate(group, .after = "env") %>%
    arrange(group)
  
  output_env <- matrix(NA, nrow(geno_identity), nrow(geno_identity),
                         dimnames = list(geno_identity$env, geno_identity$env))
  
  for(i in rownames(output_env)){
    data_i <- output_common %>% 
      filter(env == i) %>%
      pull(connect_geno)
    for(j in colnames(output_env)){
      data_j <- output_common %>%
        filter(env == j) %>%
        pull(connect_geno)
      output_env[i, j] <- length(intersect(data_i, data_j))
    }
  }
  
  groups <- unique(geno_identity$group)
  output_group <- matrix(NA, length(groups), length(groups),
                         dimnames = list(groups, groups))
  
  for(i in groups){
    env_grp_i <- geno_identity %>%
      select(env, group) %>%
      filter(group == i) %>%
      pull(env) %>%
      as.vector()
    for(j in groups){
      env_grp_j <- geno_identity %>%
        select(env, group) %>%
        filter(group == j) %>%
        pull(env) %>%
        as.vector()
      data_grp <- output_env[env_grp_i, env_grp_j]
      if(i == j){
        mean_val <- round(mean(data_grp[upper.tri(data_grp, diag = F)]), 1)
      } else if (j > i){
        mean_val <- round(mean(data_grp), 1)
      } else {
        mean_val <- NA
      }
      output_group[i, j] <- mean_val
    }
  }
  
  # Check for interesting cases
  for(i in groups){
    for(j in groups){
      if(i > j){
        diag_val <- output_group[i, i]
        mirror_val <- output_group[j, i]
        output_group[i, j] <- ifelse(mirror_val > diag_val, "l", NA)
      }
    }
  }
  
  # Make plots
  yield_gain <- output_mod %>%
    pivot_longer(cols = contains(c("delta_", "val_")), 
                 names_to = "measure", 
                 values_to = "value") %>%
    convert(chr(measure)) %>%
    arrange(measure) %>% 
    mutate(type = ifelse(grepl("delta_", measure), "delta",
                         ifelse(grepl("val_", measure), "val", NA)),
           stat = ifelse(grepl("_mean", measure), "mean",
                         ifelse(grepl("_median", measure), "median", 
                                ifelse(grepl("_prop", measure), "prop", NA))),
           measure = ifelse(measure == "val_median" | measure == "val_mean", 
                            "whole_data", 
                            gsub("delta_mean_|delta_median_|val_median_|val_mean_|val_prop_", "", measure)))
  
  # Get mean improvements
  #yield_gain %>%
  #  filter(type == "delta") %>%
  #  group_by(measure, stat) %>%
  #  summarize(mean_val = mean(value), .groups = "drop") %>% 
  #  filter(stat == "mean")
  
  # Get improvement summary
  #yield_gain %>%
  #  filter(type == "delta", stat == "mean") %>% 
  #  pivot_wider(id_cols = c("env", "group", "series", "year", "env_no"), 
  #              names_from = measure, values_from = value) %>% 
  #  mutate(diff = base_top - base_over) %>% 
  #  pull(diff) %>% summary()
  
  gain_mean <- yield_gain %>%
    filter(type == "delta") %>%
    group_by(measure, stat) %>%
    summarize(mean_val = mean(value), .groups = "drop") %>% 
    filter(stat == "mean") %>% 
    pull(mean_val) %>%
    round(1)
  
  yield_mean <- yield_gain %>%
    filter(type == "val",  stat == "mean", 
           measure == "whole_data")  %>%
    ggplot(aes(env_no, value, color = series)) +
    #geom_line(linewidth = 0.2) +
    geom_point(size = 0.75) +
    labs(y = "Mean yield\n(quintal per ha)", x = "Environment") +
    coord_cartesian(ylim = c(0, 140)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.25, "cm")) +
    scale_color_manual(values = c("#66A61E", "#E6AB02", "#A6761D", "#D95F02", 
                                  "#1B9E77", "#7570B3", "#E7298A")) +
    guides(color = guide_legend(title="Series",
                                override.aes = list(size=2)))
  
  # Check improvements per site
  yield_gain_aug <- yield_gain %>%
    left_join(connection, by = "env")
  
  order_yield_site <- yield_gain_aug %>%
    filter(type == "val",  stat == "mean", 
           measure == "whole_data") %>%
    group_by(site) %>%
    summarize(val_mean = mean(value)) %>%
    arrange(val_mean) %>%
    mutate(env_no_3 = row_number())
  
  order_yield_year <- yield_gain_aug %>%
    filter(type == "val",  stat == "mean", 
           measure == "whole_data") %>%
    group_by(site) %>%
    summarize(val_mean = median(value)) %>%
    arrange(val_mean) %>%
    mutate(env_no_3 = row_number())
  
  f <- function(y) {
    c(label=length(y), y=median(y))}
  
  yield_mean_2 <- yield_gain_aug %>%
    filter(type == "val",  stat == "mean", 
           measure == "whole_data")  %>%
    left_join(order_yield_year, by = "site") %>%
    arrange(env_no_3) %>%
    convert(fct(site, .args = list(levels = order_yield_year$site))) %>%
    ggplot(aes(interaction(env_no_3, site), value)) +
    geom_point() +
    #geom_line(linewidth = 0.2) +
    geom_boxplot(size = 0.75) +
    labs(y = "Mean yield\n(quintal per ha)", x = "Environment") +
    coord_cartesian(ylim = c(0, 140)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",
          legend.key.width = unit(0.25, "cm"),
          ggh4x.axis.nesttext.x = element_text(
            angle = 90, hjust = .5
          )) +
    guides(color = guide_legend(title="Series",
                                override.aes = list(size=2)),
           x = "axis_nested") +
    stat_summary(fun.data=f, geom="text", vjust=-0.5, col="blue")
  
  ggsave(yield_mean_2, filename = sprintf("%s/%s", write_at, "yield_gain_plot_site.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  exceptions <- yield_gain %>%
    filter(type == "val",  stat == "prop", measure == "random") %>%
    filter(value > 0) %>%
    arrange(desc(value)) %>%
    rename(prop = value) %>%
    convert(fct(prop))
  
  env_with_random_performers <- yield_gain %>%
    filter(type == "delta",  stat == "mean", measure == "base_top") %>%
    filter(env %in% exceptions$env) %>%
    left_join(exceptions %>% select(env, prop), by = "env")
    
  yield_delta <- yield_gain %>%
    filter(type == "delta",  stat == "mean")  %>%
    ggplot(aes(env_no, value, color = measure)) +
    #geom_line(linewidth = 0.5) +
    geom_point(aes(shape = measure), size = 0.75) +
    geom_hline(yintercept = gain_mean,
               linewidth = 0.1
               #, color = c("red", "blue")
    ) +
    #geom_point(aes(env_no, value, shape = prop), color = "red", 
    #           size = 0.5, 
    #           data = env_with_random_performers, 
    #           inherit.aes = FALSE, show.legend = FALSE) +
    #annotate("text", x=3, y=gain_mean -5, label=as.character(gain_mean)) +
    labs(y = "Selection gain\n(quintal per ha)", x = "Environment") +
    coord_cartesian(ylim = c(0, 15)) +
    guides(color = guide_legend(nrow=2, ncol = 1, byrow=TRUE,
                                override.aes = list(size=2)),
           shape = "none") +
    scale_shape_manual(values=c(20, 3))+
    scale_color_manual(labels = c("With candidates selected based on\noverall performance",
                                  "With candidates selected based on\nenvironment-wise performance"),
                       values = c("#666666", "#000080")) +
    theme_classic(base_size = 10) +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  joint_plot <- ggarrange(yield_mean, yield_delta, ncol = 2, nrow = 1,
                          align = c("hv"), 
                          labels = c("a", "b"), 
                          font.label = list(size = 10, color = "black", face = "plain", family = NULL))
  
  ggsave(joint_plot, filename = sprintf("%s/%s", write_at, "yield_gain_plot.png"),
         width = 16.8, height = 8.4, units = "cm", dpi = 600, bg = "white")

  # Write a section in the log file
  cat("Done with getting yield gain",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["write_at"]] <- write_at
  out[["assess_yield"]] <- output
  out[["assess_overlaps"]] <- output_group
  out[["yield_gain_data"]] <- yield_gain
  return(out)
}

make_overview_plots <- function(existing_data, data1, data2, data3, data4, data5){

  # Sequester data -------------------------------------------------------------
  log_file <- data4[["log_file"]]
  tmp_data <- data4[["tmp_data"]]
  write_at <- data4[["write_at"]]
  RD_mat <- data1$RD_out
  p1 <- data2$core_ext_plot + theme_classic(base_size = 9) + labs(x= "Yield (quintal per ha)", y = "Count")
  
  p2 <- data3$kmean_n_clust$data %>%
    filter(clusters <= 20) %>%
    ggplot(aes(x = clusters, y = wss)) +
    geom_line(size = 0.25) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 9, linetype="dotted", size = 0.25) +
    theme_classic(base_size = 9) + 
    labs(y= "WSS", x = "Clusters")
  
  p3_data <- data3$clust_plot$data %>% rename(Groups = groups)
  p3_labs <- data3$labs 
  
  p4 <- data4$mean_imp_plot$data %>%
    ggplot(aes(x=var_short, y=mean_gain_scaled)) +
    geom_bar(stat="identity") +
    labs(x = "") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(y = "Mean gain") +
    theme(axis.text.x = element_text(angle = 90, vjust=-0.001))
  
  p5 <- data4$imp_stage_plot +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 90, vjust=-0.001)) +
    labs(fill = "Gain")
  
  p6_data <- data5$training_data %>%
    filter(!is.na(blues))
  
  # Generate eigen decomp for RD mat.------------------------------------------- 
  # this is funneled to KIBREED data generation 
  geno <- qread(existing_data[["p_wtn"]]) %>%
    rename(connect_climate = Connect_at,
           connect_geno = connect_geno_data) %>%
    filter(BLUES_dt > 0) %>% 
    distinct(Series, connect_geno, .keep_all = TRUE) %>%
    select(Series, connect_geno, Type) %>%
    mutate(in_core = ifelse(connect_geno %in% data2$core$sel, TRUE, FALSE))
  
  ## Split geno data into lines and hybrids
  geno_lines <- geno %>% filter(Type != "Hybrid")
  geno_hybrids <- geno %>% filter(Type == "Hybrid")
  
  ## Process lines
  RD_sub_lines <- RD_mat[geno_lines$connect_geno, geno_lines$connect_geno]
  if (RhpcBLASctl::get_num_cores() > 60) {
    RhpcBLASctl::blas_set_num_threads(60)
  }
  pc_lines <- cmdscale(RD_sub_lines, k = 4, eig = TRUE, add = TRUE)
  pc_eig_lines <- pc_lines$eig
  pco_data_lines <- as.data.frame(cbind(pc_lines$points[, 1:4], "geno" = rownames(pc_lines$points))) %>%
    left_join(geno_lines, by = c("geno" = "connect_geno")) %>%
    group_by(geno) %>%
    mutate(Ser = ifelse(n() > 1, "multi", "unique")) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    convert(chr(Series)) %>%
    mutate(Series = ifelse(Ser == "multi", "multi", Series)) %>%
    select(-Ser) %>%
    convert(fct(Series), num(V1, V2, V3, V4))
  
  ## Process hybrids
  RD_sub_hybrids <- RD_mat[geno_hybrids$connect_geno, geno_hybrids$connect_geno]
  pc_hybrids <- cmdscale(RD_sub_hybrids, k = 4, eig = TRUE, add = TRUE)
  pc_eig_hybrids <- pc_hybrids$eig
  pco_data_hybrids <- as.data.frame(cbind(pc_hybrids$points[, 1:4], "geno" = rownames(pc_hybrids$points))) %>%
    left_join(geno_hybrids, by = c("geno" = "connect_geno")) %>%
    group_by(geno) %>%
    mutate(Ser = ifelse(n() > 1, "multi", "unique")) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    convert(chr(Series)) %>%
    mutate(Series = ifelse(Ser == "multi", "multi", Series)) %>%
    select(-Ser) %>%
    convert(fct(Series), num(V1, V2, V3, V4))
  
  pco_data <- list(pco_lines = pco_data_lines,
                   pco_hybrids = pco_data_hybrids,
                   eig_lines = pc_eig_lines,
                   eig_hybrids = pc_eig_hybrids)
  
  # Generate within and between series distance plots --------------------------
  # for lines and hybrids
  series_info <-  qread(existing_data[["p_wtn"]]) %>%
    convert(chr(Type)) %>%
    mutate(Type = ifelse(Type != "Hybrid", "Line", Type)) %>%
    distinct(Series, Type, connect_geno_data) %>%
    group_by(Series, Type) %>% 
    group_split(.keep = TRUE)
  names(series_info) <- as.vector(sapply(series_info, function(x) paste0(unique(x$Series), "@", unique(x$Type))))
  series_info <- lapply(series_info, function(x) as.vector(x$connect_geno_data))
  
  RD_data_ser <- list()
  
  for(ser_1 in names(series_info)){
    for(ser_2 in names(series_info)){
      mat <- RD_mat[series_info[[ser_1]], series_info[[ser_2]]]
      if(grepl("Hybrid", ser_1) && grepl("Hybrid", ser_2) | grepl("Line", ser_1) && grepl("Line", ser_2)){
        if(ser_1 != ser_2){
          vals <- as.vector(mat)
        } else {
          vals <- as.vector(mat[upper.tri(mat, diag = FALSE)])
          #vals <- as.vector(mat)
        }
        RD_data_ser[[sprintf("%s&%s", ser_1, ser_2)]] <- vals
      }
    }
  }
  
  # Create a vector of standardized pairs
  standardized_pairs <- sapply(names(RD_data_ser), standardize_pair)
  unique_pairs_indices <- !duplicated(standardized_pairs)
  #unique_pairs <- names(RD_data_ser)[unique_pairs_indices]
  #RD_data_ser_ug <- RD_data_ser[unique_pairs]
  unique_pairs <- names(RD_data_ser)
  RD_data_ser_ug <- RD_data_ser
  
  RD_data_ser_line <- RD_data_ser_ug[grepl("Line", unique_pairs)]
  names(RD_data_ser_line) <- gsub("@Line", "", names(RD_data_ser_line))
  RD_data_ser_hyb <- RD_data_ser_ug[grepl("Hybrid", unique_pairs)]
  names(RD_data_ser_hyb) <- gsub("@Hybrid", "", names(RD_data_ser_hyb))
  
  # Number of elements and size of the square grid
  series_line <- 1:7
  series_hybrid <- 1:5
  
  plot_line <- record_raster_grob(25.2, 25.2, 300, {
    plot_series(series_line, RD_data_ser_line)
    })
  
  plot_hybrid <- record_raster_grob(16.8, 16.8, 300, {
    plot_series(series_hybrid, RD_data_ser_hyb)
  })
  
  ggsave(wrap_elements(plot_line), filename = sprintf("%s/%s", write_at, "ser_dist_line.png"),
         width = 25.2, height = 25.2, units = "cm", dpi = 600, bg = "white")
  
  ggsave(wrap_elements(plot_hybrid), filename = sprintf("%s/%s", write_at, "ser_dist_hyb.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  # PCoA core set -----------------------------------------------------------
  
  pco_data_lines <- pco_data$pco_lines
  eig_lines <- pco_data$eig_lines
  pc_eig_lines <- round(100*(eig_lines[1:2]/sum(eig_lines)), 2)
  
  pco_data_lines_1 <- pco_data_lines %>% filter(!in_core) 
  
  pco_data_lines_2 <- pco_data_lines %>% filter(in_core)
  
  pco_plot <- pco_data_lines_1 %>%
    ggplot(aes(x = V1, y = V2))+
    geom_point(size = 0.05, alpha = .5) +
    geom_point(aes(x = V1, y = V2), color = "red", size = 0.05, shape = 4,
               data = pco_data_lines_2) +
    labs(x = paste0("PCo1 ", pc_eig_lines[1], " %"), 
         y = paste0("PCo2 ", pc_eig_lines[2], " %")) +
    coord_fixed(xlim = c(-0.2, 0.3), ylim = c(-0.2, 0.3)) +
    theme_classic(base_size = 9) +
    theme(plot.margin = margin(0,0,0,0, "pt"))
  
  # Make cluster plot ----------------------------------------------------------
  p3 <- ggscatter(p3_data,
                  x = "V1", 
                  y = "V2", 
                  xlab = paste0("PC1 = ", round(p3_labs[1], 2)," %"),
                  ylab = paste0("PC2 = ", round(p3_labs[2], 2)," %"),
                  #label = "env",
                  color = "Groups",
                  font.label = list(size = 5, face = "plain", family = NULL),
                  ellipse = TRUE, 
                  ellipse.type = "convex", 
                  repel = TRUE,
                  ggtheme = theme_classic(base_size = 9),
                  parse = TRUE)
  #eigen <- eigen(data1$RD_out)
  
  joint_plot_1 <- ggarrange(ggarrange(pco_plot, p1, 
                                      ncol = 1, align = "hv",
                                      labels = c("a", "b"),
                                      font.label = list(size = 10, color = "black", face = "plain", family = NULL)), 
                            p3, 
                            ncol = 2, widths = c(1, 2),
                            labels = c("", "c"),
                            font.label = list(size = 10, color = "black", face = "plain", family = NULL))
  
  ggsave(joint_plot_1, filename = sprintf("%s/%s", write_at, "core_set_joint.png"),
         width = 16.8, height = 12, units = "cm", dpi = 600, bg = "white")
  
  # Make feature importance plots ---------------------------------------------- 
  joint_plot_2 <- ggarrange(p4, p5, 
                            ncol = 2, widths = c(0.8, 1),
                            labels = c("d", "e"),
                            font.label = list(size = 10, color = "black", face = "plain", family = NULL))
  joint_plot_2_ann <- annotate_figure(joint_plot_2,
                  bottom = text_grob("Environment variable", size = 9))
  
  ggsave(joint_plot_2_ann, filename = sprintf("%s/%s", write_at, "feature_imp_joint.png"),
         width = 16.8, height = 8.4, units = "cm", dpi = 600, bg = "white")
  
  joint_plot <- ggarrange(joint_plot_1, joint_plot_2_ann, nrow = 2, align = "hv")
  
  ggsave(joint_plot, filename = sprintf("%s/%s", write_at, "feature_imp_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  
  # Plot for training data 2 ------------------------------------------------
  extremes <- p6_data %>% 
    filter(set_type == "extremes") %>%
    distinct(connect_geno) %>%
    pull(connect_geno) %>%
    as.vector()
  randoms <- p6_data %>% 
    filter(set_type == "randoms") %>%
    distinct(connect_geno) %>%
    pull(connect_geno) %>%
    as.vector()
  p6 <- ggplot(aes(x = blues, fill = set_type), data = p6_data %>% filter(set == "pred_set")) +
    geom_histogram(bins = 50, alpha = 0.75) +
    scale_fill_manual(labels = c("Extreme", 
                                 "Random"),
                      values = c("red", "blue")) +
    guides(fill=guide_legend(title="Type")) +
    labs(x= "Yield (quintal per ha)", y = "Count") +
    geom_rug(aes(x = blues),
             inherit.aes = FALSE,
             #show.legend = FALSE, 
             color = "#d0e700",
             data = p6_data %>% filter(set == "rest")) +
    theme_classic(base_size = 10)
  BLUEs <- qread("/proj/results/R/KIBREED_data_generation/BLUES_acr_env.qs")
  p7_data <- BLUEs %>% distinct(unique_idx, .keep_all = TRUE) %>%
    filter(Type == "Line") %>%
    filter(connect_geno_data %in% c(extremes, randoms)) %>%
    mutate(geno_cat = ifelse(connect_geno_data %in% extremes, "Extreme", 
                             ifelse(connect_geno_data %in% randoms, "Random", NA)))
  p7 <- ggplot(aes(x = BLUEs, fill = geno_cat), data = p7_data) +
    geom_histogram(bins = 50, alpha = 0.75) +
    scale_fill_manual(labels = c("Extreme", 
                                 "Random"),
                      values = c("red", "blue")) +
    guides(fill=guide_legend(title="Type")) +
    labs(x= "Yield (quintal per ha)", y = "Count") +
    theme_classic(base_size = 10) +
    theme(legend.position = "None")
  p6_final <- p6 + inset_element(p7, 0.05, 0.6, 0.4, 1)
  ggsave(p6_final, filename = sprintf("%s/%s", write_at, "training_data_2_plot.png"),
         width = 16.8, height = 16.8, units = "cm", dpi = 600, bg = "white")
  # Generate output ------------------------------------------------------------
  out <- list()
  out[["joint_plot_1"]] <- joint_plot_1
  out[["joint_plot_2"]] <- joint_plot_2
  out[["RD_data"]] <- RD_data_ser_ug
  out[["pco_data"]] <- pco_data
  
  return(out)
}