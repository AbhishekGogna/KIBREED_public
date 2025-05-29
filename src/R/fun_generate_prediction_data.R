execute <- TRUE

# Core functions
"%!in%" <- Negate("%in%")

nested_list_to_dataframe <- function(nested_list) {
  
  # Flatten the nested list
  flattened_list <- unlist(nested_list, recursive = FALSE)
  split_names <- strsplit(names(flattened_list), ".", fixed = T)
  
  # Combine the level and length values into a dataframe
  df <- data.frame(run_id = sapply(split_names, "[", 1),
                   type = sapply(split_names, "[", 2),
                   size = sapply(flattened_list, function(x) {
                      if(is.character(x)){
                        return(x)
                      } else {
                        return(length(x))
                      }
                   }),
                   stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df_0 <- df %>% filter(grepl("run_name", type)) %>% select(-type) %>% 
  rename(run_name = size)
  out <- df %>% filter(!grepl("run_name", type)) %>% 
          pivot_wider(names_from = type, values_from = size) %>% 
          left_join(df_0, by = c("run_id"))
  return(out)
}

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

get_random_string <- function(length) {
  characters <- c(letters, LETTERS, 0:9)
  random_string <- paste0(sample(characters, length, replace = TRUE), collapse = "")
  return(random_string)
}

write_run_data <- function(run_meta,
                           write_path, # relative path of base_directory
                           run_script_at, # relative path of run_script
                           user_dir){ # absolute path switch
  
  # Chunk for bash master files
  bash_cmd <- 'nohup Rscript %s "%s" "%s" "%s" "%s" "%s" > %s/cc_%s.log 2> %s/cc_%s.err &' # for the file to run from inside the container
  
  #define base directories
  master_scr_at <-  sprintf("%s/master_files", write_path)
  run_data_at <- sprintf("%s/run_data", write_path)
  
  create_dir_file(c(master_scr_at, run_data_at), file = FALSE)
  
  # Produce scripts
  scripts_data <- list()
  master_script_data <- list()
  
  for(i in 1:nrow(run_meta)){
    # Fetch variables
    cv_type <- gsub(".*/(\\w+)", "\\1", write_path, perl = TRUE)
    cv_id <- run_meta$run_id[i]
    model_name <- run_meta$model_alias[i]
    model_spec <- run_meta$model_specs[i]
    
    # Define the bash master file
    bash_file <- sprintf("%s/cc_bash_model_%s.sh", master_scr_at, model_name)
    master_script_data[[paste0("bash_", model_name)]] <- create_dir_file(bash_file) ## store its data to later make it executable
    
    # Define run_folder structure for the present row of run_meta
    save_info <- list()
    cv_run_folder <- paste0(sprintf("%s/%s", run_data_at, cv_id))
    create_dir_file(cv_run_folder, file = FALSE)
    for (sub_dirs in c("tmp_data", "preds", "logs")){
      sub_dir_path <- paste0(cv_run_folder, "/", sub_dirs)
      save_info[[sub_dirs]] <- sub_dir_path   
      create_dir_file(sub_dir_path, file = FALSE)
    }
    
    ## Store the info for later use
    scripts_data[[paste0(cv_id, "_", model_name)]] <- save_info
    
    # Write lines to the bash master file
    ## Add header
    if(identical(readLines(bash_file), character(0))) {
      cat("#!/usr/bin/env bash", 
          file = bash_file, 
          sep = "\n")
    }
    
    ## Add one line for each row of meta file
    logs_at <- paste0(cv_run_folder, "/logs")
    bash_scr_to_write <- sprintf(bash_cmd, run_script_at, cv_type, cv_id, model_name, model_spec, "TRUE", logs_at, model_name, logs_at, model_name)
    bash_script <- readLines(bash_file)
    if(bash_scr_to_write %!in% bash_script){
      cat(bash_scr_to_write,
          file = bash_file,
          sep = "\n",
          append = T)
    }
  }
  
  # Make hidden files if they do not exist
  rprofile_file <- paste0(user_dir, "/.Rprofile")
  if(!file.exists(paste0(user_dir, "/.Rprofile"))) {
    file.create(rprofile_file, recursive = TRUE)
    cat('if(getwd() == "', user_dir, '") source("renv/activate.R")\n', file = rprofile_file, sep = "")
  }
  
  for (master_script_path in names(master_script_data)) {
    file_path <- master_script_data[[master_script_path]]
    file_check <- readLines(file_path)
    to_write <- 'wait'
    if(to_write %!in% file_check){
      cat('wait\necho "All jobs done!"', file = file_path, append = T, sep = "\n")
    }
  }
  
  # Make master script executable
  system(sprintf("chmod +x %s/*.sh", master_scr_at))
  
  # prepare output
  output <- list()
  output[["run_data"]] <- scripts_data
  output[["master_scripts"]] <- master_script_data
  return(output)
}

generate_run_scripts <- function(data, 
                                 run_script_at, 
                                 input_data_at, 
                                 model_info,
                                 wtn = FALSE,
                                 project_path){
  # Sequester data
  log_file <- data[["log_file"]]
  write_at <- data[["write_at"]]
  run_info <- data[["run_info"]]
  run_type <- data[["run_type"]]
  
  # filter run_types 
  # if(length(grep("cv1|cv3", names(run_info))) != 0) run_info <- run_info[grepl("cv1|cv3", names(run_info))]
  
  # Generate metadata for the cv data
  if (wtn) {
    run_meta <- data[["run_meta"]]
    # Produce plots
    plot_cvs_long <- run_meta %>% select(connect, run, cv, train, test) %>% 
      pivot_longer(cols = c("test", "train"), names_to = "type", values_to = "val") %>%
      convert(fct(run, cv, type))
    
    plot <- ggplot(plot_cvs_long, aes(x = run, y = val, color = type))+
      geom_line(aes(group=type))+
      facet_wrap(~cv)+
      theme_classic(base_size = 14)+
      theme(axis.text.x = element_blank(),
            axis.ticks = element_blank())
    ggsave(plot = plot, filename = sprintf("%s/%s_sizes.png",  write_at, run_type), width = 16.4, height = 16.4, units = "cm", dpi = 100)
    
    # Write txt files
    write.table(run_meta, sprintf("%s/%s_meta.txt", write_at, run_type), row.names = FALSE)
    
    # Write json file
    exportJson <- toJSON(run_info)
    write(exportJson, sprintf("%s/%s.json", write_at, run_type))
  }
  
  run_meta_0 <- data.frame("run_id" = rep(names(run_info), times = nrow(model_info)), 
                         "model_alias" = rep(unique(model_info$models), each = length(run_info))) 
  
  run_meta <- run_meta_0 %>% 
    left_join(model_info, by = c("model_alias" = "models")) %>%
    left_join(nested_list_to_dataframe(run_info), by = "run_id") 
  
  to_remove_1 <- run_meta %>%
    filter(grepl("cv1", run_id) & grepl("M_9", model_alias)) # modify this to "cv1|cv2|cv4" after 15.05.25
  to_remove_2 <- run_meta %>%
    filter(grepl("cv3", run_id) & grepl("M_8", model_alias)) # modify this to "cv2|cv3|cv4" after 15.05.25
  
  if(nrow(to_remove_1) != 0){
    run_meta <- run_meta %>% anti_join(to_remove_1, by = colnames(to_remove_1))
  }
  if(nrow(to_remove_2) != 0){
    run_meta <- run_meta %>% anti_join(to_remove_2, by = colnames(to_remove_2))
  }

  # Create a results dir if it does not exist
  if(!dir.exists(write_at)){dir.create(write_at, recursive = T)}
  
  # Write a section in the log file
  cat(sprintf("Generating run scripts for %s ------------------", run_type),
      file = log_file,
      sep = "\n",
      append = TRUE)

  # Create run data
  run_data <- write_run_data(run_meta = run_meta,
                             write_path = write_at, # relative path of base_directory
                             run_script_at = run_script_at,
                             user_dir = project_path) # relative path of run_script
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- write_at
  out[["run_meta"]] <- run_meta
  out[["run_data"]] <- run_data
  return(out)
}

compute_and_save_matrix <- function(matrix_name, matrix_expr, geno, write_at) {
  matrix_save_at <- sprintf("%s/%s.qs", write_at, matrix_name)
  
  if (!file.exists(matrix_save_at)) {
    tryCatch({
      matrix_obj <- list()
      matrix_obj$kin <- matrix_expr
      matrix_obj$eigen <- eigen(matrix_expr)
      matrix_obj[["eigen"]] <- c(matrix_obj[["eigen"]], "geno" = list(geno))
      matrix_obj$name <- matrix_name
      matrix_obj$matrix_save_at <- matrix_save_at
      
      qsave(matrix_obj, matrix_save_at)
    }, error = function(e) {
      stop(paste("Error occurred while computing and saving", matrix_name))
    })
  } else {
    matrix_obj <- qload(matrix_save_at)
  }
  
  return(matrix_obj)
}

perfom_eigen_decompositions <- function(g_data_kin,
                                        g_aa_mat_kin,
                                        d_mat,
                                        log_at,
                                        write_at){
  # Put a log file
  if(!file.exists(log_at)){file.create(log_at, recursive = T)}
  
  if(!dir.exists(write_at)) {
    # Put a results folder
    dir.create(write_at, recursive = TRUE)
    
    # Compute and save G_a_mat
    cat("Eigen decomp G_a -------------------",
        file = log_at,
        sep = "\n",
        append = T)
    
    G_a_mat <- compute_and_save_matrix(
      "evd_a_mat",
      g_data_kin,  
      rownames(g_data_kin),
      write_at
    )
    
    cat("Eigen decomp G_aa -------------------",
        file = log_at,
        sep = "\n",
        append = T)
    
    # Compute and save G_aa_mat
    G_aa_mat <- compute_and_save_matrix(
      "evd_aa_mat",
      g_aa_mat_kin,
      rownames(g_aa_mat_kin),
      write_at
    )
    
    cat("Eigen decomp G_d -------------------",
        file = log_at,
        sep = "\n",
        append = T)
    
    # Compute and save G_d_mat
    G_d_mat <- compute_and_save_matrix(
      "evd_d_mat",
      d_mat,
      rownames(d_mat),
      write_at
    )
    
    cat("Eigen dedup complete -------------------",
        file = log_at,
        sep = "\n",
        append = T)
  } else {
    G_a_mat <- qread(sprintf("%s/evd_a_mat.qs", write_at))
    G_aa_mat <- qread(sprintf("%s/evd_aa_mat.qs", write_at)) 
    G_d_mat <- qread(sprintf("%s/evd_d_mat.qs", write_at))
  }
  
  
  file_list <- c(G_a_mat$matrix_save_at, G_aa_mat$matrix_save_at, G_d_mat$matrix_save_at)
  write_check <- all(file.exists(unlist(file_list)))
  
  order_check_a <- sum(G_a_mat[["eigen"]][["geno"]] == rownames(g_data_kin)) == nrow(g_data_kin)
  order_check_d <- sum(G_d_mat[["eigen"]][["geno"]] == rownames(d_mat)) == nrow(d_mat)
  order_check <- all(order_check_a,
                     order_check_d)
  
  if(write_check){
    print("EVD files exits")
    if(order_check){
      print("EVD files have same order of genotypes as provided g data")
    } else {
      print("EVD files does not have the same order of genotypes as provided g data")
    }
  } else {
    print("EVD files does not exits")
  }
    
  # Generate output
  out <- list()
  out[["write_check"]] <- write_check
  out[["evd_paths"]] <- file_list
  
  return(out)
}

# Acr
create_pred_acr_objects <- function(existing_data_path,
                                    write_at,
                                    log_at,
                                    tmp_at,
                                    subset){
  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/pred_acr.log", log_at)
  create_dir_file(log_at, file = FALSE)
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_at, run_instance)
  create_dir_file(tmp_data, file = FALSE)
  
  
  cat("Creating prediction objects for acr environment prediction ------------------",
      file = log_file,
      sep = "\n")
  
  # Put a results dir
  create_dir_file(write_at, file = FALSE)
  
  # Sequester data
  pheno_data_acr_full <- qread(sprintf("%s/pheno_data_acr.qs", existing_data_path))
  g_a_data <- qread(sprintf("%s/grm_a.qs", existing_data_path))
  g_aa_data <- qread(sprintf("%s/grm_aa.qs", existing_data_path))
  g_d_data <- qread(sprintf("%s/grm_d.qs", existing_data_path))
  # space for any last minute alliterations
  
  # Subset data
  if(subset){
    pheno_data_acr <- pheno_data_acr_full %>%
      slice_sample(n = 1000, replace = FALSE)
  } else {
    pheno_data_acr <- pheno_data_acr_full
  }
  
  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- write_at
  out[["pheno_data_acr"]] <- pheno_data_acr
  out[["g_a_data"]] <- g_a_data
  out[["g_aa_data"]] <- g_aa_data
  out[["g_d_data"]] <- g_d_data
  
  return(out)
}

## CVs

cv_acr_5f <- function(data, runs, folds){
  # Sequester data
  log_file <- data[["log_file"]]
  write_at <- data[["write_at"]]
  pheno_data <- data[["pheno_data_acr"]]
  g_a_data <- data[["g_a_data"]]
  g_aa_data <- data[["g_aa_data"]]
  g_d_data <- data[["g_d_data"]]
  
  # Produce output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- sprintf("%s/cv_acr_5f", write_at)

  # Generate write at folder
  if(!dir.exists(out[["write_at"]])){dir.create(out[["write_at"]], recursive = T)}

  # Write a section in the log file
  cat("Creating 5 fold cross validation data ------------------",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce data
  
  n_idx <- length(unique(pheno_data$unique_idx))
  check_for_presence <- sum(pheno_data$unique_idx %in% 1:n_idx)
  
  if (execute & check_for_presence) {
    save_data <- list()

    for (run in 1:runs) {
      fold_no <- 0

      if (run > 1) {
        rm(kf)
      }

      set.seed(run)
      
      kf <- cvFolds(n = n_idx, K = folds, R = runs, type = "random")

      for (n in 1:ncol(kf$subsets)){
        for (f in 1:folds){
          test_set <- kf$subsets[which(kf$which == f), n]
          train_set_raw <- kf$subsets[which(kf$which != f), n]
          test_set_geno <- pheno_data %>% filter(unique_idx %in% test_set) %>% 
            pull(Geno_new) %>% unique()
          train_set <- pheno_data %>% filter(unique_idx %in% train_set_raw) %>% 
            filter(Geno_new %!in% test_set_geno) %>% pull(unique_idx) %>%
            unique()

          val_set <- sample(train_set, as.integer(0.2 * length(train_set)))
          train_set_final <- train_set[which(train_set %!in% val_set)]
          test_set_final <- test_set[which(test_set %!in% val_set)]

          save_data[[paste0("r_", run, "_f_", f)]] <- list(
            "train" = as.integer(train_set_final),
            "val" = as.integer(val_set),
            "test" = as.integer(test_set_final),
            "run_name" = paste0("cv_acr_5f_", "r_", run, "_f_", f)
          )
        }
      }
    }
    
    # Generate eigen decomposition
    geno_order <- pheno_data %>% distinct(unique_idx, connect_geno_data) %>%
      #filter(row_number() < 101) %>%
      pull(connect_geno_data) %>% as.vector()
    eigen_paths <- perfom_eigen_decompositions(g_data_kin = g_a_data[geno_order, geno_order],
                                               g_aa_mat_kin = g_aa_data[geno_order, geno_order],
                                               d_mat = g_d_data[geno_order, geno_order],
                                               log_at = log_file,
                                               write_at = sprintf("%s/eigen_data", out[["write_at"]]))
    
    # Add generated data to output object
    out[["run_info"]] <- save_data
    out[["run_type"]] <- "cv_acr_5f"
    
    # Write json file
    exportJson <- toJSON(out[["run_info"]])
    write(exportJson, sprintf("%s/cv_acr_5f.json", out[["write_at"]]))
  }
  return(out)
}

cv_acr_sce <- function(data, take_parts){
  # Sequester data
  log_file <- data[["log_file"]]
  write_at <- data[["write_at"]]
  pheno_data <- data[["pheno_data_acr"]]
  g_a_data <- data[["g_a_data"]]
  g_aa_data <- data[["g_aa_data"]]
  g_d_data <- data[["g_d_data"]]
  
  # Produce output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- sprintf("%s/cv_acr_sce", write_at)
  
  # Generate write at folder
  if(!dir.exists(out[["write_at"]])){dir.create(out[["write_at"]], recursive = T)}
  
  # Put a log statement
  cat("Creating sce cross validation data ------------------",
      file = log_file,
      sep = "\n",
      append = TRUE)
  
  # Produce data
  if (execute) {
    save_data <- list()
    
    pheno_data_mod <- pheno_data %>% mutate(Series_mod = gsub("(Exp_\\d).*", "\\1", Series,   perl = TRUE)) 
    series <- sort(unique(pheno_data_mod$Series_mod))
    
    series_numbers <- 1:length(series)
    
    scenarios <- NULL
    
    for (i in series_numbers){
      vals <- cbind("test" = i, "train" = series_numbers)
      vals <- vals[vals[,"test"] != vals[, "train"], ]
      extension <- matrix(NA, nrow(vals), length(series_numbers))
      for(j in 1:nrow(vals)){
        extension[j, vals[j, "test"]] <- 1
        extension[j, vals[j, "train"]] <- 2
      }
      scenarios <- rbind(scenarios, cbind(1, extension))
    } # one exp series as train all others as train in a given instance
    
    colnames(scenarios) <- c("scenario", series)
    
    for (sce in series_numbers[1:(length(series_numbers)-1)]){
      scenario_train <- t(combn(series_numbers, sce))
      
      extension <- matrix(NA, nrow(scenario_train), length(series_numbers))
      
      for (i in 1:nrow(extension)){
        extension[i, scenario_train[i, ]] <- 1
        extension[i, -scenario_train[i, ]] <- 2
      }
      
      scenarios <- rbind(scenarios, cbind((sce+1), extension))
    } # 2 to total number of series experimental series as train and others as test in a  given instance. 
    
    scenarios <- as_tibble(scenarios) %>% mutate(row = row_number())
    
    for (run in 1:nrow(scenarios)){
      scenario_data <- scenarios %>% filter(row == run) %>%
        pivot_longer(!c(scenario, row), names_to = "series", values_to = "values")
      train_series <- scenario_data %>% filter(values == 1) %>% pull(series) %>% as.character ()
      test_series <- scenario_data %>% filter(values == 2) %>% pull(series) %>% as.character()
      
      test_set_data <- pheno_data_mod %>% filter(Series_mod %in% test_series)
      test_set <- test_set_data %>% pull(idx_with_series)
      train_set <- pheno_data_mod %>% filter(Series_mod %in% train_series) %>% 
        filter(Geno_new %!in% test_set_data$Geno_new) %>% pull(idx_with_series)
      
      if(take_parts){
        val_set <- c(sample(train_set, as.integer(0.16 * length(train_set))), 
                     sample(test_set, as.integer(0.04 * length(test_set))))
      } else {
        val_set <- sample(train_set, as.integer(0.2 * length(train_set)))
      }
      
      train_set_final <- train_set[which(train_set %!in% val_set)]
      test_set_final <- test_set[which(test_set %!in% val_set)]
      
      save_data[[paste0("r_", run, "_f_", run)]] <- list(
        "train" = as.integer(train_set_final),
        "val" = as.integer(val_set),
        "test" = as.integer(test_set_final),
        "run_name" = paste0("cv_acr_sce_", "r_", run, "_f_", run)
      )
    }
    
    # Generate eigen decomposition
    geno_order <- pheno_data %>% distinct(idx_with_series, connect_geno_data) %>%
      #filter(row_number() < 101) %>%
      pull(connect_geno_data) %>% as.vector()
    
    eigen_paths <- perfom_eigen_decompositions(g_data_kin = g_a_data[geno_order, geno_order],
                                               g_aa_mat_kin = g_aa_data[geno_order, geno_order],
                                               d_mat = g_d_data[geno_order, geno_order],
                                               log_at = log_file,
                                               write_at = sprintf("%s/eigen_data", out[["write_at"]]))
    
    # Add generated data to output object
    out[["run_info"]] <- save_data
    out[["run_type"]] <- "cv_acr_sce"
    
    # Write json file
    exportJson <- toJSON(out[["run_info"]])
    write(exportJson, sprintf("%s/cv_acr_sce.json", out[["write_at"]]))
  } # can be moved to generate run_scripts
  
  return(out)
}

# Wtn
create_pred_wtn_objects <- function(existing_data_path,
                                    write_at,
                                    log_at,
                                    tmp_at,
                                    subset){  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/pred_wtn.log", log_at)
  create_dir_file(log_at, file = FALSE)
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_at, run_instance)
  create_dir_file(tmp_data, file = FALSE)
  
  cat("Creating prediction objects for wtn environment prediction ------------------",
      file = log_file,
      sep = "\n")

  # Sequester data
  pheno_data_wtn_full <- qread(sprintf("%s/pheno_data_wtn.qs", existing_data_path))
  
  # Subset data
  if(subset == TRUE){
    # Get 100 unique values from the connect_geno column
    unique_genos <- pheno_data_wtn_full %>%
      pull(connect_geno) %>%
      unique() %>%
      sample(1000)
    
    # Then filter the dataframe to only include those values
    pheno_data_wtn <- pheno_data_wtn_full %>%
      filter(connect_geno %in% unique_genos)
  } else {
    pheno_data_wtn <- pheno_data_wtn_full
  }

  # Generate output
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- write_at
  out[["pheno_data_wtn"]] <- pheno_data_wtn
  return(out)
}

## CVs
cv_wtn_tra <- function(data, runs, test_prop){
    
  # Sequester data
  log_file <- data[["log_file"]]
  write_at <- data[["write_at"]]
  pheno_data <- data[["pheno_data_wtn"]]
  
  # Produce output
  output <- list()
  output[["log_file"]] <- log_file
  output[["write_at"]] <- sprintf("%s/cv_wtn_tra", write_at)

  # Generate write at folder
  if(!dir.exists(output[["write_at"]])){dir.create(output[["write_at"]], recursive = T)}
  
  # Write log
  cat("Creating wtn cross validation data ------------------",
      file = log_file,
      sep = "\n")
  
  # Produce data
  if (execute){
    geno <- unique(pheno_data$Geno_new)
    ngeno <- length(geno)
    env <- unique(pheno_data$Env)
    nenv <- length(env)

    out <- list()
    set.seed(NULL)
    out_2 <- NULL
    for(i in 1:runs){
      outdata <- list()
      # test geno
      test_geno <- as.vector(geno[sample(1:ngeno, (test_prop)*ngeno)])
      # test env
      ## figure out the distribution of test genotypes in the data
      test_geno_dist <- pheno_data %>% filter(Geno_new %in% test_geno) %>% 
        count(Env) %>%
        arrange(n) %>%
        filter(n > 50) # select only those environments which have more than 50 genotypes
      target_test_env <- test_geno_dist %>% pull(Env)
      n_target_test_env <- length(target_test_env)
      test_env <- as.vector(target_test_env[sample(1:n_target_test_env, round((test_prop) *nenv))])
      # train env and geno
      train_geno <- setdiff(geno, test_geno)
      train_env <- as.vector(env[which(env %!in% test_env)])
      
      quad_4_index <- pheno_data %>%
        filter(Geno_new %in% test_geno & Env %in% test_env) %>%
        pull(idx_cv) %>% as.vector()
      
      quad_3_index <- pheno_data %>%
        filter(Geno_new %in% test_geno & Env %in% train_env) %>%
        pull(idx_cv) %>% as.vector()
      
      quad_2_index <- pheno_data %>%
        filter(Geno_new %in% train_geno & Env %in% test_env) %>%
        pull(idx_cv) %>% as.vector()
      
      quad_1_index <- pheno_data %>%
        filter(Geno_new %in% train_geno & Env %in% train_env) %>%
        pull(idx_cv) %>% as.vector()
      
      for(j in c("cv1", "cv2", "cv3", "cv4")){
        instance_name <- paste0("run_", i, "_", j)

        # create common train data and specific test data for cv1
        test_train_split <- 0.2
        split_by_env <- pheno_data %>%
          filter(idx_cv %in% quad_1_index)
        env_to_consider <- split_by_env %>% count(Env)
        geno_to_consider <- split_by_env %>% count(Geno_new)
        env_thr <- round(test_train_split*nrow(env_to_consider))
        geno_thr <- round(test_train_split*nrow(geno_to_consider))
        critical_thr <- round(quantile(env_to_consider$n)[[3]]) #round(geno_thr/env_thr)
        split_by_env_sub <- env_to_consider %>% filter(n > critical_thr) %>% 
          pull(Env)
        valid_test_data <- split_by_env %>% filter(Env %in% split_by_env_sub)
        test_set_cv1 <- valid_test_data[sample(1:nrow(valid_test_data),   test_train_split*length(quad_1_index)), ] %>% pull(idx_cv) %>% as.vector()
        # here i count genotypes per environment, remove those environments with counts less  than 50th quantile and finally sample 20% from the data corresponding to the selected environments.
        outdata[["train"]] <- quad_1_index[quad_1_index %!in% test_set_cv1]

        test_cv1_env_subset <- NA

        if(j == "cv1"){
          ##cv1
          outdata[["test"]] <- test_set_cv1
          test_cv1_env_subset <- length(split_by_env_sub)
          #sum(unique(c(outdata$train, outdata$test)) %in% quad_1_index) # sanity check
          env_counts <- pheno_data %>% filter(idx_cv %in% test_set_cv1) %>% count(Env)
          env_with_low_geno <- env_counts %>% filter(n <  50)
          if(nrow(env_with_low_geno) >= 1){
            cat(sprintf("Test set genotype counts in cv1 %s", instance_name),
                file = log_file,
                sep = "\n",
                append = T)
            write.table(env_counts,
                        file = log_file,
                        append = T)
          }
          if(nrow(env_with_low_geno) != 0){print(paste0(instance_name, " has ", nrow  (env_with_low_geno), " environments in test set with low number of genotypes"))}
        }
        if(j == "cv2"){
          ##cv2
          outdata[["test"]] <- sort(quad_2_index)
        }
        if(j == "cv3"){
          ##cv3
          outdata[["test"]] <- sort(quad_3_index)
        }
        if(j == "cv4"){
          ##cv4
          outdata[["test"]] <- sort(c(quad_4_index))
        }
        outdata[["run_name"]] <- instance_name
        
        

        ## generate output
        meta_info <- data.frame("connect" = instance_name,
                                "cv" = j,
                                "run" = i, 
                                "train" = length(outdata[["train"]]),
                                "test" = length(outdata[["test"]]),
                                "train_env" = pheno_data %>% filter(idx_cv %in% outdata[["train"]]) %>% pull(Env)  %>% unique() %>% length(),
                                "train_geno" = pheno_data %>% filter(idx_cv %in% outdata[["train"]]) %>% pull(Geno_new) %>% unique() %>% length(),
                                "test_env" = pheno_data %>% filter(idx_cv %in% outdata[["test"]]) %>% pull(Env)  %>% unique() %>% length(),
                                "test_geno" = pheno_data %>% filter(idx_cv %in% outdata[["train"]]) %>% pull(Geno_new) %>% unique() %>% length())
        out[[instance_name]] <- outdata
        out_2 <- rbind(out_2, meta_info)
      }
    }
    
    # Produce output
    output[["run_info"]] <- out
    output[["run_type"]] <- "cv_wtn_tra"
    output[["run_meta"]] <- out_2

    return(output)
  } else {
    print("No cross-validation data was generated")
    return(NULL)
  }
}