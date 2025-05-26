# core functions
"%!in%" <- Negate("%in%")

load_and_process_pred_files_R <- function(path){
  data <- qread(path)
  out <- list()
  
  if(grepl("cv_acr", path)){
    
    pattern <- "(\\S+)\\_(r\\_\\d+)\\_(f\\_\\d+)\\_(M\\_\\d+)"
    out[["acc"]] <- data$preds %>% filter(is.na(blues)) %>% 
      mutate(run_type = gsub(pattern, "\\1", run_name, perl = T),
             run_idx = gsub(pattern, "\\2", run_name, perl = T),
             fold_idx = gsub(pattern, "\\3", run_name, perl = T),
             model_idx = gsub(pattern, "\\4", run_name, perl = T),
             cv_idx = NA) %>%
      select(series, geno, type, idx_col, run_type, run_idx, fold_idx, model_idx, cv_idx, obs, pred)
    if(grepl("cv_acr_str|cv_acr_sce", path)){
      out[["acc"]] <- out[["acc"]] %>%
        group_by(run_type, model_idx, run_idx, type) %>% 
        summarize(cor_val = cor(pred, obs),
                  rmse_val = rmse(pred,obs),
                  .groups = "drop") %>%
        pivot_longer(cols = c("cor_val", "rmse_val"), names_to = "val_type", values_to = "value") %>%
        arrange(val_type)}
  } else if(grepl("cv_wtn", path)){
    if(grepl("tra", path)){
      pattern <- "^(cv\\_\\S+)\\_(M\\S+)\\_(run\\S+)\\_(cv\\d)$"
      data_0 <- data$preds %>% filter(is.na(blues)) %>%
        mutate(run_type = gsub(pattern, "\\1", run_name, perl = TRUE),
               model_idx = gsub(pattern, "\\2", run_name, perl = TRUE),
               run_idx = gsub(pattern, "\\3", run_name, perl = TRUE),
               cv_idx = gsub(pattern, "\\4", run_name, perl = TRUE),
               idx_col = NA,
               fold_idx = NA)
    } else {
      pattern <- "^(cv\\_\\S+)\\_(M\\S+)\\_(\\S+)\\_(run\\_\\d\\d?)$"
      data_0 <- data$preds %>% filter(is.na(blues)) %>%
        mutate(run_type = gsub(pattern, "\\1", run_name, perl = TRUE),
               model_idx = gsub(pattern, "\\2", run_name, perl = TRUE),
               run_idx = gsub(pattern, "\\4", run_name, perl = TRUE),
               cv_idx = gsub(pattern, "\\3", run_name, perl = TRUE),
               idx_col = NA,
               fold_idx = NA)
    }
  
    out[["acc"]] <- data_0 %>%
      select(env, geno, type, idx_col, run_type, run_idx, fold_idx, model_idx, cv_idx, obs, pred) %>%
      group_by(run_type, model_idx, run_idx, cv_idx, env, type) %>% # within environment
      summarize(cor_val = cor(obs, pred),
                rmse_val = rmse(obs, pred),
                .groups = "drop") %>%
      pivot_longer(cols = c("cor_val", "rmse_val"), names_to = "val_type", values_to = "value_0") %>%
      group_by(run_type, model_idx, run_idx, cv_idx, type, val_type) %>% # mean across environments
      summarize(value = mean(value_0, na.rm = T), .groups = "drop") %>%
      arrange(val_type)
  }
  
  out[["vars"]] <- data$vars
  out[["name"]] <- unique(data$preds$run_name)
  
  return(out)
}

load_and_process_pred_files_Py <- function(path){
  data <- read.csv(path) 
  out <- list()
  
  if (grepl("acr_cv", path)){
    out[["acc"]] <- data %>% select(-index, -BLUEs_raw, -BLUEs_scaled)
  } else if(grepl("acr_st|acr_sce", path)){
    out[["acc"]] <- data %>%
      group_by(run_type, model_idx, run_idx, type) %>% 
      summarize(cor_val = cor(pred, obs),
                rmse_val = rmse(pred,obs),
                .groups = "drop") %>%
      pivot_longer(cols = c("cor_val", "rmse_val"), names_to = "val_type", values_to = "value") %>%
      arrange(val_type)
  } else if(grepl("wtn", path)){
    out[["acc"]] <- data %>%
      group_by(run_type, model_idx, run_idx, cv_idx, env, type) %>% # within environment
      summarize(cor_val = cor(obs, pred),
                rmse_val = rmse(obs, pred),
                .groups = "drop") %>%
      pivot_longer(cols = c("cor_val", "rmse_val"), names_to = "val_type", values_to = "value_0") %>%
      group_by(run_type, model_idx, run_idx, cv_idx, type, val_type) %>% # mean across environments
      summarize(value = mean(value_0, na.rm = T), .groups = "drop") %>%
      arrange(val_type) %>%
      mutate(run_idx = ifelse(grepl("LoO|cvL", run_idx), gsub("LoO_|cvL_", "", run_idx), run_idx),
             cv_idx = ifelse(grepl("LoO", cv_idx), "LoO", 
                             ifelse(grepl("cvL", cv_idx), "cvL", cv_idx)))
    }
  
  return(out)
}

load_files_parallel <- function(r_files, py_files, num_cores = 10){
  registerDoParallel(cores = num_cores)
  
  # Parallelize the loading and processing of files for R
  res_R <- foreach(file = r_files) %dopar% {
                     load_and_process_pred_files_R(file)
                   }
  names(res_R) <- r_files
  
  # Parallelize the loading and processing of files for Python
  res_Py <- foreach(file =  py_files) %dopar% {
                      load_and_process_pred_files_Py(file)
                    }
  names(res_Py) <- py_files
  
  # Combine the results
  res <- c(res_R, res_Py)
  
  stopImplicitCluster()
  
  return(res)
}

load_and_process_pred_files_bahareh <- function(path, geno_mapping){
  data <- read.csv(path)
  out <- data %>% 
    left_join(geno_mapping, by = c("Genotype" = "Geno_new_mod")) %>%
    mutate(Type = as.character(Type)) %>%    # Convert 'Type' column to character type
    mutate(Type = ifelse(Type != "Hybrid", "Non_hybrid", Type)) %>%  # Replace non-"Hybrid" values with "Non_hybrid"
    mutate(run_type = "cv_wtn_LoO", 
           model_idx = "CGM_only", 
           run_idx = gsub("test_(\\d\\d?)\\.csv", "run_\\1", basename(path)), 
           cv_idx = "LoO",
           env = Env, 
           type = Type,
           obs = Yield_dt/10,
           pred = Simulated_Yield) %>%
    group_by(run_type, model_idx, run_idx, cv_idx, env, type) %>% # within environment
    summarize(cor_val = cor(obs, pred),
              rmse_val = rmse(obs, pred),
              .groups = "drop") %>%
    pivot_longer(cols = c("cor_val", "rmse_val"), names_to = "val_type", values_to = "value_0") %>%
    group_by(run_type, model_idx, run_idx, cv_idx, type, val_type) %>% # mean across environments
    summarize(value = mean(value_0, na.rm = T), .groups = "drop") %>%
    arrange(val_type) %>%

  return(out)
}

load_and_process_var_files <- function(path){
  data <- qread(path)
  vars <- data[["vars"]]
  vars$run_name <- rownames(vars)
  rownames(vars) <- NULL
  
  if(grepl("tra", path)){
    pattern <- "^(cv\\_\\S+)\\_(M\\S+)\\_(run\\S+)\\_(cv\\d)$"
    data_0 <- vars %>%
      mutate(run_type = gsub(pattern, "\\1", run_name, perl = TRUE),
             model_idx = gsub(pattern, "\\2", run_name, perl = TRUE),
             run_idx = gsub(pattern, "\\3", run_name, perl = TRUE),
             cv_idx = gsub(pattern, "\\4", run_name, perl = TRUE)) %>%
      relocate(starts_with("run"), ends_with("idx"))
  } else {
    pattern <- "^(cv\\_\\S+)\\_(M\\S+)\\_(\\S+)\\_(run\\_\\d\\d?)$"
    data_0 <- vars %>%
      mutate(run_type = gsub(pattern, "\\1", run_name, perl = TRUE),
             model_idx = gsub(pattern, "\\2", run_name, perl = TRUE),
             run_idx = gsub(pattern, "\\4", run_name, perl = TRUE),
             cv_idx = gsub(pattern, "\\3", run_name, perl = TRUE)) %>%
      relocate(starts_with("run"), ends_with("idx"))
  }
  
  return(data_0)
}

load_files_parallel_var <- function(files, num_cores = 10){
  # Create connections
  registerDoParallel(cores = num_cores)
  
  # Parallelize the loading and processing of files for R
  res <- foreach(file = files) %dopar% {
    load_and_process_var_files(file)
  }
  names(res) <- files
  
  # close connections
  stopImplicitCluster()
  
  return(res)
}

make_cv_plots <- function(data, data_means, 
                          plot_type, model_type,
                          x_lab, y_lab, 
                          #col_lab, 
                          ylim, 
                          #facet_at, 
                          text_at, text_at_y, 
                          labs_x, labs_y){
  data <- data %>% 
    droplevels() %>%
    filter(label %in% plot_type,
           name %in% model_type)
  data_means <- data_means %>%
    droplevels() %>%
    filter(label %in% plot_type,
           name %in% model_type)
  
  wtn_cv_plot <- ggplot(data, 
                        aes_string(x = x_lab, y = y_lab)) +
    geom_half_boxplot(lwd = 0.25, width = 0.5, 
                      side = "l", errorbar.draw = TRUE,
                      outlier.color = NA) +
    geom_half_point(size = 0.125, side = "r", alpha = 0.25) +
    theme_classic() +
    coord_cartesian(ylim = ylim)+
    facet_grid(rows = vars(type), cols = vars(cv_idx), scales = "free_x") + 
    #geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    #geom_hline(yintercept  = 0.5, linetype ="dashed", color = "blue") + 
    geom_text(aes_string(label = text_at, y = text_at_y), data = data_means, 
              angle = 67.5, size = 3, position = position_dodge(1), show.legend = FALSE)+
    labs(x = labs_x, y = labs_y, color = "Type") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  #ggsave(plot = wtn_cv_plot, filename = sprintf("%s/%s.png", write_at, plot_name),
  #       width = 16.8, height = 16.8, units = "cm", dpi = 600)
  
  return(wtn_cv_plot)
}

# specific functions
get_pred_file_paths <- function(existing_data,
                                log_at,
                                tmp_at){
  # Log time
  current_time <- Sys.time()
  time_stamp <- format(current_time, "D_%d_%m_T_%H_%M")
  
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/process_R_pred_data.log", log_at)
  
  if(!dir.exists(log_at)){dir.create(log_at, recursive = T)}
  
  # Put a temp_dir
  tmp_data <- sprintf("%s/%s", tmp_at, run_instance)
  
  if(!dir.exists(tmp_data)){dir.create(tmp_data, recursive = T)}
  
  # write log
  cat(sprintf("Sequestering the pred_file paths and checking their status ------------------\nTime stamp = %s", time_stamp),
      file = log_file,
      sep = "\n")
  
  # convert existing metainfo to a dataframe
  existing_meta_acr <- do.call(bind_rows, lapply(existing_data[grepl("cv_acr", names(existing_data))], function(x) x[["run_meta"]])) %>% rename(mem = memory)
  existing_meta_wtn <- do.call(bind_rows, lapply(existing_data[grepl("cv_wtn", names(existing_data))], function(x) x[["run_meta"]]))
  
  ## filters 
  to_remove_1 <- existing_meta_wtn %>%
    filter(grepl("cv1|cv2|cv4", run_id) & grepl("M_9", model_alias)) # modify this to "cv1|cv2|cv4" after 15.05.25
  to_remove_2 <- existing_meta_wtn %>%
    filter(grepl("cv2|cv3|cv4", run_id) & grepl("M_8", model_alias)) # modify this to "cv2|cv3|cv4" after 15.05.25
  
  existing_meta <- existing_meta_acr %>% bind_rows(existing_meta_wtn) %>%
    anti_join(to_remove_1, by = colnames(to_remove_1)) %>%
    anti_join(to_remove_1, by = colnames(to_remove_2))
  
  # add bahareh`s data
  by_bahareh_meta <- existing_meta_wtn %>%
    filter(grepl("LoO_run", run_id)) %>%
    distinct(run_id, .keep_all = TRUE) %>%
    mutate("model_alias" = "CGM_alone",
           "model_specs" = NA,
           "cpu" = NA,
           "mem" = NA,
           "time" = NA)
  
  if(!is.null(by_bahareh_meta)){
    existing_meta <- existing_meta %>% 
      bind_rows(by_bahareh_meta)
  }
  
  # check pred file status
  check_0 <- lapply(existing_data[names(existing_data) %!in% c("pred_paths_py", "data_bahareh_at")], function(x) x[["run_data"]][["run_data"]])
  check_vec <- unlist(check_0)
  pred_data_meta_R <- data.frame("idx" = names(check_vec), "path" = check_vec, row.names = NULL) %>%
    rowwise() %>%
    mutate(cv_type = strsplit(idx, "\\.")[[1]][1],
           run_type = strsplit(idx, "\\.")[[1]][2],
           dir_type = strsplit(idx, "\\.")[[1]][3]) %>%
    filter(dir_type == "preds") %>%
    mutate(model_name = gsub("\\S+\\_(M\\.*)", "\\1", run_type, perl = TRUE), 
           pred_paths = sprintf("%s/%s.qs", path, model_name), 
           pred_check = file.exists(pred_paths),
           run_env = "R") %>%
    select(-path, -idx, -dir_type) %>%
    ungroup()
  
  to_remove_3 <- pred_data_meta_R %>%
    filter(grepl("cv1|cv2|cv4", run_type) & grepl("M_9", model_name)) # modify this to "cv1|cv2|cv4" after 15.05.25
  to_remove_4 <- pred_data_meta_R %>%
    filter(grepl("cv2|cv3|cv4", run_type) & grepl("M_8", model_name)) # modify this to "cv2|cv3|cv4" after 15.05.25
  
  pred_data_meta_R <- pred_data_meta_R %>%
    anti_join(to_remove_3, by = colnames(to_remove_3)) %>%
    anti_join(to_remove_4, by = colnames(to_remove_4))
  
  pattern_py <- "(\\S+)\\#(\\S+)\\#(\\S+)"
  files_py <- existing_data[["pred_paths_py"]]
  if(!is.null(files_py)){
    pred_data_meta_Py <- as_tibble(cbind("file" = names(files_py), 
                                         do.call(rbind, lapply(files_py, function(x) rbind(unlist(x)))))) %>%
      pivot_longer(cols = all_of(c("base_folder", "tb_cb", "mc_cb", "pred_at", "model_at", "tmp_at")),
                   names_to = "dir_type", values_to = "path") %>%
      filter(dir_type == "pred_at") %>%
      mutate(cv_type = gsub(pattern_py, "\\2", file, perl = T),
             run_type = gsub(pattern_py, "\\3", file, perl = T),
             model_name = gsub(pattern_py, "\\1", file, perl = T),
             pred_paths = paste0(path, "/output.csv"),
             pred_check = file.exists(pred_paths),
             run_env = "Py") %>%
      select(-file, -path, -dir_type)
    
    pred_data_meta <- pred_data_meta_R %>%
      bind_rows(pred_data_meta_Py)
  } else {
    pred_data_meta <- pred_data_meta_R
  }
  
  # add bahareh`s data
  by_bahareh_status <- data.frame(pred_paths = list.files(existing_data[["data_bahareh_at"]], full.names = T),
                                  run_type = list.files(existing_data[["data_bahareh_at"]])) %>%
    mutate(run_type = gsub("test\\_(\\d\\d?)\\.csv", "LoO_run_\\1_CGM_only", run_type, perl = TRUE),
           cv_type = "cv_wtn_LoO",
           model_name = "CGM_only",
           pred_check = file.exists(pred_paths),
           run_env = "Py") %>%
    select(all_of(colnames(pred_data_meta)))
  
  if(!is.null(by_bahareh_status)){
    pred_data_meta <- pred_data_meta %>% 
      bind_rows(by_bahareh_status)
  }
  
  pred_data_meta_sum <- pred_data_meta %>%
    group_by(run_env, cv_type, model_name) %>%
    summarise(files_present = sum(pred_check),
              files_pending = n() - sum(pred_check),
              total_files = n(),
              .groups = "drop") %>%
    arrange(cv_type, model_name)
  write.table(pred_data_meta_sum,
              file = log_file,
              row.names = FALSE,
              col.names = FALSE,
              append = T)
  
  pred_data_meta_sum %>% print(n = 40)
  # Generate output
  cat("Data prepared",
      file = log_file,
      sep = "\n",
      append = T)
  
  out <- list()
  out[["log_file"]] <- log_file
  out[["tmp_data"]] <- tmp_data
  out[["pred_data_meta"]] <- pred_data_meta
  out[["existing_meta"]] <- existing_meta
  out[["pred_data_meta_overview"]] <- pred_data_meta_sum
  return(out)
}

load_pred_data <- function(data, write_at, geno_mapping_from){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  pred_data_meta <- data[["pred_data_meta"]]
  existing_meta <- data[["existing_meta"]]
  # Make write directory
  if(!dir.exists(write_at)) dir.create(write_at, recursive = TRUE)
  
  model_renames_acr <- data.frame("model_idx"  = paste0("M_", c(1:3)),
                                  "name" = c("GBLUP", "GBLUP_D", "E-GBLUP_D"))
  
  # Load data
  ## acr_cv
  pred_data_files_acr_cv <- pred_data_meta %>% filter(pred_check == TRUE) %>%
    filter(grepl("cv_acr_5f|acr_cv", pred_paths)) %>%
    pull(pred_paths) %>% as.vector()
  #acr_cv_res <- c(lapply(pred_data_files_acr_cv[grepl("results/R", pred_data_files_acr_cv)],
  #                load_and_process_pred_files_R),
  #                lapply(pred_data_files_acr_cv[grepl("results/Py", pred_data_files_acr_cv)],
  #                load_and_process_pred_files_Py))
  
  acr_cv_res <- load_files_parallel(r_files = pred_data_files_acr_cv[grepl("results/R", pred_data_files_acr_cv)],
                                    py_files = pred_data_files_acr_cv[grepl("results/Py", pred_data_files_acr_cv)])
  acr_cv_res_0 <- do.call(rbind, lapply(acr_cv_res, function(x) x[["acc"]]))
  
  if(!is.null(acr_cv_res_0)){
    acr_cv_res_df <- acr_cv_res_0 %>%
      group_by(run_type, model_idx, run_idx, type) %>% 
      summarize(cor_val = cor(pred, obs),
                rmse_val = rmse(pred, obs),
                .groups = "drop") %>%
      pivot_longer(cols = ends_with("val"), names_to = "val_type", values_to = "value") %>%
      left_join(model_renames_acr, by = "model_idx") %>%
      mutate(name = ifelse(is.na(name), model_idx, name))
    uq_names <- unique(acr_cv_res_df$name)
    uq_renames_R <- unique(model_renames_acr$name[c(1:3)])
    uq_renames_Py <- uq_names[uq_names %!in% uq_renames_R]
    acr_cv_res_df$name <- factor(acr_cv_res_df$name, levels = c(uq_renames_R, uq_renames_Py))
    
    acr_cv_res_df_means <- acr_cv_res_df %>% filter(val_type == "cor_val") %>%
      aggregate(value ~ run_type + name + type, ., function(i) round(mean(i), 3)) %>%
      convert(fct(run_type, name, type))
    acr_cv_plot <- ggplot(acr_cv_res_df %>% filter(val_type == "cor_val"),
                          aes(x = name, y = value, color = type)) +
      geom_boxplot()+
      theme_classic()+
      coord_cartesian(ylim = c(0, 1))+
      geom_hline(yintercept  = 0.75, linetype ="dashed", color = "red") + 
      geom_text(aes(label = value, y = 0.1), data = acr_cv_res_df_means, 
                angle = 90, size = 3, position = position_dodge(1), show.legend = FALSE)+
      labs(x = "", y = "correlation")
  } else {
    acr_cv_plot <- NULL
  }

  ## acr_st
  #pred_data_files_acr_st <- pred_data_meta %>% filter(pred_check == TRUE) %>%
  #  filter(grepl("cv_acr_str|acr_st", pred_paths)) %>%
  #  pull(pred_paths) %>% as.vector()
  ##acr_st_res <- c(lapply(pred_data_files_acr_st[grepl("results/R", pred_data_files_acr_st)],
  ##                load_and_process_pred_files_R),
  ##                lapply(pred_data_files_acr_st[grepl("results/Py", pred_data_files_acr_st)],
  ##                load_and_process_pred_files_Py))
  #
  #acr_st_res <- load_files_parallel(r_files = pred_data_files_acr_st[grepl("results/R", pred_data_files_acr_st)],
  #                                   py_files = pred_data_files_acr_st[grepl("results/Py", pred_data_files_acr_st)])
  #acr_st_res_df_0 <- do.call(rbind, lapply(acr_st_res, function(x) x[["acc"]]))
  #
  #if(!is.null(acr_st_res_df_0)){
  #  acr_st_res_df <- acr_st_res_df_0 %>%
  #    left_join(model_renames_acr, by = "model_idx") %>%
  #    mutate(name = ifelse(is.na(name), model_idx, name))
  #  uq_names <- unique(acr_st_res_df$name)
  #  uq_renames_R <- unique(model_renames_acr$name[c(1:3)])
  #  uq_renames_Py <- uq_names[uq_names %!in% uq_renames_R]
  #  acr_st_res_df$name <- factor(acr_st_res_df$name, levels = c(uq_renames_R, uq_renames_Py))
  #  
  #  acr_st_res_df_means <- acr_st_res_df %>% filter(val_type == "cor_val") %>%
  #    aggregate(value ~ run_type + name + type, ., function(i) round(mean(i), 3)) %>%
  #    convert(fct(run_type, name, type))
  #  acr_st_plot <- ggplot(acr_st_res_df %>% filter(val_type == "cor_val"),
  #                        aes(x = name, y = value, color = type)) +
  #    geom_boxplot()+
  #    theme_classic()+
  #    coord_cartesian(ylim = c(0, 1)) + 
  #    geom_hline(yintercept  = 0.75, linetype ="dashed", color = "red") + 
  #    geom_text(aes(label = value, y = 0.1), data = acr_st_res_df_means, 
  #              angle = 90, size = 3, position = position_dodge(1), show.legend = FALSE)+
  #    labs(x = "", y = "correlation")
  #} else {
  #  acr_st_plot <- NULL
  #}
  
  ## acr_sce
  
  scenario_overview <- read.csv("~/ext_dir/KIBREED/results_plots/scenarios_pred.csv") %>%
    mutate(idx = paste0("r_", row_number())) %>%
    select(idx, scenario)
  
  acr_sce_meta <- existing_meta %>% filter(grepl("cv_acr_sce", run_name)) %>%
    mutate(run_name = gsub("\\S+\\_(r\\_\\d+)\\_f\\_\\d+", "\\1", run_name, perl = T)) %>%
    select(run_name, test, train, val) %>% convert(num(train, val, test)) %>% distinct()
  pred_data_files_acr_sce <- pred_data_meta %>% filter(pred_check == TRUE) %>%
    filter(grepl("cv_acr_sce|acr_sce", pred_paths)) %>%
    pull(pred_paths) %>% as.vector()
  #acr_sce_res <- c(lapply(pred_data_files_acr_sce[grepl("results/R", pred_data_files_acr_sce)],
  #                       load_and_process_pred_files_R),
  #                lapply(pred_data_files_acr_sce[grepl("results/Py", pred_data_files_acr_sce)],
  #                       load_and_process_pred_files_Py))
  acr_sce_res <- load_files_parallel(r_files = pred_data_files_acr_sce[grepl("results/R", pred_data_files_acr_sce)],
                                     py_files = pred_data_files_acr_sce[grepl("results/Py", pred_data_files_acr_sce)])
  
  acr_sce_res_df_0 <- do.call(rbind, lapply(acr_sce_res, function(x) x[["acc"]]))
  
  if(!is.null(acr_sce_res_df_0)){
    acr_sce_res_df <- acr_sce_res_df_0 %>%
      left_join(model_renames_acr, by = "model_idx") %>%
      mutate(name = ifelse(is.na(name), model_idx, name)) %>%
      left_join(acr_sce_meta, by = c("run_idx" = "run_name")) %>%
      mutate(train_val = train+val) %>%
      arrange(run_type, train_val) %>%
      left_join(scenario_overview, by = c("run_idx" = "idx"))
    uq_names <- unique(acr_sce_res_df$name)
    uq_renames_R <- unique(model_renames_acr$name[c(1:3)])
    uq_renames_Py <- uq_names[uq_names %!in% uq_renames_R]
    acr_sce_res_df$name <- factor(acr_sce_res_df$name, levels = c(uq_renames_R, uq_renames_Py))
    acr_sce_res_df$scenario <- factor(acr_sce_res_df$scenario, levels = as.character(1:7))
    acr_sce_res_df$type <- factor(acr_sce_res_df$type, levels = c("Non_hybrid", "Hybrid"))
    
    plot_acr_sce_1 <- ggplot(acr_sce_res_df %>% filter(val_type == "cor_val", scenario != "1"),
                             aes(x = train_val, y = value)) +
      geom_point(size = 0.1) +
      coord_cartesian(ylim = c(-0.25, 1)) +
      facet_grid(name~type) +
      geom_smooth(method = "lm") +
      theme_classic() +
      stat_cor(p.accuracy = 0.05, r.accuracy = 0.01, 
               label.x = 0, label.y = 1,  size = 1.5)
    
    plot_acr_sce_2 <- ggplot(acr_sce_res_df %>% filter(val_type == "cor_val", scenario != "1"),
                             aes(x = train_val, y = value, color = scenario)) +
      geom_point(size = 0.1) +
      coord_cartesian(ylim = c(-0.25, 1)) +
      facet_grid(name~type) +
      geom_smooth(method = "lm") +
      theme_classic() +
      stat_cor(aes(color = scenario), 
               p.accuracy = 0.05, r.accuracy = 0.01, 
               label.x = 0, label.y = seq(0.75, 1, 0.05), size = 1.5)
    
  } else {
    plot_acr_sce_1 <- NULL
    plot_acr_sce_2 <- NULL
  }
  
  ## wtn_cv
  # Function to calculate percent differences
  calculate_percent_differences <- function(data, cols_to_calc, base_col) {
    base_col_sym <- sym(base_col)
    data %>%
      mutate(across(all_of(cols_to_calc), ~ (. - !!base_col_sym) / !!base_col_sym * 100, .names = "{col}_pct_diff"))
  }
  
  model_renames <- data.frame("model_idx"  = paste0("M_", c(1:10)), 
                              "name" = c("e&g", "e&G", 
                                         "E&G", "E&G&G_E", 
                                         "E_nl&G", "E_nl&G&G_E_nl", 
                                         "E&G&g_s", "E&G&G_S",
                                         "E&G&G_S_p", "E&G&G_S_al_p"),
                              "pred_type" = c(rep("climate_based", times = 6), rep("cgm_based", times = 4)))
  
  pred_data_file_wtn_cv <- pred_data_meta %>% filter(pred_check == TRUE) %>%
    filter(grepl("wtn", pred_paths)) %>%
    pull(pred_paths) %>% as.vector()
  
  wtn_cv_res <- load_files_parallel(r_files = pred_data_file_wtn_cv[grepl("results/R", pred_data_file_wtn_cv)],
                                    py_files = pred_data_file_wtn_cv[grepl("results/Py", pred_data_file_wtn_cv)])
  
  #wtn_cv_res <- c(lapply(pred_data_file_wtn_cv[grepl("results/R", pred_data_file_wtn_cv)],
  #                load_and_process_pred_files_R),
  #         lapply(pred_data_file_wtn_cv[grepl("results/Py", pred_data_file_wtn_cv)],
  #                load_and_process_pred_files_Py))
  
  wtn_cv_res_df_00 <- do.call(rbind, lapply(wtn_cv_res, function(x) x[["acc"]]))
  
  # add bahareh data
  data_bahareh <- do.call(rbind, lapply(pred_data_meta %>% filter(pred_check == TRUE) %>%
                                          filter(grepl("ZALF", pred_paths)) %>%
                                          pull(pred_paths) %>% as.vector(), 
                                        load_and_process_pred_files_bahareh, 
                                        geno_mapping = qread(geno_mapping_from) %>% 
                                          convert(chr(Geno_new)) %>%
                                          rowwise() %>%
                                          mutate(Geno_new_mod = ifelse(Series == "Exp_3" & Type == "Hybrid",
                                                                       paste0(MaleNew, FemaleNew),
                                                                       Geno_new)) %>%
                                          ungroup() %>%
                                          distinct(Geno_new, Geno_new_mod, Type) %>%
                                          mutate(Geno_new_mod = gsub(" ", ",", Geno_new_mod))))
  
  if(!is.null(data_bahareh)){
    wtn_cv_res_df_00 <- wtn_cv_res_df_00 %>% 
      bind_rows(data_bahareh)
  }
  
  if(!is.null(wtn_cv_res_df_00)){
    wtn_cv_res_df_0 <- wtn_cv_res_df_00 %>% 
      filter(val_type == "cor_val") %>% 
      left_join(model_renames, by = "model_idx") %>%
      mutate(name = ifelse(is.na(name), model_idx, name),
             pred_type = ifelse(is.na(pred_type) & run_type == "wtn_tra" & grepl("CGM", name), "cgm_based", 
                                ifelse(is.na(pred_type) & run_type == "wtn_tra" & !grepl("CGM", name), "climate_based", 
                                       ifelse(model_idx == "CGM_only", "cgm_based", pred_type))),
             cv_idx = ifelse(cv_idx == "LoO", "cv2",
                             ifelse(cv_idx == "cvL", "cv4", cv_idx)),
             run_type = ifelse(run_type == "wtn_tra", "cv_wtn_tra", run_type))
    
    wtn_cv_res_df_climate_based <- wtn_cv_res_df_0 %>% filter(pred_type == "climate_based",
                                                              run_type %in% c("cv_wtn_tra"),
                                                              name %in% c("e&g", "e&G", 
                                                                          "E&G", "E&G&G_E", 
                                                                          "E_nl&G", "E_nl&G&G_E_nl",
                                                                          "wtn_CNN_EC")) %>%
      mutate(label = "no_cgm_input") %>%
      select(-pred_type)
    
    cv1_bench <- wtn_cv_res_df_0 %>% filter(pred_type == "climate_based", 
                                                run_type %in% c("cv_wtn_tra"), 
                                                name %in% c("e&g", "e&G"), 
                                                cv_idx %in% c("cv1")) %>%
      bind_rows(wtn_cv_res_df_0 %>% filter(pred_type == "cgm_based", 
                                           run_type %in% c("cv_wtn_tra"), 
                                           name %in% c("E&G&g_s", "E&G&G_S", "wtn_CNN_CGM_EC"), 
                                           cv_idx %in% c("cv1")))
    
    cv2_bench <- wtn_cv_res_df_0 %>% filter(pred_type == "climate_based", 
                                            run_type %in% c("cv_wtn_LoO"), 
                                            name %in% c("e&g", "e&G", 
                                                        "E&G", "E&G&G_E", 
                                                        "E_nl&G", "E_nl&G&G_E_nl"),
                                            cv_idx %in% c("cv2")) %>%
      bind_rows(wtn_cv_res_df_0 %>% filter(pred_type == "cgm_based", 
                                           run_type %in% c("cv_wtn_LoO"), 
                                           name %in% c("E&G&g_s", "CGM_only"), 
                                           cv_idx %in% c("cv2")))
    
    cv3_bench <- wtn_cv_res_df_0 %>% filter(pred_type == "climate_based", 
                                            run_type %in% c("cv_wtn_tra"), 
                                            name %in% c("e&g", "e&G"), 
                                            cv_idx %in% c("cv3")) %>%
      bind_rows(wtn_cv_res_df_0 %>% filter(pred_type == "cgm_based", 
                                           run_type %in% c("cv_wtn_tra"), 
                                           name %in% c("E&G&g_s", "E&G&G_S_p"), 
                                           cv_idx %in% c("cv3")))
    
    cv4_bench <- wtn_cv_res_df_0 %>% filter(pred_type == "climate_based", 
                                            run_type %in% c("cv_wtn_cvL"), 
                                            name %in% c("e&g", "e&G", 
                                                        "E&G", "E&G&G_E", 
                                                        "E_nl&G", "E_nl&G&G_E_nl"), 
                                            cv_idx %in% c("cv4")) %>%
      bind_rows(wtn_cv_res_df_0 %>% filter(pred_type == "cgm_based", 
                                           run_type %in% c("cv_wtn_cvL"), 
                                           name %in% c("E&G&g_s", "E&G&G_S_al_p"), 
                                           cv_idx %in% c("cv4")))
    
    
    wtn_cv_res_df_cgm_based <- bind_rows(cv1_bench) %>%
      bind_rows(cv2_bench) %>%
      bind_rows(cv3_bench) %>%
      bind_rows(cv4_bench) %>%
      mutate(label = "with_cgm_input") %>%
      select(-pred_type)
    
    #wtn_cv_res_df_0 %>% filter(pred_type == "climate_based", 
    #run_type %in% c("cv_wtn_LoO",  "cv_wtn_cvL"), 
    #name %in% c("e&g", "e&G",
    #            "E&G&g_s", "E&G&G_S",
    #            "E&G&G_S_p", "E&G&G_S_al_p",
    #            "wtn_CNN_CGM_EC"))
    
    filter_data <- wtn_cv_res_df_0 %>% 
      filter(run_type == "cv_wtn_tra",
             cv_idx %in% c("cv2", "cv4")) %>%
      select(-pred_type)
    
    wtn_cv_res_df <- wtn_cv_res_df_climate_based %>% 
      bind_rows(wtn_cv_res_df_cgm_based) %>%
      #anti_join(filter_data, by = colnames(filter_data)) %>%
      convert(fct(label, type)) %>%
      mutate(model_idx = case_when(
        model_idx == "wtn_CNN_CGM_EC" ~ "CNN_GS",
        model_idx == "wtn_CNN_EC" ~ "CNN_EC",
        .default = as.character(model_idx)
      ),
      type = case_when(
        type == "Non_hybrid" ~ "Lines",
        .default = as.character(type)
      ))
    
    wtn_cv_res_df$cv_idx <- factor(wtn_cv_res_df$cv_idx, levels = c("cv1", "cv2", "cv3", "cv4"))
    uq_names <- unique(wtn_cv_res_df$name)
    uq_renames_R <- as.vector(unique(model_renames$name))
    uq_renames_Py <- as.vector(uq_names[uq_names %!in% uq_renames_R])
    uq_renames <- c(uq_renames_R, uq_renames_Py)
    wtn_cv_res_df$name <- factor(wtn_cv_res_df$name, levels = uq_renames)
    wtn_cv_res_df$model_idx <- factor(wtn_cv_res_df$model_idx, levels = c("M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "CNN_EC", "M_7", "M_8", "M_9", "M_10", "CNN_GS", "CGM_only"))
    base_model <- "e&G"
    to_remove <- "e&g"
    #int_models <- model_renames$name[c(3:10)]
    int_models <- uq_renames[which(uq_renames %!in% c(base_model, to_remove))]
    int_models <- int_models[which(int_models %in% wtn_cv_res_df$name)] # subset only those which are present in the loaded data
    
    wtn_cv_res_df_means <- wtn_cv_res_df %>%
      aggregate(value ~ run_type + name + model_idx + cv_idx + type + label, ., function(i) round(mean(i), 3)) %>%
      mutate(value_mod = ifelse(cv_idx == "cv2" | cv_idx == "cv4", round(value, 3), round(value, 2))) %>% 
      convert(fct(run_type, name, model_idx, cv_idx, type))
    
    idx_cols <- c("run_type", "run_idx", "cv_idx", "type", "label")
    
    wtn_cv_res_df_2 <- wtn_cv_res_df %>% 
      pivot_wider(id_cols = idx_cols, names_from = "name", values_from = "value") %>%
      calculate_percent_differences(int_models, base_model) %>%
      select(all_of(idx_cols), contains("diff")) %>%
      pivot_longer(cols = contains("diff"), names_to = "name", values_to = "value") %>%
      mutate(name = factor(gsub("_pct_diff", "", name), levels = int_models)) %>%
      left_join(wtn_cv_res_df %>% distinct(name, model_idx), by = "name")
    
    wtn_cv_res_df_means_2 <- wtn_cv_res_df_2 %>%
      aggregate(value ~ run_type + name + model_idx + cv_idx + type + label, ., function(i) round(mean(i), 2)) %>%
      mutate(value_mod = ifelse(value > 1, round(value), NA)) %>% 
      convert(fct(run_type, name, model_idx, cv_idx, type))
    
    plot_no_cgm <- make_cv_plots(data = wtn_cv_res_df, 
                                 data_means = wtn_cv_res_df_means,
                                 plot_type = c("no_cgm_input"), 
                                 model_type = c("e&g", "e&G", 
                                                "E&G", "E&G&G_E", 
                                                "E_nl&G", "E_nl&G&G_E_nl",
                                                "wtn_CNN_EC"),
                                 x_lab = "model_idx", y_lab = "value",
                                 ylim = c(-0.2, 1),
                                 text_at = "value_mod", text_at_y = 0.9,
                                 labs_x = "", labs_y = "Mean correlation")
    
    plot_no_cgm_diff <- make_cv_plots(data = wtn_cv_res_df_2, 
                                   data_means = wtn_cv_res_df_means_2,
                                   plot_type = c("no_cgm_input"), 
                                   model_type = c("e&g", "e&G", 
                                                  "E&G", "E&G&G_E", 
                                                  "E_nl&G", "E_nl&G&G_E_nl",
                                                  "wtn_CNN_EC"),
                                   x_lab = "model_idx", y_lab = "value",
                                   ylim = c(-175, 175),
                                   text_at = "value_mod", text_at_y = 175,
                                   labs_x = "", labs_y = sprintf("percent_diff_from _%s", base_model))
    
    plot_cgm <- make_cv_plots(data = wtn_cv_res_df,
                              data_means = wtn_cv_res_df_means,
                              plot_type = c("with_cgm_input"), 
                              model_type = c("e&g", "e&G",
                                             "E&G", "E&G&G_E",
                                             "E_nl&G", "E_nl&G&G_E_nl",
                                             "E&G&g_s", "E&G&G_S",
                                             "E&G&G_S_p", "E&G&G_S_al_p",
                                             "wtn_CNN_CGM_EC", "CGM_only"),
                              x_lab = "model_idx", y_lab = "value",
                              ylim = c(-0.2, 1),
                              text_at = "value", text_at_y = 0.9,
                              labs_x = "", labs_y = "Mean correlation")
    
    plot_cgm_diff <- make_cv_plots(data = wtn_cv_res_df_2 %>%
                                     filter(!is.na(value)),
                                data_means = wtn_cv_res_df_means_2,
                                plot_type = c("with_cgm_input"), 
                                model_type = c("e&g", "e&G",
                                               "E&G", "E&G&G_E",
                                               "E_nl&G", "E_nl&G&G_E_nl",
                                               "E&G&g_s", "E&G&G_S",
                                               "E&G&G_S_p", "E&G&G_S_al_p",
                                               "wtn_CNN_CGM_EC", "CGM_only"),
                                x_lab = "model_idx", y_lab = "value",
                                ylim = c(-175, 175),
                                text_at = "value_mod", text_at_y = 175,
                                labs_x = "", labs_y = sprintf("percent_diff_from _%s", base_model))
    
    #wtn_cv_plot <- ggplot(wtn_cv_res_df, 
    #                      aes(x = name, y = value, color = type)) +
    #  geom_boxplot(lwd = 0.5)+
    #  theme_classic()+
    #  coord_cartesian(ylim = c(-0.2, 1))+
    #  facet_grid(label ~ cv_idx, scales = "free_x") + 
    #  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    #  geom_hline(yintercept  = 0.5, linetype ="dashed", color = "blue") + 
    #  geom_text(aes(label = value, y = 0.9), data = wtn_cv_res_df_means, 
    #            angle = 67.5, size = 3, position = position_dodge(1), show.legend = FALSE)+
    #  labs(x = "", y = "mean_correlation") +
    #  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    #
    #wtn_cv_plot_2 <- wtn_cv_res_df_2 %>%
    #  ggplot(aes(x = name, y = value, color = type)) +
    #  geom_boxplot(lwd = 0.5)+
    #  theme_classic()+
    #  coord_cartesian(ylim = c(-75, 175))+
    #  facet_wrap(~cv_idx, ncol = 2) + 
    #  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    #  geom_text(aes(label = value_mod, y = 175), data = wtn_cv_res_df_means_2, 
    #            angle = 0, size = 3, position = position_dodge(1), show.legend = FALSE)+
    #  labs(x = "", y = sprintf("percent_diff_from _%s", base_model)) +
    #  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    #ggsave(plot = wtn_cv_plot, filename = "/proj/tmp_data/fig_1.png", width = 16.8, height = 8.4, units = "cm", dpi = 600)
    #ggsave(plot = wtn_cv_plot_2, filename = "/proj/tmp_data/fig_2.png", width = 16.8, height = 8.4, units = "cm", dpi = 600)

  } else {
    plot_no_cgm <- NULL
    plot_no_cgm_diff <- NULL
    plot_cgm <- NULL
    plot_cgm_diff <- NULL
  }
  
  # for variances
  wtn_cv_vars <- load_files_parallel_var(files = pred_data_file_wtn_cv[grepl("results/R", pred_data_file_wtn_cv)])
  wtn_cv_vars_df <- do.call(bind_rows, wtn_cv_vars) %>% 
    pivot_longer(!c(starts_with("run"), ends_with("idx")), names_to = "type", values_to = "value") %>%
    filter(grepl("var", type, ignore.case = T), !is.na(value)) %>%
    group_by(run_type, model_idx, cv_idx, type) %>%
    summarize(avg_val = mean(value)) %>%
    mutate(p_cent = avg_val/sum(avg_val)) %>%
    ungroup()
  
  ggplot(wtn_cv_vars_df,
         aes(x = model_idx, y = p_cent, fill = type)) +
    geom_bar(stat = "identity", position = "fill", colour = "black") +
    # geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5)) + # if labels are desired
    facet_grid(run_type ~ cv_idx, scales = "free_x")
  
  # data for predicted params
  param_names <- data.frame("orig" = c("Tsum"          , "BaseT"     , "SLA"         , "MAR"                        
                                       , "Vernalization", "DayLength" , "GR"          , "MR"                         
                                       , "IRD"          , "RPR"       , "RF"          , "KC"                         
                                       , "DST"          , "CTT"       , "LSenescence" , "AssimilatePartitioningCoeff"
                                       , "IPB"),
                            "param_names" = c("Tsum"   , "BaseT"    , "SLA"  ,"MAR"                        
                                              , "VRN"  , "DL"       , "GR"   , "MR"                         
                                              , "IRD"  , "RPR"      , "RF"   , "KC"                         
                                              , "DST"  , "CTT"      , "LSen" ,  "APC"
                                              , "IPB"))
  pred_param_paths_00 <- pred_data_file_wtn_cv[grepl("results/R/\\S+/cv_wtn_tra/\\S+/M_9.qs", pred_data_file_wtn_cv, perl = T)]
  pred_param_paths_0 <- list.files(gsub("preds/M_9.qs", "tmp_data/M_9_tmp/par_pred" , pred_param_paths_00), full.names = T)
  pred_param_paths <- pred_param_paths_0[grepl(".csv", pred_param_paths_0)] 
  pred_param <- do.call(rbind, lapply(pred_param_paths, read.table, header = T)) %>%
    pivot_longer(cols = c("rep", ends_with("cor_avg")), names_to = "type", values_to = "value") %>%
    mutate(param_stripped_1 = gsub("([A-Za-z]+)\\d\\d?", "\\1", param), 
           param_stripped_2 = gsub("[A-Za-z]+(\\d\\d?)", "\\1", param)) %>%
    left_join(param_names, by = c("param_stripped_1" = "orig")) %>%
    mutate(param_new = ifelse(param_stripped_1 != param_stripped_2, 
                              paste0(param_names, "-", param_stripped_2),
                              param_stripped_1)) %>%
    select(-param_names, -param_stripped_1, -param_stripped_2)
  
  pred_param$type <- factor(pred_param$type, levels = c("rep", "train_set_cor_avg", "test_set_cor_avg"))
  pred_param_means <- pred_param %>%
    aggregate(value ~ param_new + type, ., function(i) format(round(mean(i), 3), nsmall = 2))
  type.labs <- c("Repetabilities", "Mean train set correlation", "Mean test set correlation")
  names(type.labs) <- c("rep", "train_set_cor_avg", "test_set_cor_avg")
  pred_param_plot <- pred_param %>%
    ggplot(aes(x = param_new, y = value)) +
    coord_cartesian(ylim = c(-0.5, 1.4))+
    geom_boxplot(lwd = 0.65, outlier.size = 0.5) +
    geom_text(aes(label = value, y = 1.25), data = pred_param_means, 
              angle = 90, size = 3, position = position_dodge(1), show.legend = FALSE) +
    facet_wrap(~type, 
               ncol = 1,
               labeller = labeller(type = type.labs)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
    labs(x = "Parameters", y = "Values")

  # Produce output
  out <- list()
  out[["log_file"]] <- log_file
  out[["pred_data_meta"]] <- pred_data_meta 
  out[["param_meta"]] <- pred_param_means 
  out[["plot_acr_cv"]] <- acr_cv_plot
  #out[["plot_acr_st"]] <- acr_st_plot
  out[["plot_acr_sce_1"]] <- plot_acr_sce_1
  out[["plot_acr_sce_2"]] <- plot_acr_sce_2
  out[["wtn_cor_data"]] <- wtn_cv_res_df
  out[["plot_wtn_no_cgm"]] <- plot_no_cgm
  out[["plot_wtn_no_cgm_diff"]] <- plot_no_cgm_diff
  out[["plot_wtn_with_cgm"]] <- plot_cgm
  out[["plot_wtn_with_cgm_diff"]] <- plot_cgm_diff
  out[["plot_wtn_pred_param"]] <- pred_param_plot
  
  # write_files
  na_check <- is.na(out)
  plots <- grep("plot", names(na_check[which(!na_check)]), value = T)
  for (plot_name in plots){
    if (grepl("wtn", plot_name)) {
        ht = 16.8; wh = 16.8
    } else if (grepl("acr_sce", plot_name)) {
      ht = 16.8; wh = 16.8
    } else {ht = 8.4; wh = 16.8}
    plot_data <- out[[plot_name]]
    ggsave(plot = plot_data, filename = sprintf("%s/%s.png", write_at, plot_name),
           width = wh, height = ht, units = "cm", dpi = 600)
  }
  return(out)
}

make_plots <- function(data, write_at){
  # Sequester data
  log_file <- data[["log_file"]]
  data_no_cgm <- data[["plot_wtn_no_cgm"]]$plot_env$data %>%
    convert(chr(type, model_idx)) %>%
    mutate(type = ifelse(type == "Hybrid", "Hybrids", type),
           model_idx = ifelse(model_idx == "CNN_EC", "CNN_EV", model_idx)) %>%
    convert(fct(type)) %>%
    convert(fct(model_idx, 
                .args = list(levels = c("M_1", "M_2", "M_3", 
                                        "M_4", "M_5", "M_6", 
                                        "M_7", "M_8", "CNN_EV",
                                        "CNN_GS"))))
  
  data_no_cgm_means <- data[["plot_wtn_no_cgm"]]$plot_env$data_means %>%
    convert(chr(type, model_idx)) %>%
    mutate(type = ifelse(type == "Hybrid", "Hybrids", type),
           model_idx = ifelse(model_idx == "CNN_EC", "CNN_EV", model_idx)) %>%
    convert(fct(type)) %>%
    convert(fct(model_idx, 
                .args = list(levels = c("M_1", "M_2", "M_3", 
                                        "M_4", "M_5", "M_6", 
                                        "M_7", "M_8", "CNN_EV",
                                        "CNN_GS"))))
  
  data_cgm <- data[["plot_wtn_with_cgm"]]$plot_env$data %>%
    mutate(type = ifelse(type == "Hybrid", "Hybrids", type))
  data_cgm_means <- data[["plot_wtn_with_cgm"]]$plot_env$data_means %>%
    convert(chr(type)) %>%
    mutate(type = ifelse(type == "Hybrid", "Hybrids", type)) %>%
    convert(fct(type))
  
  data_param <- data[["plot_wtn_pred_param"]]$data

  # Make write directory
  if(!dir.exists(write_at)) dir.create(write_at, recursive = TRUE)
  
  # Make plots
  cv1_data <- data_no_cgm %>%
    convert(chr(cv_idx, model_idx)) %>%
    filter(cv_idx == "cv1") %>%
    bind_rows(data_cgm %>%
                filter(cv_idx == "cv1",
                       model_idx %!in% c("M_1", "M_2"))) %>%
    convert(fct(cv_idx)) %>%
    convert(fct(model_idx, 
                .args = list(levels = c("M_1", "M_2", "M_3", 
                                        "M_4", "M_5", "M_6", 
                                        "M_7", "M_8", "CNN_EV",
                                        "CNN_GS"))))
  cv1_means <- data_no_cgm_means %>%
    convert(chr(cv_idx, model_idx)) %>%
    filter(cv_idx == "cv1") %>%
    bind_rows(data_cgm_means %>%
                filter(cv_idx == "cv1",
                       model_idx %!in% c("M_1", "M_2"))) %>%
    convert(fct(cv_idx)) %>%
    convert(fct(model_idx, 
                .args = list(levels = c("M_1", "M_2", "M_3", 
                                        "M_4", "M_5", "M_6", 
                                        "M_7", "M_8", "CNN_EV",
                                        "CNN_GS"))))
  
  plot_main <- ggplot(cv1_data, aes(x = model_idx, y = value)) +
    #geom_half_boxplot(lwd = 0.25, width = 0.5, 
    #                  side = "l", errorbar.draw = TRUE,
    #                  outlier.color = NA) +
    #geom_half_point(size = 0.125, side = "r", alpha = 0.25) +
    geom_jitter(size = 0.125, width = 0.25, alpha = 0.25) +
    theme_classic() +
    coord_cartesian(ylim = c(-0.2, 1))+
    facet_grid(cols = vars(type), 
               #cols = vars(cv_idx), 
               scales = "free_x") + 
    geom_text(aes_string(label = "value_mod", y = 0.9), 
              data = cv1_means, 
              angle = 67.5, size = 3, 
              position = position_dodge(1), show.legend = FALSE)+
    labs(x = "", y = "Mean prediction ability", color = "Type") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  plot_main_2 <- ggplot(data_no_cgm %>%
                          filter(cv_idx != "cv1")
                        , aes(x = model_idx, y = value)) +
    #geom_half_boxplot(lwd = 0.25, width = 0.5, 
    #                  side = "l", errorbar.draw = TRUE,
    #                  outlier.color = NA) +
    #geom_half_point(size = 0.125, side = "r", alpha = 0.25) +
    geom_jitter(size = 0.125, width = 0.25, alpha = 0.25) +
    theme_classic() +
    coord_cartesian(ylim = c(-0.2, 1))+
    facet_grid(rows = vars(type), 
               cols = vars(cv_idx), 
               scales = "free_x") + 
    geom_text(aes_string(label = "value_mod", y = 0.9), 
              data = data_no_cgm_means %>%
                filter(cv_idx != "cv1"), 
              angle = 67.5, size = 3, 
              position = position_dodge(1), show.legend = FALSE)+
    labs(x = "", y = "Mean prediction ability", color = "Type") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  # Create an empty placeholder plot for '#'
  plot_hash <- ggplot() + theme_void()
  
  # Define layout with areas, 3 columns, 2 rows:
  layout <- "
#AA
BBB
"
  
  # Since 'B' should span all 3 bottom columns, redefine layout as area layout:
  layout_area <- c(
    area(t = 1, l = 1, b = 1, r = 1),  # # at (1,1)
    area(t = 1, l = 2, b = 1, r = 3),  # A spanning cols 2 & 3 at row 1 (two copies of plot_main combined)
    area(t = 2, l = 1, b = 2, r = 3)   # B spanning all 3 columns row 2
  )
  
  # Compose the patchwork with repeated plot_main for two top cols
  plot_main_join <- plot_hash + plot_main + plot_main_2 +
    plot_layout(design = layout_area)
    #plot_annotation(tag_levels = 'a') & 
    #theme(
      #plot.tag = element_text(size = 10, color = "black", face = "plain", family = NULL),
      #plot.tag.position = c(0, 1)
    #)
  
  ggsave(plot = plot_main_join, filename = sprintf("%s/%s.png", write_at, "cv_joint"),
         width = 17, height = 20, units = "cm", dpi = 600)
  
  plot_supp <- ggplot(data_cgm %>%
                        filter(cv_idx %!in% c("cv1", "cv3", "cv4")) %>%
                        filter(model_idx %!in% c("M_7", "M_10")) %>%
                        convert(chr(cv_idx)) %>%
                        mutate(cv_idx = ifelse(cv_idx == "cv2", "LoO", 
                                               ifelse(cv_idx == "cv4", "cvL", cv_idx))) %>%
                        convert(fct(cv_idx, .args = list(levels = c("LoO", "cvL"))))
                      , aes(x = model_idx, y = value)) +
    geom_half_boxplot(lwd = 0.25, width = 0.5, 
                      side = "l", errorbar.draw = TRUE,
                      outlier.color = NA) +
    geom_half_point(size = 0.125, side = "r", alpha = 0.25) +
    theme_classic() +
    coord_cartesian(ylim = c(-0.2, 1))+
    facet_grid(rows = vars(type), 
               cols = vars(cv_idx),
               scales = "free_x") + 
    geom_text(aes_string(label = "value_mod", y = 0.9), 
              data = data_cgm_means %>%
                filter(cv_idx %!in% c("cv1", "cv3", "cv4")) %>%
                filter(model_idx %!in% c("M_7", "M_10")) %>%
                convert(chr(cv_idx)) %>%
                mutate(cv_idx = ifelse(cv_idx == "cv2", "LoO", 
                                       ifelse(cv_idx == "cv4", "cvL", cv_idx))) %>%
                convert(fct(cv_idx, .args = list(levels = c("LoO", "cvL"))))
              , angle = 67.5, size = 3, 
              position = position_dodge(1), show.legend = FALSE)+
    labs(x = "", y = "Mean prediction ability", color = "Type") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  
  ggsave(plot = plot_supp, filename = sprintf("%s/%s.png", write_at, "cv_joint_supp"),
         width = 12.75, height = 12.75, units = "cm", dpi = 600)
  
  # plot for poster IWC 2024
  tags <- c("no G\nor E", "with G\n", " with G\n+ E\n+ GxE")
  poster_data <- data_no_cgm %>%
    filter(cv_idx %!in% c("cv1", "cv2", "cv4"),
           model_idx %in% c("M_1", "M_2", "M_6")) %>%
    convert(chr(model_idx)) %>%
    mutate(model_idx = case_when(
      model_idx == "M_2" ~ tags[2],
      model_idx == "M_6" ~ tags[3],
      .default = tags[1]
    )) %>%
    convert(fct(model_idx,
                .args = list(levels = tags)))
  
  poster_data_means <- data_no_cgm_means %>%
    filter(cv_idx %!in% c("cv1", "cv2", "cv4"),
           model_idx %in% c("M_1", "M_2", "M_6")) %>%
    convert(chr(model_idx)) %>%
    mutate(model_idx = case_when(
      model_idx == "M_2" ~ tags[2],
      model_idx == "M_6" ~ tags[3],
      .default = tags[1]
    )) %>%
    convert(fct(model_idx,
                .args = list(levels = tags)))
  
  plot_poster <- ggplot(data = poster_data
                      , aes(x = model_idx, y = value)) +
    geom_half_boxplot(lwd = 0.25, width = 0.5, 
                      side = "l", errorbar.draw = TRUE,
                      outlier.color = NA) +
    geom_half_point(size = 0.25, side = "r", alpha = 0.25) +
    theme_classic() +
    coord_cartesian(ylim = c(-0.2, 1))+
    facet_grid(cols = vars(type), 
               scales = "free_x") + 
    geom_text(aes_string(label = "value_mod", y = 0.9), 
              data = poster_data_means, 
              angle = 0, size = 3, 
              position = position_dodge(1), show.legend = FALSE)+
    labs(x = "", y = "Mean prediction ability", color = "Type") +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8))
  
  ggsave(plot = plot_poster, filename = sprintf("%s/%s.png", write_at, "cv_poster_IWC"),
         width = 8.4, height = 8.4, units = "cm", dpi = 600)
  
  # plot params
  plot_data_mod <- data_param %>% filter(type == "rep") %>% 
    mutate(process = case_when(
      grepl("APC", param_new) ~ "Assimilate partitioning",
      grepl("BaseT", param_new) ~ "Phenology &\nGrowth",
      grepl("CTT", param_new) ~ "Heat\nstress",
      grepl("DL", param_new) ~ "Phenology",
      grepl("DST", param_new) ~ "Drought stress",
      grepl("GR", param_new) ~ "Canopy growth",
      grepl("IPB", param_new) ~ "Canopy growth",
      grepl("IRD", param_new) ~ "Root\ngrowth",
      grepl("KC", param_new) ~ "Drought stress",
      grepl("LSen", param_new) ~ "Canopy growth",
      grepl("MAR", param_new) ~ "Canopy growth",
      grepl("MR", param_new) ~ "Canopy growth",
      grepl("RF", param_new) ~ "Canopy growth",
      grepl("RPR", param_new) ~ "Root\ngrowth",
      grepl("SLA", param_new) ~ "Canopy growth",
      grepl("Tsum", param_new) ~ "Phenology",
      grepl("VRN", param_new) ~ "Phenology"
    )
    ) %>%
    group_by(cv_id, param_new, process, type) %>%
    summarise(mean_value = mean(value), .groups = "drop") %>%
    convert(fct(type, process, param_new))

  create_plot <- function(data, param_range) {
    data %>%
      filter(param_new %in% levels(param_new)[param_range]) %>%
      ggplot(aes(x = interaction(param_new, process), y = mean_value)) +
      coord_cartesian(ylim = c(0, 0.5)) +
      #geom_boxplot(lwd = 0.25, outlier.size = 0.125) +
      geom_half_boxplot(lwd = 0.25, width = 1, 
                        side = "l", errorbar.draw = TRUE,
                        outlier.color = NA) +
      labs(x = "", y = "Mean repeatability") +
      guides(x = "axis_nested") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
            ggh4x.axis.nesttext.y = element_text(
              angle = 45),
            plot.margin = margin(c(0,0,0,0, unit = "cm")))
  }
  
  param_ranges <- list(1:28, 29:56)
  
  new_arrangement <- plot_data_mod %>% 
    distinct(process, param_new) %>%
    arrange(process, param_new) %>%
    pull(param_new) %>%
    as.vector()
  
  plot_data_mod_leveled <- plot_data_mod %>% 
    convert(fct(param_new, .args = list(levels = new_arrangement))) 
  
  plots <- lapply(param_ranges, function(range) create_plot(plot_data_mod_leveled, range))
  
  p1 <- plots[[1]]
  p2 <- plots[[2]]
  #p3 <- plots[[3]]
  #p4 <- plots[[4]]
  
  library(patchwork) #looks like a bug between ggh4x and patchwork
  joint_plot <- p1/p2 +
    plot_layout(axis_titles = "collect")
  
  ggsave(plot = joint_plot, filename = sprintf("%s/%s.png", write_at, "param_supp"),
         width = 17, height = 12.75, units = "cm", dpi = 600)
  
  # Create output
  out <- list()
  out[["plot_main"]] <- plot_main
  out[["plot_main_2"]] <- plot_main_2
  out[["plot_supp"]] <- plot_supp
  out[["plot_poster"]] <- plot_poster
  out[["plot_param_supp"]] <- plot_data_mod

  return(out)
}

change_job_meta_data <- function(data, new_time, 
                                 cv_type_vec, model_name_vec){
  # Sequester data
  log_file <- data[["log_file"]]
  tmp_data <- data[["tmp_data"]]
  pred_data_meta <- data[["pred_data_meta"]]
  
  # Load data and filter using cv_type_vec and model_name_vec
  pred_data_files <- pred_data_meta %>%
    filter(cv_type %in% cv_type_vec, model_name %in% model_name_vec) %>%
    mutate(
      base_dir = gsub("(\\S+)\\/cv\\S+", "\\1", pred_paths, perl = TRUE),
      master_file_srun = sprintf("%s/%s/master_files/srun_model_%s.sh", base_dir, cv_type, model_name),
      master_file_cc = sprintf("%s/%s/master_files/cc_bash_model_%s.sh", base_dir, cv_type, model_name),
      time = new_time,
      model_name = gsub("(.*)\\_M\\_\\d", "\\1", run_type, perl = TRUE)
    )
  
  # For srun master files
  unique_master_files_srun <- unique(pred_data_files$master_file_srun)
  
  out <- list()
  
  for (file in unique_master_files_srun){
    pred_data_files_subset <- pred_data_files %>% filter(master_file_srun == file)
    out[[file]]$complete <- pred_data_files_subset %>% filter(pred_check) %>% pull(model_name) %>% as.vector()
    out[[file]]$incomplete <- pred_data_files_subset %>% filter(!pred_check) %>% pull(model_name) %>% as.vector()
    raw_file_srun <- readLines(file)
    
    for(idx in out[[file]]$complete){
      target_line <- grep(idx, raw_file_srun, fixed = TRUE)
      if (!startsWith(raw_file_srun[target_line], "#")) {
        raw_file_srun[target_line] <- paste0("#", raw_file_srun[target_line])
      }
    }
    
    for (idx in out[[file]]$incomplete){
      target_line <- grep(idx, raw_file_srun, fixed = TRUE)
      raw_file_srun[target_line] <-  gsub("\\d+-\\d+:\\d+:\\d+", new_time, raw_file_srun[target_line])
    }
    cat(raw_file_srun, file = file, sep = "\n")
  }
  
  # For cc bash master files
  unique_master_files_cc <- unique(pred_data_files$master_file_cc)
  
  for (file in unique_master_files_cc){
    pred_data_files_subset <- pred_data_files %>% filter(master_file_cc == file)
    out[[file]]$complete <- pred_data_files_subset %>% filter(pred_check) %>% pull(model_name) %>% as.vector()
    out[[file]]$incomplete <- pred_data_files_subset %>% filter(!pred_check) %>% pull(model_name) %>% as.vector()
    
    raw_file_cc <- readLines(file)
    
    for(idx in out[[file]]$complete){
      target_line <- grep(idx, raw_file_cc, fixed = TRUE)
      if (!startsWith(raw_file_cc[target_line], "#")) {
        raw_file_cc[target_line] <- paste0("#", raw_file_cc[target_line])
      }
    }
    cat(raw_file_cc, file = file, sep = "\n")
  }
  
  # Produce output
  out[["log_file"]] <- log_file
  out[["change_log"]] <- out
}

# Make plots 
#plot_5f_str <- ggplot(pred_data_acc %>% 
#                        filter(grepl("cv_acr_5f|cv_acr_str", run_type, perl = TRUE),
#                               val_type == "cor_val"), 
#                      aes(x = run_type, y = value, color = Type)) + # for str calculate correlation in groups of series?
#  geom_boxplot() +
#  theme_classic(base_size = 11) +
#  coord_cartesian(ylim = c(0, 1))
#
##plot_sce <- ggplot(pred_data_acc %>% filter(grepl("cv_acr_sce", run_type, perl = TRUE)), 
##                      aes(x = run_type, y = value, color = Type)) +
##  geom_boxplot() +
##  facet_wrap(~val_type, scales = "free") + 
##  theme_classic(base_size = 11)
#
#wtn_data <- pred_data_acc %>% 
#  filter(grepl("cv_wtn_5f", run_type, perl = TRUE),
#         val_type == "cor_val") %>%
#  mutate(run_id = gsub("(run\\_\\d+)\\_(cv\\d)", "\\1", run_idx, perl = TRUE),
#         cv_id = gsub("(run\\_\\d+)\\_(cv\\d)", "\\2", run_idx, perl = TRUE)) %>%
#  convert(fct(cv_id))
#
#wtn_data_means <- wtn_data %>%
#  aggregate(value ~ run_type + model_idx + cv_id + Type, ., function(i) round(mean(i), 3)) %>%
#  convert(fct(run_type, model_idx, cv_id, Type))
#
#plot_wtn <- ggplot(wtn_data, aes(x = model_idx, y = value, color = Type)) +
#  geom_boxplot() +
#  coord_cartesian(ylim = c(0,1)) +
#  facet_wrap(~cv_id, scales = "free") + 
#  theme_classic(base_size = 11)+
#  geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
#  geom_hline(yintercept  = 0.5, linetype ="dashed", color = "blue")+
#  geom_text(aes(label = value, y = 0.1), data = wtn_data_means, angle = 90, size = 3, position = position_dodge(1), show.legend = FALSE)