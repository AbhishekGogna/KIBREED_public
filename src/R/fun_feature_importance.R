# Core functions
"%!in%" <- Negate("%in%")

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