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
                          ggtheme = theme_classic(base_size = 9),
                          title = sprintf("Clusters = %s", n_clust))
  
  ggsave(plot = clust_plot, filename = sprintf("%s/clust_plot_%s_%s.png", write_at, n_clust, ifelse(is.null(label), "no_label", label)), 
         width = 16.8, height = 16.8, units = "cm", dpi = 600)
  return(clust_plot)
}

get_env_clusters <- function(existing_data_path, write_at, log_at){
  # Put a log file
  run_instance <- as.character(format(Sys.time(),  format = "%d_%m_%Y_%H_%M"))
  log_at <- sprintf("%s/%s", log_at, run_instance)
  log_file <- sprintf("%s/pred_acr.log", log_at)
  create_dir_file(log_at, file = FALSE)
  
  # Put a results dir
  create_dir_file(write_at, file = FALSE)
  
  # Sequester data
  residual_vals_wide <- qread(sprintf("%s/GxE_patterns.qs", existing_data_path))
  
  # Write a section in the log file
  cat("Data loaded",
      file = log_file,
      sep = "\n",
      append = TRUE)

  # Derive relationship matrix -------------------------------------------------
  mat <- as.matrix(residual_vals_wide[, 2:ncol(residual_vals_wide)])
  row.names(mat) <- residual_vals_wide$ids_mod

  cat("Relationship matrix derived",
      file = log_file,
      sep = "\n",
      append = TRUE)
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
  
  clust_plot_0 <- create_cluster_plot(mds.df = mds.df, label = NULL,
                                      wss_values = wss_values,
                                      write_at = write_at, 
                                      n_clust = n_clust, 
                                      labs = labs, 
                                      col_vector = col_vector)
  
  # Write a section in the log file
  cat("Done with deriving env_clusters",
      file = log_file,
      sep = "\n",
      append = TRUE)

  # Generate output  -----------------------------------------------------------
  out <- list()
  out[["log_file"]] <- log_file
  out[["write_at"]] <- write_at
  out[["input_mat"]] <- mat
  out[["kmean_n_clust"]] <- kmean_n_clust
  out[["clust_plot"]] <- clust_plot_0

  return(out)
}