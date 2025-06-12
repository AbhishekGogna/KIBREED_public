# Get the project directory
project_dir <- getwd()
lockfile_path <- file.path(project_dir, "renv.lock")

cat("Checking for renv.lock in:", lockfile_path, "\n")

# Check if renv.lock exists
if (!file.exists(lockfile_path)) {
  stop("Error: renv.lock file not found. Please provide a valid renv.lock to proceed.\n")
}

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv in bare mode if not already done
if (!file.exists(file.path(project_dir, "renv"))) {
  renv::init(project = project_dir, bare = TRUE, force = TRUE)
}

# Restore the environment from the lock file
cat("Restoring R packages from renv.lock...\n")
Sys.setenv(LC_ALL = "C")
options(install.packages.compile.from.source = "never")
renv::restore(lockfile = lockfile_path, prompt = FALSE)

cat("All packages restored successfully from renv.lock.\n")

# Some house keeping
if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd()
}

# Content to write to .Rprofile
rprofile_content <- 'if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}
source(sprintf("%s/renv/activate.R", project_path))
Sys.setenv(LC_ALL = "C")
options(install.packages.compile.from.source = "never")'

# Write to .Rprofile in the project directory
rprofile_path <- file.path(project_path, ".Rprofile")
writeLines(rprofile_content, rprofile_path)