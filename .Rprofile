if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}

source(sprintf("%s/renv/activate.R", project_path))
Sys.setenv(LC_ALL = "C")
options(install.packages.compile.from.source = "never")
