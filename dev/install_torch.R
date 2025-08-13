
# increasing timeout is recommended since we will be downloading a 2GB file.
options(timeout = 600) 

# For Windows and Linux: "cpu", "cu124" are the only currently supported
# For MacOS the supported are: "cpu-intel" or "cpu-m1"
kind <- "cu124"
version <- available.packages()["torch","Version"]
options(repos = c(
  torch = sprintf("https://torch-cdn.mlverse.org/packages/%s/%s/", kind, version),
  CRAN = "https://cloud.r-project.org" # or any other from which you want to install the other R dependencies.
))
install.packages("torch")
