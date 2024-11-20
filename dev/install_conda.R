## ---------------------------------------------------------------------
## Install Kaleido for plotly
## ---------------------------------------------------------------------

## Install a clean reticulate and miniconda
if(!require("reticulate"))  install.packages('reticulate') 
unlink("~/.local/share/r-miniconda", recursive = TRUE)
reticulate::install_miniconda()
reticulate::conda_install("r-reticulate", "python-kaleido")
reticulate::conda_install("r-reticulate", "plotly", channel = "plotly")
reticulate::use_miniconda("r-reticulate")

## clean up downloaded packages
unlink("~/.local/share/r-miniconda/pkgs", recursive = TRUE)
