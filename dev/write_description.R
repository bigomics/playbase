##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "playbase") {
    stop("This file is supposed to run from the root Playbase folder")
}

if(!require("renv")) install.packages("renv")
if(!require("remotes")) install.packages("remotes")
if(!require("devtools")) install.packages("devtools")

source("dev/functions.R")
pkg <- scan_packages(path='R') 

if(file.exists("DESCRIPTION")) {
  hash <- substring(tempfile(),nchar(tempfile())-3, nchar(tempfile()))
  file.copy("DESCRIPTION",paste0("DESCRIPTION.save.",hash))
  message("WARNING: overwriting existing DESCRIPTION file.")
}

pkg.remotes <- pkg$remotes
pkg.imports <- pkg$imports

desc.header <- readLines("dev/DESCRIPTION.header")
desc.file <- "DESCRIPTION"

write(desc.header, file=desc.file)

write("Imports:", file=desc.file, append=TRUE)
write(paste0("    ",sort(pkg.imports),","), file=desc.file, append=TRUE)

write("Remotes:", file=desc.file, append=TRUE)
write(paste0("    ",sort(pkg.remotes),","), file=desc.file, append=TRUE)


