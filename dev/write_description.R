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
pkg <- scan_packages(path=c('R'))

desc.file <- "DESCRIPTION"

if(!file.exists(desc.file)) {
    stop("DESCRIPTION file not found. Please create one first.")
}

desc.lines <- readLines(desc.file)

imports.start <- grep("^Imports:", desc.lines)
if(length(imports.start) == 0) {
    remotes.line <- grep("^Remotes:", desc.lines)
    if(length(remotes.line) > 0) {
        imports.start <- remotes.line[1]
        imports.end <- imports.start - 1
    } else {
        imports.start <- length(desc.lines) + 1
        imports.end <- imports.start - 1
    }
} else {
    imports.start <- imports.start[1]
    subsequent.fields <- grep("^[A-Za-z]", desc.lines[(imports.start+1):length(desc.lines)])
    if(length(subsequent.fields) > 0) {
        imports.end <- imports.start + subsequent.fields[1] - 1
    } else {
        imports.end <- length(desc.lines)
    }
}

new.desc <- c()
if(imports.start > 1) {
    new.desc <- desc.lines[1:(imports.start-1)]
}

pkg.imports <- sort(pkg$imports)
new.desc <- c(new.desc, "Imports:")
new.desc <- c(new.desc, paste0("    ", pkg.imports, ","))

if(imports.end < length(desc.lines)) {
    new.desc <- c(new.desc, desc.lines[(imports.end+1):length(desc.lines)])
}

writeLines(new.desc, desc.file)
message("DESCRIPTION file updated with ", length(pkg.imports), " imports.")
