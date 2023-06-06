##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

# check the pgx matches the expected output from example data

filespath <- "dev//pgx"
filesname <- dir(filespath)

lapply(filesname, function(x){
  tmp <- playbase::pgx.load(file.path(filespath,x))
}) -> pgx_files

# check special characters

names(pgx_files) <- filesname

example_data_truth <- pgx_files[["example_data.pgx"]]

example_gx <- example_data_truth$gx.meta$meta$"act_vs_notact"$"avg.1"

example_gset <- as.numeric(example_data_truth$gset.meta$meta$"act_vs_notact"[,"fc"])


lapply(pgx_files, function(x){
  #x = pgx_files[[1]]
  gx <- cor(x$gx.meta$meta$"act_vs_notact"$"avg.1", example_gx)
  gset <- cor(as.numeric(x$gset.meta$meta$"act_vs_notact"[,"fc"]), example_gset)

  return(data.frame(gx = gx, gset = gset))

}) -> pgx_files_check

pgx_files_check <- do.call(rbind,pgx_files_check)

# save pgx_files_check as csv

write.csv(pgx_files_check, file = "dev//pgx_files_check.csv", row.names = TRUE)