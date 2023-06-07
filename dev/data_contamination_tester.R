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

example_data_truth <- pgx_files[[ "data-test__example_data.pgx"]]

example_gx <- example_data_truth$gx.meta$meta$"act_vs_notact"$"avg.1"

example_gset <- as.numeric(example_data_truth$gset.meta$meta$"act_vs_notact"[,"fc"])

contrast_to_test <- "act_vs_notact"

pgx_files_check <- data.frame(filesname = filesname)

lapply(1:length(pgx_files), function(idx){
  #idx = 28
  print(idx)
  x = pgx_files[[idx]]
  filesname[idx]
  contrast_default <- contrast_to_test

  if(!contrast_to_test %in% names(x$gx.meta$meta)) {
    contrast_default <- grep(contrast_default, names(x$gx.meta$meta), value = TRUE, ignore.case = TRUE)
  }
  if(!as.character(contrast_default)[1] %in% names(x$gx.meta$meta)) {
    contrast_default <- names(x$gx.meta$meta)[length(x$gx.meta$meta)]
  }
  gx <- cor(x$gx.meta$meta[[contrast_default]]$"avg.1", example_gx)
  gset <- cor(as.numeric(x$gset.meta$meta[[contrast_default]][,"fc"]), example_gset)

  return(data.frame(gx = gx, gset = gset))

}) -> res

res <- do.call(rbind,res)


pgx_files_check <- cbind(pgx_files_check, res)

# save pgx_files_check as csv

write.csv(pgx_files_check, file = "dev//pgx_files_check.csv", row.names = TRUE)
