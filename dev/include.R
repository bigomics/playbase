if(!file.exists("include.R")) {
  message("ERROR: please call with source(..., chdir=TRUE)")
} else {
  message("sourcing playbase code...")
  rfiles <- dir("../R", pattern="*.r$|*.R$")
  for(r in rfiles) source(file.path("../R",r),local=TRUE,encoding='UTF-8',echo=FALSE,verbose=FALSE)
}
