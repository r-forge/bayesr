if(file.exists("~/svn/bayesr/pkg/minify"))
  file.remove("~/svn/bayesr/pkg/minify", recursive = TRUE)
dir.create("~/svn/bayesr/pkg/minify")
file.copy("~/svn/bayesr/pkg/bamlss", "~/svn/bayesr/pkg/minify/",
  overwrite = TRUE, recursive = TRUE)
setwd("~/svn/bayesr/pkg/minify/bamlss/R/")
Rfiles <- dir()
for(f in Rfiles) {
  rf <- readLines(f)
  rf <- as.character(parse(text = rf))
  rf <- gsub("        ", " ", rf)
  rf <- gsub("    ", "  ", rf)
  rf <- paste(rf, collapse = "\n")
  writeLines(rf, f)
}

