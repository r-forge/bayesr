dir.create("src/bayesxsrc")
dirs <- c("adaptiv", "alex", "andrea", "bib", "dag",
  "graph", "leyre", "mcmc", "psplines", "samson", "structadd")
files <- c("export_type.h", "main.cpp", "values.h")

tdir <- tempfile()
dir.create(tdir)
owd <- getwd()
setwd(tdir)
system("git clone https://gitlab.gwdg.de/bayesx/bayesx.git")
for(d in dirs) {
  file.copy(file.path("bayesx", d), file.path(owd, "src", "bayesxsrc"), overwrite = TRUE, recursive = TRUE)
}
for(f in files) {
  file.copy(file.path("bayesx", f), file.path(owd, "src", "bayesxsrc"), overwrite = TRUE, recursive = TRUE)
}
setwd(owd)
file.copy(file.path("src", "dev-Makefile"), file.path("src", "Makefile"))
file.copy(file.path("src", "dev-Makefile.win"), file.path("src", "Makefile.win"))
file.copy(file.path("src", "dev-Makefile.inner.win"), file.path("src", "Makefile.inner.win"))
