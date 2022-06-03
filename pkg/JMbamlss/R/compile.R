compile <- function(dir = NULL, tdir = NULL) {
  hold <- getwd()
  on.exit(setwd(getwd()))
  if(is.null(dir)) dir <- "~/svn/bayesr/pkg/JMbamlss/src"
  dir <- path.expand(dir)
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir), add = TRUE)
  }
  tdir <- path.expand(tdir)
  setwd(tdir)
  if(file.exists(file.path(dir, "Makevars")))
    file.copy(file.path(dir, "Makevars"), file.path(tdir, "Makevars"))
  cf <- grep(".c", dir(dir), value = TRUE, fixed = TRUE)
  cf <- cf[!grepl(".c~", cf, fixed = TRUE)]
  cf <- cf[!grepl("init.c", cf, fixed = TRUE)]
  for(j in cf) {
    file.copy(file.path(dir, j), file.path(tdir, j))
    system(paste("R CMD SHLIB", j))
    dyn.load(gsub(".c", ".so", j, fixed = TRUE))
  }
  setwd(hold)
}

compile_alex <- function(location = c("workstation", "laptop", "server_linux",
                                      "server_windows"), 
                         dir = NULL, tdir = NULL) {
  location <- match.arg(location)
  hold <- getwd()
  on.exit(setwd(getwd()))
  if(is.null(dir)) {
    dir <- switch(
      location, 
      "workstation" = "~/Documents/joint_models/bayesr/pkg/JMbamlss/src",
      "laptop" = "~/Dokumente/Arbeit/Greven/JM/JMbamlss/src",
      "server_linux" = "H:/volkmana.hub/JMbamlss/src",
      "server_windows" = "H:/JMbamlss/src")
  }
  dir <- path.expand(dir)
  if(is.null(tdir)) {
    dir.create(tdir <- tempfile())
    on.exit(unlink(tdir), add = TRUE)
  }
  tdir <- path.expand(tdir)
  setwd(tdir)
  if(file.exists(file.path(dir, "Makevars")))
    file.copy(file.path(dir, "Makevars"), file.path(tdir, "Makevars"))
  cf <- grep(".c", dir(dir), value = TRUE, fixed = TRUE)
  cf <- cf[!grepl(".c~", cf, fixed = TRUE)]
  cf <- cf[!grepl("init.c", cf, fixed = TRUE)]
  for(j in cf) {
    file.copy(file.path(dir, j), file.path(tdir, j))
    system(paste("R CMD SHLIB", j))
    dyn.load(gsub(".c", ".so", j, fixed = TRUE))
  }
  setwd(hold)
}
