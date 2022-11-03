REPOS=http://svn.gwdg.de/svn/bayesx/trunk
USER=guest
PASSWD=guest
DIRS="adaptiv alex andrea bib dag graph leyre mcmc psplines samson structadd"
FILES="export_type.h main.cpp values.h"
mkdir -p src/bayesxsrc
cd src/bayesxsrc
for i in $DIRS ; do
  svn checkout --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
for i in $FILES ; do
  svn export --username "${USER}" --password "${PASSWD}" $REPOS/$i $i
done
cd ..
cp dev-Makefile Makefile
cp dev-Makefile.win Makefile.win
cp dev-Makefile.inner.win Makefile.inner.win

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
  file.copy(file.path("bayesx", d), file.path(owd, "src", "bayesxsrc", d), overwrite = TRUE, recursive = TRUE)
}
for(f in files) {
  file.copy(file.path("bayesx", f), file.path(owd, "src", "bayesxsrc", f), overwrite = TRUE, recursive = TRUE)
}

