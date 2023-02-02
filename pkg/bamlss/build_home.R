library("bamlss")

if(!file.exists("../../www/misc")) {
  dir.create("../../www/misc")
}
figs <- dir("vignettes/figures", full.names = TRUE)
file.copy(figs, file.path("../../www/misc/", basename(figs)), overwrite = TRUE)

png("../../www/misc/mcycle1.png", units = "in", res = 150, width = 9, height = 3.5)
data("mcycle", package = "MASS")
load("vignettes/figures/toymodel.rda")
par(mar = c(4.1, 4.1, 0, 1.5), mfrow = c(1, 2), oma = c(0, 0, 2.5, 0))
p <- predict(b, model = "mu", FUN = c95)
plot2d(p ~ mcycle$times, fill.select = c(0, 1, 0, 1), scheme = 2, grid = 100,
  xlab = "times", ylab = expression(paste("Effect on ", mu)))
points(mcycle, cex = 0.6)
plot(b, model = "sigma", term = "s(times)", spar = FALSE, scheme = 2, grid = 100,
  ylab = expression(paste("Effect on ", log(sigma))))
mtext(expression(paste(accel, "~ N(", mu, "=", f(times), ",", log(sigma), "=", f(times), ")")),
  side = 3, line = 1, outer = TRUE, cex = 1.2)
dev.off()

png("../../www/misc/mcycle2.png", units = "in", res = 150, width = 8.5, height = 3.5)
par(mar = c(4.1, 4.1, 2.5, 1.1), mfrow = c(1, 2))
plot(b, which = 3, spar = FALSE, col = "lightgray")
plot(b, which = 4, spar = FALSE)
dev.off()

rmarkdown::render("index.Rmd", output_format = "md_document")

if(file.exists("index_files")) {
  ipath <- "index_files/figure-markdown_strict"
  figs <- dir(ipath, full.names = TRUE)
  bn <- paste0("bamlssindex-", basename(figs))
  nf <- file.path("misc", bn)

  file.copy(figs, nf, overwrite = TRUE)
  system("rm -rf index_files")

  md <- readLines("index.md")
  for(j in seq_along(figs))
    md <- gsub(figs[j], nf[j], md, fixed = TRUE)
  writeLines(md, "index.md")
}

pkgdown::build_home()

system("rm index.md")

