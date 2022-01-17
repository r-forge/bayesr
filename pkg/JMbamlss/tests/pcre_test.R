# New test file doing the following steps
# Create a PCA-Basis
# Evaluate PCA-Basis on random time points for different subjects
# Construct a smooth object
# Use a new data set to predict smooth


source("R/simMultiJM.R")
source("R/preprocessing.R")
source("R/eval_mfun.R")
source("R/pcre_smooth.R")

library(funData)
library(survival)
library(bamlss)
library(refund)

n <- 4
M <- 2
nmarker <- 1
nobs <- 5

mfpca <- create_true_MFPCA(M = M, nmarker = nmarker)
tpoints <- runif(nobs*n, min = 0, max = 120)
wfpc <- eval_mfpc(mfpca, tpoints)
dat <- data.frame(obstime = tpoints, id = factor(rep(seq_len(n), each = nobs)),
                  wfpc)

pcre <- smoothCon(s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = mfpca)),
                  dat)[[1]]
class(pcre) <- "pcre2.random.effect"
pcre$timevar <- "obstime"
pcre$term <- c(pcre$term, "obstime")

unc_pcre <- smoothCon(s(id, wfpc.1, wfpc.2, bs = "unc_pcre", 
                        xt = list("mfpc" = mfpca)), dat)[[1]]
unc_pcre$timevar <- "obstime"
unc_pcre$term <- c(unc_pcre$term, "obstime")

# alternatively 
# smoothCon(s(id, wfpc.1, wfpc.2, bs = "unc_pcre"), dat, scale.penalty = FALSE)
unc_pcre$S.scale*unc_pcre$S[[1]]

plot(dat$obstime, rowSums(unc_pcre$X[, c(2, 4, 6, 8)]) / sqrt(mfpca$values[2]))
plot(mfpca$functions)

newdat <- dat[rep(2*nobs, 11), ]

newdat$obstime <- seq(0, 120, length.out = 11)

pcre_pred <- PredictMat(pcre, newdat)
unc_pcre_pred <- PredictMat(unc_pcre, newdat)
plot(newdat$obstime, rowSums(pcre_pred[, c(2, 4, 6)]))

all.equal(unc_pcre_pred[1, c(3, 4)], 
          mfpca$functions@.Data[[1]]@X[, 1] * sqrt(mfpca$values))
