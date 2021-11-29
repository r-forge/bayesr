
# Comparison between data examples ----------------------------------------

# Data Generation
source("R/simMultiJM.R")

# Helperfunction PCRE
source("R/eval_mfun.R")
source("R/pcre_smooth.R")

# Family Construction
source("R/mjm_bamlss.R")
source("R/MJM_transform.R")
source("R/opt_MJM.R")
source("R/opt_updating.R")


# Data Generation ---------------------------------------------------------

library("refund")
library("bamlss")
library("MFPCA")

Surv2 <- bamlss:::Surv2

set.seed(1808)
d_new <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121), 
                lambda = function(t, x) {
                  1.4*log((120*t + 10)/1000)
                })

marker_dat <- split(d_new, d_new$marker)
marker_dat <- lapply(marker_dat, function (mark) {
  mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
  mark
})
m_irregFunData <- lapply(marker_dat, function (mark) {
  mark <- mark[order(mark$obstime), ]
  irregFunData(argvals = split(mark$obstime, mark$id), 
               X = split(mark$res, mark$id))
})
FPCA <- lapply(m_irregFunData, function(mark) {
  PACE(mark)
})
mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
uniExpansions <- lapply(FPCA, function (mark) {
  list(type = "given", functions = mark$functions)
})
MFPCA_new <- MFPCA(mFData = mFData, M = 2, uniExpansions = uniExpansions)



# PCRE Model --------------------------------------------------------------

f_new <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA_new)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_new <- bamlss(f_new, family = mjm_bamlss, data = d_new, timevar = "obstime",
                sampler = FALSE)


par(mfrow = c(1,2))
curve(1.4*log((120*x + 10)/1000), from = 0, to = 1)

pred_data <- data.frame(seq(0, 1, length.out = 100))
colnames(pred_data) <- b_new$x$lambda$smooth.construct[[1]]$term

k <- b_new$x$lambda$smooth.construct[[1]]$bs.dim
b_it <- b_new$x$lambda$smooth.construct[[1]]$state$parameters[-k]

pred_mat <- PredictMat(b_new$x$lambda$smooth.construct[[1]], pred_data)
plot(pred_data[, 1], pred_mat%*%b_new$parameters$lambda$s$`s(survtime)`[1:9])



# Old example -------------------------------------------------------------


d_old <- simMultiJM(nsub = 50)

marker_dat <- split(d_old, d_old$marker)
marker_dat <- lapply(marker_dat, function (mark) {
  mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
  mark
})
m_irregFunData <- lapply(marker_dat, function (mark) {
  mark <- mark[order(mark$obstime), ]
  irregFunData(argvals = split(mark$obstime, mark$id), 
               X = split(mark$res, mark$id))
})
FPCA <- lapply(m_irregFunData, function(mark) {
  PACE(mark)
})
mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
uniExpansions <- lapply(FPCA, function (mark) {
  list(type = "given", functions = mark$functions)
})
MFPCA_old <- MFPCA(mFData = mFData, M = 2, uniExpansions = uniExpansions)

f_old <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA_old)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)

# Model fit with PCRE and plot
b_old <- bamlss(f_old, family = mjm_bamlss, data = d_old, timevar = "obstime",
                sampler = FALSE)

curve(1.4*log((x + 10)/1000), from = 0, to = 120)

pred_data <- data.frame(seq(0, 120, length.out = 100))
colnames(pred_data) <- b_old$x$lambda$smooth.construct[[1]]$term

k <- b_old$x$lambda$smooth.construct[[1]]$bs.dim
b_it <- b_old$x$lambda$smooth.construct[[1]]$state$parameters[-k]

pred_mat <- PredictMat(b_old$x$lambda$smooth.construct[[1]], pred_data)
plot(pred_data[, 1], pred_mat%*%b_old$parameters$lambda$s$`s(survtime)`[1:9])


# Model fit with FRE and plot
f_re <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + ti(obstime, by = marker) +
    ti(id, bs = "re") +
    ti(id, obstime, bs = c("re", "cr"), k = c(50, 5)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_re <- bamlss(f_re, family = mjm_bamlss, data = d_old, timevar = "obstime",
               sampler = FALSE)

curve(1.4*log((x + 10)/1000), from = 0, to = 120)
plot(pred_data[, 1], pred_mat%*%b_re$parameters$lambda$s$`s(survtime)`[1:9])



# Old data set with first reducing time scale -----------------------------


d_old_1 <- d_old
d_old_1$survtime <- d_old_1$survtime/120
marker_dat <- split(d_old_1, d_old_1$marker)
marker_dat <- lapply(marker_dat, function (mark) {
  mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
  mark
})
m_irregFunData <- lapply(marker_dat, function (mark) {
  mark <- mark[order(mark$obstime), ]
  irregFunData(argvals = split(mark$obstime, mark$id), 
               X = split(mark$res, mark$id))
})
FPCA <- lapply(m_irregFunData, function(mark) {
  PACE(mark)
})
mFData <- multiFunData(lapply(FPCA, "[[", "fit"))
uniExpansions <- lapply(FPCA, function (mark) {
  list(type = "given", functions = mark$functions)
})
MFPCA_old_1 <- MFPCA(mFData = mFData, M = 2, uniExpansions = uniExpansions)

# List element has to be updated
f_old_1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA_old_1)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)

# Model fit with pcre and plot
b_old_1 <- bamlss(f_old_1, family = mjm_bamlss, data = d_old_1, timevar = "obstime",
                sampler = FALSE)

curve(1.4*log((x*120 + 10)/1000), from = 0, to = 1)

pred_data <- data.frame(seq(0, 1, length.out = 100))
colnames(pred_data) <- b_old_1$x$lambda$smooth.construct[[1]]$term

k <- b_old_1$x$lambda$smooth.construct[[1]]$bs.dim
b_it <- b_old_1$x$lambda$smooth.construct[[1]]$state$parameters[-k]

pred_mat <- PredictMat(b_old_1$x$lambda$smooth.construct[[1]], pred_data)
plot(pred_data[, 1], pred_mat%*%b_old_1$parameters$lambda$s$`s(survtime)`[1:9])


# Model fit with fre and plot
b_re_1 <- bamlss(f_re, family = mjm_bamlss, data = d_old_1, timevar = "obstime",
                 sampler = FALSE)

curve(1.4*log((x*120 + 10)/1000), from = 0, to = 1)
plot(pred_data[, 1], pred_mat%*%b_re_1$parameters$lambda$s$`s(survtime)`[1:9])
