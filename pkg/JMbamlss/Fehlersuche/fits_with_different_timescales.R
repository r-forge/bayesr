
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

# Load Packages
library("refund")
library("bamlss")
library("MFPCA")

Surv2 <- bamlss:::Surv2

par(mfrow = c(1,2))

set.seed(1808)



# PCRE Model for new example ----------------------------------------------


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


curve(1.4*log((120*x + 10)/1000), from = 0, to = 1)

pred_data <- data.frame(seq(0, 1, length.out = 100))
colnames(pred_data) <- b_new$x$lambda$smooth.construct[[1]]$term

k <- b_new$x$lambda$smooth.construct[[1]]$bs.dim
b_it <- b_new$x$lambda$smooth.construct[[1]]$state$parameters[-k]

pred_mat <- PredictMat(b_new$x$lambda$smooth.construct[[1]], pred_data)
plot(pred_data[, 1], pred_mat%*%b_new$parameters$lambda$s$`s(survtime)`[1:9])


# Change the time-scale ---------------------------------------------------

d_new120 <- d_new
d_new120$survtime <- d_new120$survtime*120
d_new120$obstime <- d_new120$obstime*120

marker_dat120 <- split(d_new120, d_new120$marker)
marker_dat120 <- lapply(marker_dat120, function (mark) {
  mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
  mark
})
m_irregFunData120 <- lapply(marker_dat120, function (mark) {
  mark <- mark[order(mark$obstime), ]
  irregFunData(argvals = split(mark$obstime, mark$id), 
               X = split(mark$res, mark$id))
})
FPCA120 <- lapply(m_irregFunData120, function(mark) {
  PACE(mark)
})
mFData120 <- multiFunData(lapply(FPCA120, "[[", "fit"))
uniExpansions120 <- lapply(FPCA120, function (mark) {
  list(type = "given", functions = mark$functions)
})
MFPCA_new120 <- MFPCA(mFData = mFData120, M = 2, 
                      uniExpansions = uniExpansions120)

f_new120 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA_new120)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b_new120 <- bamlss(f_new120, family = mjm_bamlss, data = d_new120, 
                   timevar = "obstime", sampler = FALSE)

pred_data120 <- data.frame(seq(0, 120, length.out = 100))
colnames(pred_data120) <- b_new120$x$lambda$smooth.construct[[1]]$term

k120 <- b_new120$x$lambda$smooth.construct[[1]]$bs.dim
b_it120 <- b_new120$x$lambda$smooth.construct[[1]]$state$parameters[-k120]

pred_mat120 <- PredictMat(b_new120$x$lambda$smooth.construct[[1]], pred_data120)
curve(1.4*log((x + 10)/1000), from = 0, to = 120)
plot(pred_data120[, 1],
     pred_mat120%*%b_new120$parameters$lambda$s$`s(survtime)`[1:9])



# Old example -------------------------------------------------------------


d_old <- simMultiJM(nsub = 50)

marker_dat_old <- split(d_old, d_old$marker)
marker_dat_old <- lapply(marker_dat_old, function (mark) {
  mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
  mark
})
m_irregFunData_old <- lapply(marker_dat_old, function (mark) {
  mark <- mark[order(mark$obstime), ]
  irregFunData(argvals = split(mark$obstime, mark$id), 
               X = split(mark$res, mark$id))
})
FPCA_old <- lapply(m_irregFunData_old, function(mark) {
  PACE(mark)
})
mFData_old <- multiFunData(lapply(FPCA_old, "[[", "fit"))
uniExpansions_old <- lapply(FPCA_old, function (mark) {
  list(type = "given", functions = mark$functions)
})
MFPCA_old <- MFPCA(mFData = mFData_old, M = 2, 
                   uniExpansions = uniExpansions_old)

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

pred_data_old <- data.frame(seq(0, 120, length.out = 100))
colnames(pred_data_old) <- b_old$x$lambda$smooth.construct[[1]]$term

k_old <- b_old$x$lambda$smooth.construct[[1]]$bs.dim
b_it_old <- b_old$x$lambda$smooth.construct[[1]]$state$parameters[-k_old]

pred_mat_old <- PredictMat(b_old$x$lambda$smooth.construct[[1]], pred_data_old)
plot(pred_data_old[, 1], 
     pred_mat_old%*%b_old$parameters$lambda$s$`s(survtime)`[1:9])


# Old data set with FREs on 120 -------------------------------------------


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
pred_data_re <- data.frame(seq(0, 120, length.out = 100))
colnames(pred_data_re) <- b_re$x$lambda$smooth.construct[[1]]$term
k_re <- b_re$x$lambda$smooth.construct[[1]]$bs.dim
b_it_re <- b_re$x$lambda$smooth.construct[[1]]$state$parameters[-k_re]
pred_mat_re <- PredictMat(b_re$x$lambda$smooth.construct[[1]], pred_data_re)
plot(pred_data_re[, 1],
     pred_mat_re%*%b_re$parameters$lambda$s$`s(survtime)`[1:9])



# Old data set with first reducing time scale -----------------------------


d_old1 <- d_old
d_old1$survtime <- d_old1$survtime/120
d_old1$obstime <- d_old1$obstime/120
marker_dat1 <- split(d_old1, d_old1$marker)
marker_dat1 <- lapply(marker_dat1, function (mark) {
  mark$res <- bam(y ~ s(obstime) + s(x2), data = mark)$residuals
  mark
})
m_irregFunData1 <- lapply(marker_dat1, function (mark) {
  mark <- mark[order(mark$obstime), ]
  irregFunData(argvals = split(mark$obstime, mark$id), 
               X = split(mark$res, mark$id))
})
FPCA1 <- lapply(m_irregFunData1, function(mark) {
  PACE(mark)
})
mFData1 <- multiFunData(lapply(FPCA1, "[[", "fit"))
uniExpansions1 <- lapply(FPCA1, function (mark) {
  list(type = "given", functions = mark$functions)
})
MFPCA_old1 <- MFPCA(mFData = mFData1, M = 2, uniExpansions = uniExpansions1)

# List element has to be updated
f_old1 <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA_old1)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)

# Model fit with pcre and plot
b_old1 <- bamlss(f_old1, family = mjm_bamlss, data = d_old1, timevar = "obstime",
                sampler = FALSE)

curve(1.4*log((x*120 + 10)/1000), from = 0, to = 1)

pred_data1 <- data.frame(seq(0, 1, length.out = 100))
colnames(pred_data1) <- b_old1$x$lambda$smooth.construct[[1]]$term

k1 <- b_old1$x$lambda$smooth.construct[[1]]$bs.dim
b_it1 <- b_old1$x$lambda$smooth.construct[[1]]$state$parameters[-k1]

pred_mat1 <- PredictMat(b_old1$x$lambda$smooth.construct[[1]], pred_data1)
plot(pred_data1[, 1], pred_mat1%*%b_old1$parameters$lambda$s$`s(survtime)`[1:9])



# Old data set with FREs on 1 ---------------------------------------------


# Model fit with fre and plot
b_re1 <- bamlss(f_re, family = mjm_bamlss, data = d_old1, timevar = "obstime",
                 sampler = FALSE)

curve(1.4*log((x*120 + 10)/1000), from = 0, to = 1)

pred_data1re <- data.frame(seq(0, 1, length.out = 100))
colnames(pred_data1re) <- b_re1$x$lambda$smooth.construct[[1]]$term
k1re <- b_re1$x$lambda$smooth.construct[[1]]$bs.dim
b_it1re <- b_re1$x$lambda$smooth.construct[[1]]$state$parameters[-k1re]
pred_mat1re <- PredictMat(b_re1$x$lambda$smooth.construct[[1]], pred_data1re)

plot(pred_data1re[, 1], pred_mat1re%*%b_re1$parameters$lambda$s$`s(survtime)`[1:9])
