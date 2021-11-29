
# Multivariate Joint Model Example ----------------------------------------

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

if(!exists("d")) {
  
  library("refund")
  library("bamlss")
  library("MFPCA")
  
  Surv2 <- bamlss:::Surv2
  
  set.seed(1808)
  d <- simMultiJM(nsub = 50, times = seq(0, 1, length.out = 121), 
                  lambda = function(t, x) {
                    1.4*log((120*t + 10)/1000)
                  })
  
  marker_dat <- split(d, d$marker)
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
  MFPCA <- MFPCA(mFData = mFData, M = 2, uniExpansions = uniExpansions)
  
}




# PCRE Model --------------------------------------------------------------

f <- list(
  Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
  gamma ~ 1,
  mu ~ -1 + marker + s(obstime, by = marker) +
    s(id, wfpc.1, wfpc.2, bs = "pcre", xt = list("mfpc" = MFPCA)),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker + s(survtime, by = marker)
)
b <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime",
            sampler = FALSE)



# RE Model ----------------------------------------------------------------

# f_re <- list(
#   Surv2(survtime, event, obs = y) ~ -1 + s(survtime),
#   gamma ~ 1, 
#   mu ~ -1 + marker + ti(obstime, by = marker) +
#     ti(id, bs = "re") +
#     ti(id, obstime, bs = c("re", "cr"), k = c(50, 5)),
#   sigma ~ -1 + marker,
#   alpha ~ -1 + marker + s(survtime, by = marker)
# )
# b <- bamlss(f_re, family = mjm_bamlss, data = d, timevar = "obstime")


