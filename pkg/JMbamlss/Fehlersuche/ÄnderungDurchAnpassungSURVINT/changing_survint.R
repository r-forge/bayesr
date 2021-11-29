
# Try out model fit with different versions survint_gq --------------------

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

curve(1.4*log((120*x + 10)/1000), from = 0, to = 1)


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
it_param_new <- it_param


# Old Version survint_gq --------------------------------------------------

survint_gq <- function(pred = c("lambda", "gamma", "long"), pre_fac,
                       pre_vec = NULL, omega, int_fac = NULL,
                       int_vec = NULL, weights, survtime) {
  
  if (sum(c(is.null(pre_vec), is.null(int_vec))) != 1) {
    stop("Either pre_vec or int_vec must be specified.")
  }
  if (pred != "long" & !is.null(int_fac)) {
    stop("Argument int_fac is only used for longitudinal predictors.")
  }
  # int_vec <- (t(rep(1, nmarker)) %x% diag(length(eta_timegrid))) %*%
  #   (eta_timegrid_alpha * x$Xgrid)
  # Integration as weighted sum of evaluated points
  n <- length(survtime)
  gq_mat <- diag(n)%x%t(weights)
  
  if (pred == "long") {
    nmarker <- nrow(int_vec)/(n*length(weights))
    if (nmarker %% 1 != 0) {
      stop("Dimensions of longitudinal design matrix do not match.")
    }
    pre_fac <- rep(pre_fac, nmarker)
    omega <- rep(omega, nmarker)
    gq_mat <- diag(n*nmarker)%x%t(weights)
    
    # Alternativ mit Matrixmultiplikation statt verlängertem Vektor
    # dim_mat <- t(rep(1, nmarker)) %x% diag(n*length(weights))
    # score_int <- pre_fac*gq_mat %*% (omega*dim_mat %*% (int_fac*int_vec))
    # hess_int <- pre_fac*gq_mat %*% 
    #   (omega*dim_mat %*% t(apply(int_fac*int_vec, 1, tcrossprod)))
  }
  switch(pred, 
         "lambda" =, "long" = {
           if (is.null(int_fac)) {
             int_fac <- rep(1, n*length(weights))
           }
           score_int <- pre_fac*gq_mat %*% (omega*int_fac*int_vec)
           hess_int <- pre_fac*gq_mat %*% (omega*int_fac^2*t(apply(int_vec, 1,
                                                                   tcrossprod)))
         },
         "gamma" = {
           pre_fac <- c(pre_fac*gq_mat %*% omega)
           score_int <- pre_fac*pre_vec
           hess_int <- pre_fac*t(apply(pre_vec, 1, tcrossprod))
         })
  
  if (dim(score_int)[1] != dim(hess_int)[1]) {
    hess_int <- t(hess_int)
    if (dim(score_int)[1] != dim(hess_int)[1]) {
      stop("Problem with dimensions in gauss quadrature.")
    }
  }
  # Multiply with borders of integration for Legendre
  list(score_int = survtime/2 * score_int, hess_int = survtime/2 * hess_int)
  # hess_int: each row corresponds to one individual 
  # each row has ncol(vec)^2 elements -> is a matrix of derivatives
}


# New function call -------------------------------------------------------


b_old <- bamlss(f, family = mjm_bamlss, data = d, timevar = "obstime",
                sampler = FALSE)

all.equal(b, b_old)
all.equal(it_param_new, it_param)
identical(b$parameters, b_old$parameters)
