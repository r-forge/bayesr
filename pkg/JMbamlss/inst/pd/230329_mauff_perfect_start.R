# Set up R session --------------------------------------------------------



# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}
server_wd <- switch(location,
  "workstation" = paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu",
                         "-berlin.de,share=volkmana.hub/JMbamlss/simulation/"),
  "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation/",
  "server_windows" = "H:/JMbamlss/simulation/")


# Always
library(survival)
library(JMbayes2)
library(bamlss)
library(MFPCA)
library(tidyverse)
library(parallel)
library(Rcpp)
library(Matrix)
library(sparseFLMM)
library(JMbamlss)

setting <- "scen_mauff"
i <- 1
set.seed(i)

# Load the data
load(paste0(server_wd, setting, "/data/data", i, ".Rdata"))


# Set the correct parameters ----------------------------------------------

# Get the true lambda parameters
phi <- 1.65
hazard_dat <- data.frame(
  id = dat.id$id,
  Time = dat.id$Time,
  hazard = log(phi) + (phi - 1) * log(dat.id$Time) - 5.8
)

m <- gam(hazard ~  s(Time, k = 20, bs = "ps"), data = hazard_dat)
lambda_pars <- c(m$coefficients[2:20], m$sp, sum(m$edf[2:20]))
names(lambda_pars) <- c(paste0("b", 1:19), "tau21", "edf")


# Get the true RE parameters
library(MASS)
set.seed(123)
random.seed <- sample(c(5001:6000), 500, replace = FALSE)
datanum <- 1
n <- 500
K <- 15
t.max <- 25
random_seed_inuse <- random.seed[datanum]
set.seed(random_seed_inuse)
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max)))))
group <- sample(rep(0:1, each = n/2))
D_num <- matrix(NA, nrow = 12, ncol = 12)
D_num[1:4, 1:4] <- matrix(c(0.919,  0.003,  0.220, -0.008,
                            0.003,  0.006,  0.005, -0.002,
                            0.220,  0.005,  0.421, -0.013,
                            -0.008, -0.002, -0.013,  0.008), byrow = TRUE)
D_num[5:8, 5:8] <- matrix(c(0.551,  0.007, -0.141,  0.015,
                            0.007,  0.014, -0.005, -0.005,
                            -0.141, -0.005,  0.204, -0.042,
                            0.015, -0.005, -0.042,  0.043), byrow = TRUE)  
D_num[9:12, 9:12] <- matrix(c(0.110, -0.028,  0.176,  0.021,
                              -0.028,  0.035, -0.026, -0.003,
                              0.176, -0.026,  3.580,  0.040,
                              0.021, -0.003,  0.040,  0.197), byrow = TRUE)  
D_num[1:4, 5:8] <- matrix(c(0.308,  0.005,  0.265, -0.013,
                            0.012,  0.001,  0.012,  0.002,
                            -0.073, -0.007, -0.172,  0.018,
                            0.007,  0.005,  0.012, -0.013), byrow = TRUE)
D_num[1:4, 9:12] <- matrix(c(0.049,  0.005,  0.124, -0.013,
                             -0.007, -0.005, -0.010,  0.012,
                             0.698,  0.006,  0.680, -0.027,
                             0.056,  0.006,  0.034,  0.004), byrow = TRUE)  
D_num[5:8, 9:12] <- matrix(c(0.095,  0.004, -0.144,  0.032,
                             -0.013,  0.005,  0.035, -0.037,
                             0.826,  0.018, -0.263,  0.024,
                             0.077,  0.019, -0.015,  0.006), byrow = TRUE)
D_num[5:8, 1:4] <- t(D_num[1:4, 5:8])
D_num[9:12, 1:4] <- t(D_num[1:4, 9:12])
D_num[9:12, 5:8] <- t(D_num[5:8, 9:12])
D_numsym <- forceSymmetric(D_num) 
D <- D_numsym
b <- mvrnorm(n, rep(0, nrow(D)), D)
seq1 <- seq(0, max(dat.id$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)
mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs, scores = b)
# names(b1) <- c(paste0("b", 1:500), "tau21", "edf")

for (l in 1:12) {
  assign(paste0("b", l), c(mfpca_tru$scores[, l], 
                           "tau21" = 1/mfpca_tru$values[l], 
                           "edf" = 500))
}

# Full parameter list
par <- list("lambda" = list("s(Time)" = lambda_pars),
            "gamma" = list("p" = c(m$coefficients[1], "group1" = 0.5)),
            "mu" = list("s" = list("s(id,fpc.1)" = b1,
                                   "s(id,fpc.2)" = b2,
                                   "s(id,fpc.3)" = b3,
                                   "s(id,fpc.4)" = b4,
                                   "s(id,fpc.5)" = b5,
                                   "s(id,fpc.6)" = b6,
                                   "s(id,fpc.7)" = b7,
                                   "s(id,fpc.8)" = b8,
                                   "s(id,fpc.9)" = b9,
                                   "s(id,fpc.10)" = b10,
                                   "s(id,fpc.11)" = b11,
                                   "s(id,fpc.12)" = b12),
                        "p" = c("markerm1" = 4.93,
                                "markerm2" = 3.58,
                                "markerm3" = 1.46,
                                "markerm4" = 1.78,
                                "markerm5" = 0.31,
                                "markerm6" = 1.71,
                                "markerm1:year" = 0.4,
                                "markerm2:year" = 1,
                                "markerm3:year" = -0.7,
                                "markerm4:year" = 0.5,
                                "markerm5:year" = -1.2,
                                "markerm6:year" = -0.3)),
            "sigma" = c("markerm1" = log(0.86),
                        "markerm2" = log(0.39),
                        "markerm3" = log(0.94),
                        "markerm4" = log(0.4),
                        "markerm5" = log(0.36),
                        "markerm6" = log(0.57)),
            "alpha" = c("markerm1" = 0.1,
                        "markerm2" = -0.6,
                        "markerm3" = -0.1,
                        "markerm4" = -1.41,
                        "markerm5" = -1.81,
                        "markerm6" = 0.75)
)



# Estimate Model with True Parameters -------------------------------------

# Load the data
load(paste0(server_wd, setting, "/data/data", i, ".Rdata"))

d_rirs <- pivot_longer(dat, y1:y6, names_to = "marker", values_to = "y") %>%
  mutate(marker = factor(marker, labels = paste0("m", 1:6)),
         id = factor(id)) %>%
  arrange(marker, id, year) %>%
  as.data.frame()

mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })
d_rirs_tru <- JMbamlss:::attach_wfpc(mfpca_tru, d_rirs, n = 12,
                                     obstime = "year")

nfpc <- 12
f_tru <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  as.formula(paste0(
    "mu ~ -1 + marker + year:marker +",
    paste0(lapply(seq_len(nfpc), function(x) {
      paste0("s(id, fpc.", x, ", bs = 'unc_pcre', xt = list('mfpc' = ",
             "mfpca_list[[", x, "]]))")
    }), collapse = " + "))),
  sigma ~ -1 + marker,
  alpha ~ -1 + marker
)

set.seed(1)
b_tru <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = d_rirs_tru, 
                timevar = "year", optimizer = FALSE,# n.iter = 1500,
                #burnin = 500, thin = 5,
                verbose  = TRUE, start = par)



# Estimate the univariate model with true parameters ----------------------

# Full parameter list
par <- list("lambda" = list("s(Time)" = lambda_pars),
            "gamma" = list("p" = c(m$coefficients[1], "group1" = 0.5)),
            "mu" = list("s" = list("s(id)" = b[, 1],
                                   "s(year,id)" = b[, 2]),
                        "p" = c("(Intercept)" = 4.93,
                                "year" = 3.58)),
            "sigma" = list(c("(Intercept)" = log(0.86))),
            "alpha" = list(c("(Intercept)" = 0.1))
)
f_uni <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year +
    s(id, bs = "re") +
    s(year, id, bs ="re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(1)
b_uni <- bamlss(f_uni, family = jm_bamlss, data = d_rirs_tru %>%
                  filter(marker == "m1") %>% droplevels() %>% as.data.frame(), 
                timevar = "year", optimizer = FALSE,
                verbose  = TRUE, start = par, idvar = "id")

