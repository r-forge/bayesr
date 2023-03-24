# Code from MAUFF et al. (2021)

# Specify location
location <- "workstation"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if(location == "server_linux") "./simulation/"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
}
setting <- "scen_mauff"

#-----------------------------------------------------------------------------

# Creating data for simulation scenario E, (joint longitudinal and survival model, 6 longitudinal outcomes) ###################

# loading libraries
library(MASS)
library(JMbayes)
library(Matrix)

set.seed(123)
# creating datasets
random.seed <- sample(c(5001:6000), 500, replace = FALSE)

for (datanum in c(1:200)) {
  
  random_seed_inuse <- random.seed[datanum]
  set.seed(random_seed_inuse)
  
  n <- 500
  K <- 15
  t.max <- 25
  
  # parameter values (longitudinal)  
  betas1 <- c(4.93 , 0.4) 
  betas2 <- c(3.58 , 1.0)
  betas3 <- c(1.46 , -0.70)
  betas4 <- c(1.78 , 0.5)
  betas5 <- c(0.31 , -1.2)
  betas6 <- c(1.71 , -0.3)
  
  sigma.y1 <- 0.86
  sigma.y2 <- 0.39
  sigma.y3 <- 0.94 
  sigma.y4 <- 0.40
  sigma.y5 <- 0.36 
  sigma.y6 <- 0.57 
  
  # parameter values (survival)
  gammas <- c(-5.8, 0.5)
  alpha1 <- 0.10 
  alpha2 <- -0.6 
  alpha3 <- -0.1   
  alpha4 <- -1.41   
  alpha5 <- -1.81 
  alpha6 <- 0.75  
  
  phi <- 1.65
  mean.Cens <- 15
  
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
  
  ################################################
  
  times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max)))))
  group <- sample(rep(0:1, each = n/2))
  DF <- data.frame(year = times, group = factor(rep(group, each = K)))
  X1 <- model.matrix(~ year , data = DF)
  Z1 <- model.matrix(~ year, data = DF)
  X2 <- model.matrix(~ year, data = DF)
  Z2 <- model.matrix(~ year, data = DF)
  X3 <- model.matrix(~ year, data = DF)
  Z3 <- model.matrix(~ year, data = DF)
  
  X4 <- model.matrix(~ year , data = DF)
  Z4 <- model.matrix(~ year, data = DF)
  X5 <- model.matrix(~ year, data = DF)
  Z5 <- model.matrix(~ year, data = DF)
  X6 <- model.matrix(~ year, data = DF)
  Z6 <- model.matrix(~ year, data = DF)
  
  W <- cbind(1, group)
  
  ################################################
  
  b <- mvrnorm(n, rep(0, nrow(D)), D)
  
  id <- rep(1:n, each = K)
  eta.y1 <- as.vector(X1 %*% betas1 + rowSums(Z1 * b[id, 1:2])) # linear predictor
  eta.y2 <- as.vector(X2 %*% betas2 + rowSums(Z2 * b[id, 3:4])) # linear predictor
  eta.y3 <- as.vector(X3 %*% betas3 + rowSums(Z3 * b[id, 5:6])) # linear predictor
  eta.y4 <- as.vector(X4 %*% betas4 + rowSums(Z4 * b[id, 7:8])) # linear predictor
  eta.y5 <- as.vector(X5 %*% betas5 + rowSums(Z5 * b[id, 9:10])) # linear predictor
  eta.y6 <- as.vector(X6 %*% betas6 + rowSums(Z6 * b[id, 11:12])) # linear predictor
  
  y1 <- rnorm(n * K, eta.y1, sigma.y1)
  y2 <- rnorm(n * K, eta.y2, sigma.y2)
  y3 <- rnorm(n * K, eta.y3, sigma.y3)
  y4 <- rnorm(n * K, eta.y4, sigma.y4)
  y5 <- rnorm(n * K, eta.y5, sigma.y5)
  y6 <- rnorm(n * K, eta.y6, sigma.y6)
  
  eta.t <- as.vector(W %*% gammas)
  invS <- function (t, u, i) {
    h <- function (s) {
      group.i <- group[i]
      XX1 <- cbind(1, s)
      ZZ1 <- cbind(1, s)
      XX2 <- cbind(1, s)
      ZZ2 <- cbind(1, s)
      XX3 <- cbind(1, s)
      ZZ3 <- cbind(1, s)
      XX4 <- cbind(1, s)
      ZZ4 <- cbind(1, s)
      XX5 <- cbind(1, s)
      ZZ5 <- cbind(1, s)
      XX6 <- cbind(1, s)
      ZZ6 <- cbind(1, s)
      f1 <- as.vector(XX1 %*% betas1 + rowSums(ZZ1 * b[rep(i, nrow(ZZ1)), 1:2]))
      f2 <- as.vector(XX2 %*% betas2 + rowSums(ZZ2 * b[rep(i, nrow(ZZ2)), 3:4]))
      f3 <- as.vector(XX3 %*% betas3 + rowSums(ZZ3 * b[rep(i, nrow(ZZ3)), 5:6]))
      f4 <- as.vector(XX4 %*% betas4 + rowSums(ZZ4 * b[rep(i, nrow(ZZ4)), 7:8]))
      f5 <- as.vector(XX5 %*% betas5 + rowSums(ZZ5 * b[rep(i, nrow(ZZ5)), 9:10]))
      f6 <- as.vector(XX6 %*% betas6 + rowSums(ZZ6 * b[rep(i, nrow(ZZ6)), 11:12]))
      exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha1 + f2 * alpha2 + f3 * alpha3 + f4*alpha4 +
            f5*alpha5 + f6*alpha6)
      
    }
    integrate(h, lower = 0, upper = t)$value + log(u)
  }
  u <- runif(n)
  trueTimes <- numeric(n)
  for (i in 1:n) {
    Up <- 50
    tries <- 5
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    while(inherits(Root, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 200
      Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
    }
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
  }
  trueTimes[is.na(trueTimes)] <- 1e06
  
  # simulate censoring times from an exponential distribution,
  # and calculate the observed event times, i.e., min(true event times, censoring times)
  Ctimes <- runif(n, 0, 2 * mean.Cens)
  Time <- pmin(trueTimes, Ctimes)
  event <- as.numeric(trueTimes <= Ctimes) # event indicator
  table(event)
  
  ################################################
  
  
  # keep the nonmissing cases, i.e., drop the longitudinal measurements
  # that were taken after the observed event time for each subject.
  ind <- times <= rep(Time, each = K)
  y1 <- y1[ind]
  y2 <- y2[ind]
  y3 <- y3[ind]
  y4 <- y4[ind]
  y5 <- y5[ind]
  y6 <- y6[ind]
  
  id <- id[ind]
  id <- match(id, unique(id))
  
  dat <- DF[ind, ]
  dat$id <- id
  dat$y1 <- y1
  dat$y2 <- y2
  dat$y3 <- y3
  dat$y4 <- y4
  dat$y5 <- y5
  dat$y6 <- y6
  
  dat$Time <- Time[id]
  dat$event <- event[id]
  
  dat <- dat[c("id", "year", "y1", "y2", "y3", 'y4', 'y5', 'y6', "Time", "event", "group")]
  
  dat.id <- dat[!duplicated(dat$id), ]
  
  summary(tapply(id, id, length))
  table(event)
  n
  mean(event)
  summary(Time)
  
  # delete all unused objects
  which_to_delete <- function () {
    objs <- ls(pos = 1L)
    objs[!objs %in% c("dat", "dat.id", "rootIn", "rootOut", "args", "nsub", 
                      "scenario", "datanum", "random.seed", "results_wd",
                      "setting")]
  }
  
  rm(list = which_to_delete())
  
  fname <- paste(results_wd, setting, '/data/data', datanum, '.RData', sep = "")
  save.image(fname)
  print(datanum)
  
}


