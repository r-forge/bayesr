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
                    "workstation" = paste0("/run/user/1000/gvfs/smb-share:",
                                           "server=clapton.wiwi.hu-berlin.de,",
                                           "share=volkmana.hub/JMbamlss/",
                                           "simulation/"),
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


# Load data
i <- 1
set.seed(i)
setting <- "scen_mauff/"
load(paste0(server_wd, setting, "/data/data", i, ".Rdata"))

simdat <- dat %>% 
  pivot_longer(y1:y6, names_to = "marker", values_to = "y") %>%
  mutate(marker = factor(marker, labels = paste0("m", 1:6)),
         id = factor(id))

ggplot(simdat,
       aes(x = year, y = y, colour = id)) + 
  facet_wrap(~marker, scales = "free_y") +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") + 
  ggtitle("Longitudinal Trajectories",
          "Simulated Data Set 1")

plot(survfit(Surv(Time, event) ~ group, data = dat.id))



# Univariate JM -----------------------------------------------------------


f_uni <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year + s(id, bs = "re") + s(year, id, bs = "re"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

b_uni1 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m1"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE) 
b_uni2 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m2"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni3 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m3"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE) 
b_uni4 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m4"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni5 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m5"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE)
b_uni6 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m6"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE) 

save(b_uni1, b_uni2, b_uni3, b_uni4, b_uni5, b_uni6, 
     file = paste0(server_wd, setting, "/uni/uni_", i, ".Rdata"))


summary(b_uni3$samples[[1]][,grep("\\.alpha", colnames(b_uni3$samples[[1]]))])
summary(b_uni4$samples[[1]][,grep("\\.alpha", colnames(b_uni4$samples[[1]]))])
summary(b_uni5$samples[[1]][,grep("\\.alpha", colnames(b_uni5$samples[[1]]))])





# Step width adaption -----------------------------------------------------

set.seed(i)
b_nu_3_all <- bamlss(f_uni, family = "jm", data = simdat %>%
                       filter(marker == "m3"),
                     timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                     burnin = 500, thin = 5, verbose = TRUE, update.nu = TRUE)
summary(b_nu_3_all$samples[[1]][,grep("\\.alpha", 
                                      colnames(b_nu_3_all$samples[[1]]))])

set.seed(i)
b_nu_5_all <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m5"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE, update.nu = TRUE)
summary(b_nu_5_all$samples[[1]][,grep("\\.alpha", 
                                      colnames(b_nu_5_all$samples[[1]]))])


set.seed(i)
b_nu_5 <- bamlss(f_uni, family = "jm", data = simdat %>%
                   filter(marker == "m5"),
                 timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
                 burnin = 500, thin = 5, verbose = TRUE, update.nu = TRUE,
                 start = parameters(b_uni5), optimizer = FALSE)
summary(b_nu_5$samples[[1]][,grep("\\.alpha", colnames(b_nu_5$samples[[1]]))])
save(b_nu_5, file = paste0(server_wd, setting, "/uni/nu5_", i, ".Rdata"))


# Functional Random Intercepts --------------------------------------------

f_fri <- list(
  Surv2(Time, event, obs = y) ~ -1 + s(Time, k = 20, bs = "ps"),
  gamma ~ 1 + group,
  mu ~  year + ti(id, bs = "re") + 
    ti(year, id, bs = c("cr", "re"), k = c(8, nlevels(simdat$id))),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)

set.seed(i)
# b_fri5 <- bamlss(f_fri, family = "jm", data = simdat %>%
#                    filter(marker == "m5"),
#                  timevar = "year", idvar = "id", maxit = 1500, n.iter = 5500,
#                  burnin = 500, thin = 5, verbose = TRUE)
# Takes too long to calculate!
#save(b_fri5, file = paste0(server_wd, setting, "/uni/uni_fri", i, ".Rdata"))

# Debugging ---------------------------------------------------------------


debug(bamlss:::sam_JM)
# Change nu <- 1 to nu <- 0.5 in sampler function
b_0.5 <- bamlss(f_uni, family = "jm", data = simdat %>% filter(marker == "m5"),
                timevar = "year", idvar = "id", optimizer = FALSE, 
                n.iter = 5500, burnin = 500, thin = 5, verbose = TRUE, 
                start = parameters(b_uni5))
summary(b_0.5$samples[[1]][,grep("\\.alpha", colnames(b_0.5$samples[[1]]))])



# Only Longitudinal Model -------------------------------------------------

f_long <- list(
  y ~  year + s(id, bs = "re") + s(year, id, bs = "re"),
  sigma ~ 1
)

set.seed(i)
b <- bamlss(f_long, data = simdat %>% filter(marker == "m5"),
            family = "gaussian")
summary(b$samples[[1]][,grep("\\.alpha", colnames(b$samples[[1]]))])



# Compare nu=1 and nu=0.5 -------------------------------------------------

set.seed(i)
b_nu05_s <- bamlss(f_uni, family = "jm",
                   data = simdat %>% filter(marker == "m5"),
                   timevar = "year", idvar = "id", optimizer = FALSE, 
                   n.iter = 10, burnin = 0, thin = 1, verbose = TRUE, 
                   start = parameters(b_uni5))
summary(b_nu05_s$samples[[1]][,grep("\\.alpha", colnames(b_nu05_s$samples[[1]]))])

set.seed(i)
b_nu1_s <- bamlss(f_uni, family = "jm",
                   data = simdat %>% filter(marker == "m5"),
                   timevar = "year", idvar = "id", optimizer = FALSE, 
                   n.iter = 10, burnin = 0, thin = 1, verbose = TRUE, 
                   start = parameters(b_uni5))
summary(b_nu1_s$samples[[1]][,grep("\\.alpha", colnames(b_nu1_s$samples[[1]]))])



# Step by Step Calculation ------------------------------------------------

set.seed(i)
debug(bamlss:::propose_jm_sigma)
b_nu05_s <- bamlss(f_uni, family = "jm",
                   data = simdat %>% filter(marker == "m5"),
                   timevar = "year", idvar = "id", optimizer = FALSE, 
                   n.iter = 10, burnin = 0, thin = 1, verbose = TRUE, 
                   start = parameters(b_uni5))
# In each MCMC step, do 
# before first Hessian calculation
# xgrad1 <- xgrad
# and
# assign("step2", list("xgrad1" = xgrad1, "xgrad2" = xgrad,
#                      "Sigma1" = Sigma, "Sigma2" = Sigma2, 
#                      "g1" = g0, "g2" = g, 
#                      "mu1" = mu, "mu2" = mu2, 
#                      "pibetaprop" = pibetaprop, "pibeta" = pibeta,
#                      "qbetaprop" = qbetaprop, "qbeta" = qbeta, 
#                      "p1" = p1, "p2" = p2),
#        envir = .GlobalEnv)
save(list = paste0("step", 1:10), 
     file = paste0(server_wd, setting, "/uni/steps_", i, ".Rdata"))


# Plot Proposal densities
xseq <- seq(-0.6, -0.3, length.out = 201)

# Step1
yseq_old <- dnorm(xseq, mean = step1$mu1, sd = sqrt(drop(step1$Sigma1)),
                  log = TRUE)
yseq_prop <- dnorm(xseq, mean = step1$mu2, sd = sqrt(drop(step1$Sigma2)),
                   log = TRUE)
plot(xseq, yseq_old, type = "l", ylab = "Log(Dnorm)", xlab = "eta_sigma",
     col = "blue", main = "MCMC Step 1")
abline(v = step1$mu1, col = "blue")
abline(v = step1$g1)
lines(xseq, yseq_prop, col = "red")
abline(v = step1$mu2, col = "red")
abline(v = step1$g2, col = "darkgreen", lty = "dotted")
legend("topleft",
       legend = c("eta0", "etaprop", "mu1", "mu2"),
       col = c("black", "darkgreen", "blue", "red"),
       lty = c("solid", "dotted", "solid", "solid"))

# Step2
yseq_old <- dnorm(xseq, mean = step2$mu1, sd = sqrt(drop(step2$Sigma1)),
                  log = TRUE)
yseq_prop <- dnorm(xseq, mean = step2$mu2, sd = sqrt(drop(step2$Sigma2)),
                   log = TRUE)
plot(xseq, yseq_old, type = "l", ylab = "Log(Dnorm)", xlab = "eta_sigma",
     col = "blue", main = "MCMC Step 2")
abline(v = step2$mu1, col = "blue")
abline(v = step2$g1)
lines(xseq, yseq_prop, col = "red")
abline(v = step2$mu2, col = "red")
abline(v = step2$g2, col = "darkgreen", lty = "dotted")
legend("topleft",
       legend = c("eta0", "etaprop", "mu1", "mu2"),
       col = c("black", "darkgreen", "blue", "red"),
       lty = c("solid", "dotted", "solid", "solid"))

# Step3
yseq_old <- dnorm(xseq, mean = step3$mu1, sd = sqrt(drop(step3$Sigma1)),
                  log = TRUE)
yseq_prop <- dnorm(xseq, mean = step3$mu2, sd = sqrt(drop(step3$Sigma2)),
                   log = TRUE)
plot(xseq, yseq_old, type = "l", ylab = "Log(Dnorm)", xlab = "eta_sigma",
     col = "blue", main = "MCMC Step 3")
abline(v = step3$mu1, col = "blue")
abline(v = step3$g1)
lines(xseq, yseq_prop, col = "red")
abline(v = step3$mu2, col = "red")
abline(v = step3$g2, col = "darkgreen", lty = "dotted")
legend("topright",
       legend = c("eta0", "etaprop", "mu1", "mu2"),
       col = c("black", "darkgreen", "blue", "red"),
       lty = c("solid", "dotted", "solid", "solid"))

# Step4
xseq <- seq(-1, -0.3, length.out = 201)
yseq_old <- dnorm(xseq, mean = step4$mu1, sd = sqrt(drop(step4$Sigma1)),
                  log = TRUE)
yseq_prop <- dnorm(xseq, mean = step4$mu2, sd = sqrt(drop(step4$Sigma2)),
                   log = TRUE)
plot(xseq, yseq_old, type = "l", ylab = "Log(Dnorm)", xlab = "eta_sigma",
     col = "blue", main = "MCMC Step 4")
abline(v = step4$mu1, col = "blue")
abline(v = step4$g1)
lines(xseq, yseq_prop, col = "red")
abline(v = step4$mu2, col = "red")
abline(v = step4$g2, col = "darkgreen", lty = "dotted")
legend("topright",
       legend = c("eta0", "etaprop", "mu1", "mu2"),
       col = c("black", "darkgreen", "blue", "red"),
       lty = c("solid", "dotted", "solid", "solid"))

# Gradient is extreme after the second step
sapply(paste0("step", 1:10), function (x) get(x)$xgrad1)
sapply(paste0("step", 1:10), function (x) get(x)$xgrad2)
sapply(paste0("step", 1:10), function (x) get(x)$Sigma1)
sapply(paste0("step", 1:10), function (x) get(x)$Sigma2)
sapply(paste0("step", 1:10), function (x) get(x)$g1)
sapply(paste0("step", 1:10), function (x) get(x)$g2)
sapply(paste0("step", 1:10), function (x) get(x)$pibeta)
sapply(paste0("step", 1:10), function (x) get(x)$pibetaprop)
