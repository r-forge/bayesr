#  Univariate Example from Meike ------------------------------------------

library(bamlss)
library(tidyverse)
library(RColorBrewer)
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=stats_thesis/output/mcmc_b_functional_constant_zero_300_17",
            "7.RData"))
dat <- b$model.frame
dat_id <- dat %>% group_by(id) %>% filter(row_number() == n()) %>% ungroup() %>%
  as.data.frame
times <- seq(0, 72, 1)
newdat <- dat_id[rep(1:300, each = length(times)), ]
newdat$obstime <- rep(times, 300)
newdat <- newdat[newdat$survtime > newdat$obstime, ]
newdat$mu <- predict(b, newdat, model = "mu")
newdat$fri <- predict(b, newdat, model = "mu", term = c(3, 5))
newdat$mean <- newdat$mu - newdat$fri

newdat <- left_join(newdat,
                    dat %>% mutate(y = b$y[[1]][, 3]) %>% 
                      select("id", "obstime", "y"),
                    by = c("id", "obstime"))
set.seed(81891)
sampl <- sample(seq_len(300), 5)
p <- ggplot(data = newdat %>% filter(id %in% sampl), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mu), size = 1.2) +
  geom_point(aes(y = y, shape = id), size = 3) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  ggtitle("Functional RI")
pdf("../Fig1_2.pdf", width = 6, height = 6)
p
dev.off()

p <- ggplot(data = newdat %>% filter(id %in% sampl), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fri), size = 1.2) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time")
pdf("../Fig2.pdf", width = 6, height = 6)
p
dev.off()


p <- ggplot(data = newdat %>% filter(id %in% sampl), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mu), size = 1.2) +
  geom_point(aes(y = y, shape = id), size = 3) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  ggtitle("Data and Fit")
pdf("../riga_Fig2_1.pdf", width = 6, height = 6)
p
dev.off()

p <- ggplot(data = newdat %>% filter(id %in% sampl), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mean), size = 1.2) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  ggtitle("Fixed Effects")
pdf("../riga_Fig2_2.pdf", width = 6, height = 6)
p
dev.off()

p <- ggplot(data = newdat %>% filter(id %in% sampl), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fri), size = 1.2) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  ggtitle("Functional Random Intercept")
pdf("../riga_Fig2_3.pdf", width = 6, height = 6)
p
dev.off()

f <- list(
  Surv2(survtime, event, obs = y) ~ s(survtime, k = 6, bs = "ps"),
  gamma ~ s(x1, k = 6, bs = "ps"),
  mu ~ ti(id, bs = "re") + ti(obstime, bs = "ps", k = 8) + 
    ti(obstime, id, bs = "re") + 
    s(x2, k = 8, bs = "ps"),
  sigma ~ 1,
  alpha ~ 1,
  dalpha ~ -1
)
dat <- cbind(data.frame(survtime = b$y$`Surv2(survtime, event, obs = y)`[, 1],
                        event = b$y$`Surv2(survtime, event, obs = y)`[, 2],
                        y = b$y$`Surv2(survtime, event, obs = y)`[, 3]),
             b$model.frame)

b_rirs <- bamlss(formula = f, family = "jm", data = dat, 
                 maxit = 1500, timevar = "obstime", 
                 idvar = "id", subdivisions = 25, n.iter = 23000, 
                 burnin = 3000, thin = 20, nodf = TRUE, update.nu = TRUE)
newdat$mu_rirs <- predict(b_rirs, newdat, model = "mu")

p <- ggplot(data = newdat %>% filter(id %in% sampl), aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = mu_rirs), size = 1.2) +
  geom_point(aes(y = y, shape = id), size = 3) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  ggtitle("RI - RS")
pdf("../Fig1_1.pdf", width = 6, height = 6)
p
dev.off()




# Riga --------------------------------------------------------------------

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/B/sim100.Rdata"))
sampl <- c(16, 46, 95, 118, 126)
p <- ggplot(data = out$simdat$data %>% filter(id %in% ids), 
            aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = y), size = 1.2) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_grid(marker~., scales = "free_y")
pdf("../riga_Fig1.pdf", width = 6, height = 12)
p
dev.off()


ids <- c(13, 16, 27, 38, 43)

simdat$data_full$fit_est <- predict(out$b_est, out$simdat$data_full,
                                    model = "mu")
simdat$data_full$fit_tru <- predict(out$b_tru, out$simdat$data_full,
                                    model = "mu")
p <- ggplot(data = simdat$data_full %>% filter(id %in% ids), 
            aes(x = obstime, colour = id)) +
  theme_bw(base_size = 22) +
  geom_line(aes(y = fit_tru), size = 1.2) +
  geom_point(data = simdat$data %>% filter(id %in% ids), aes(y = y, shape = id), size = 3) +
  theme(legend.position = "none") +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  facet_wrap(marker~., nrow = 3, ncol = 2, scales = "free_y") +
  ggtitle("Data and Fit")
pdf("../riga_Fig4.pdf", width = 12, height = 12)
p
dev.off()

mfpca_es <- preproc_MFPCA(out$simdat$data, 
                          uni_mean = "y ~ 1 + obstime + x3 + obstime:x3",
                          nbasis = 10)
c(mfpca_es$values/sum(mfpca_es$values))[1:5]
dat <- data.frame(t = argvals(mfpca_es$functions)[[1]][[1]],
                  y = c(c(sapply(mfpca_es$functions, function(x) {
                    x@X[1, ]
                  })), c(sapply(mfpca_es$functions, function(x) {
                    x@X[2, ]
                  })), c(sapply(mfpca_es$functions, function(x) {
                    x@X[3, ]
                  }))),
                  MFPC = factor(rep(c("1 (51%)", "2 (31%)", "3 (15%)"), each = 6*101)),
                  marker = factor(rep(paste0("m", 1:6), each = 101)))
p <- ggplot(dat, aes(x = t, y = y, colour = MFPC)) +
  geom_line() +
  theme_bw(base_size = 22) +
  scale_y_continuous("Marker") +
  scale_x_continuous("Time") +
  #theme(legend.position = "none") +
  facet_wrap(marker~., nrow = 3, ncol = 2, scales = "free_y")
pdf("../riga_Fig3.pdf", width = 12, height = 12)
p
dev.off()
