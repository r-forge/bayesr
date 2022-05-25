library(manipulate)
load("inst/objects/find_sim030522.Rdata")


# Lambda
lambda_pars <- data.frame(
  t(sapply(b_sim1$model.stats$optimizer$par_trace, function (i) {
  i$lambda$s[[1]]
})))
l_pars <- lambda_pars %>%
  rownames_to_column(var = "It") %>% 
  mutate(It = as.numeric(It)) %>%
  pivot_longer(-It)
ggplot(l_pars, aes(x = It, y = value)) + 
  geom_line() +
  facet_wrap(~name, scales = "free") +
  geom_vline(xintercept = 967, color = "red") +
  geom_vline(xintercept = 227, color = "blue", linetype = "dashed") +
  ggtitle("Lambda")

# Gamma
gamma_pars <- data.frame(
  t(sapply(b_sim1$model.stats$optimizer$par_trace, function (i) {
    i$gamma$p
  })))
g_pars <- lambda_pars %>%
  rownames_to_column(var = "It") %>% 
  mutate(It = as.numeric(It)) %>%
  pivot_longer(-It)
ggplot(g_pars, aes(x = It, y = value)) + 
  geom_line() +
  facet_wrap(~name, scales = "free") +
  geom_vline(xintercept = 967, color = "red") +
  geom_vline(xintercept = 227, color = "blue", linetype = "dashed") +
  ggtitle("Gamma")

# Alpha
alpha_pars <- data.frame(
  t(sapply(b_sim1$model.stats$optimizer$par_trace, function (i) {
    i$alpha$p
  })))
a_pars <- alpha_pars %>%
  rownames_to_column(var = "It") %>% 
  mutate(It = as.numeric(It)) %>%
  pivot_longer(-It)
ggplot(a_pars, aes(x = It, y = value)) + 
  geom_line() +
  facet_wrap(~name, scales = "free") +
  geom_vline(xintercept = 967, color = "red") +
  geom_vline(xintercept = 227, color = "blue", linetype = "dashed") +
  ggtitle("Alpha")

# Mu
mu_p_pars <- data.frame(
  t(sapply(b_sim1$model.stats$optimizer$par_trace, function (i) {
    i$mu$p
  })))
m_p_pars <- mu_p_pars %>%
  rownames_to_column(var = "It") %>% 
  mutate(It = as.numeric(It)) %>%
  pivot_longer(-It)
ggplot(m_p_pars, aes(x = It, y = value)) + 
  geom_line() +
  facet_wrap(~name, scales = "free") +
  geom_vline(xintercept = 967, color = "red") +
  geom_vline(xintercept = 227, color = "blue", linetype = "dashed") +
  ggtitle("Mu")

mu_s_pars <- data.frame(
  t(sapply(b_sim1$model.stats$optimizer$par_trace, function (i) {
    i$mu$s[[1]]
  })))
m_s_pars <- mu_s_pars %>%
  rownames_to_column(var = "It") %>% 
  mutate(It = as.numeric(It)) %>%
  pivot_longer(-It)
manipulate(plot(1:5, cex=size), size = slider(1,150))
manipulate(
  {ggplot(m_s_pars %>% filter(name %in% paste0("b", (x-1)*8+1:8)), 
          aes(x = It, y = value)) + 
    geom_line() +
    facet_wrap(~name, scales = "free", nrow = 2, ncol = 4) +
    geom_vline(xintercept = 967, color = "red") +
    geom_vline(xintercept = 227, color = "blue", linetype = "dashed") +
    ggtitle("Mu")},
  x = slider(1, 150)
)

# Sigma
sigma_pars <- data.frame(
  t(sapply(b_sim1$model.stats$optimizer$par_trace, function (i) {
    i$sigma$p
  })))
s_pars <- sigma_pars %>%
  rownames_to_column(var = "It") %>% 
  mutate(It = as.numeric(It)) %>%
  pivot_longer(-It)
ggplot(s_pars, aes(x = It, y = value)) + 
  geom_line() +
  facet_wrap(~name, scales = "free") +
  geom_vline(xintercept = 967, color = "red") +
  geom_vline(xintercept = 227, color = "blue", linetype = "dashed") +
  ggtitle("Sigma")

debug(opt_MJM)
b_sim750 <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
                   timevar = "obstime", maxit = 2, sampler = FALSE,
                   start = b_sim1$model.stats$optimizer$par_trace[[750]],
                   par_trace = TRUE)
b_sim999 <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
                   timevar = "obstime", maxit = 2, sampler = FALSE,
                   start = b_sim1$model.stats$optimizer$par_trace[[999]],
                   par_trace = TRUE)
b_sim2 <- bamlss(f, family = mjm_bamlss, data = d_indepri$data, 
                   timevar = "obstime", maxit = 2, sampler = FALSE,
                   start = b_sim1$model.stats$optimizer$par_trace[[2]],
                   par_trace = TRUE)
