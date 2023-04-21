get_to_wd <- "../JMbamlss/"

source(paste0(get_to_wd, "inst/simulation/scenarioI/create_problematic_data.R"))
library(JMbayes2)
library(bamlss)
library(funData)
server_wd <- paste0("/run/user/1000/gvfs/smb-share:",
                    "server=clapton.wiwi.hu-berlin.de,",
                    "share=volkmana.hub/JMbamlss/simulation/")



# Set up all objects ------------------------------------------------------


# Basis functions on each dimension
seq1 <- seq(0, max(simdat$Time), length.out = 101)
b_funs <- rep(list(funData(argvals = seq1,
                           X = matrix(c(rep(1, length(seq1)), seq1),
                                      byrow = TRUE, ncol = length(seq1)))), 6)

mfpca_tru <- JMbamlss:::MFPCA_cov(cov = as.matrix(D), basis_funs = b_funs,
                                  scores = b)
mfpca_list <- lapply(seq_along(mfpca_tru$values), 
                     function (i, mfpca = mfpca_tru) {
                       list(functions = extractObs(mfpca$functions, i),
                            values = mfpca$values[i])
                     })

simdat_tru <- JMbamlss:::attach_wfpc(mfpca_tru, simdat,
                                     n = length(mfpca_tru$values),
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

# For faster problems, do not optimize nu
# set.seed(1)
# b_mul <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = simdat_tru, 
#                 timevar = "year", maxit = 1500, verbose  = TRUE,
#                 par_trace = TRUE)
# Crashes relaitvely late

# set.seed(1)
# b_mul <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = simdat_tru, 
#                 timevar = "year", maxit = 493, verbose  = TRUE, 
#                 sampler = FALSE, par_trace = TRUE)
# saveRDS(b_mul, file = paste0(server_wd, "scen_mauff/mul/bamlss_mul_",
                             # "nonuup.Rds"))
b_mul <- readRDS(paste0(server_wd, "scen_mauff/mul/bamlss_mul_",
                        "nonuup.Rds"))


# Analyse the traceplots --------------------------------------------------

# Traceplots
lambda <- sapply(b_mul$model.stats$optimizer$par_trace, 
                 function(x) x$lambda$s[[1]])
gamma <- sapply(b_mul$model.stats$optimizer$par_trace, 
                function(x) x$gamma$p)
alpha <- sapply(b_mul$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
mus1 <- sapply(b_mul$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[1]])
mus2 <- sapply(b_mul$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[2]])
mup <- sapply(b_mul$model.stats$optimizer$par_trace, 
              function(x) x$mu$p)
sigma <- sapply(b_mul$model.stats$optimizer$par_trace, 
                function(x) x$sigma$p)

ggplot(data = data.frame(t(lambda)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Lambda Parameters")

ggplot(data = data.frame(t(gamma)) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = -5.8, linetype = "dotted", color = "#00BFC4") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Gamma Parameters")

ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters") +
  ylim(c(-10,10))

ggplot(data = data.frame(t(mus1)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Random Intercept Parameters")

ggplot(data = data.frame(t(mus2)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Random Slope Parameters")

ggplot(data = data.frame(t(mup)) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = -0.3, linetype = "dotted", color = "#00BFC4") +
  geom_hline(yintercept = 1.71, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Parameters")

ggplot(data = data.frame(t(sigma)) %>%
         mutate(it = 1:493) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = log(0.57), linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Sigma Parameters")


# Adapt variances ---------------------------------------------------------

debugonce(JMbamlss:::MJM_transform)
set.seed(1)
b_multau <- bamlss(f_tru, family = JMbamlss:::mjm_bamlss, data = simdat_tru,
                timevar = "year", maxit = 448, verbose  = TRUE,
                par_trace = TRUE, sampler = FALSE,
                tau = list("mu" = list("s(id,fpc.1)" = mfpca_tru$values[1],
                                       "s(id,fpc.2)" = mfpca_tru$values[2],
                                       "s(id,fpc.3)" = mfpca_tru$values[3],
                                       "s(id,fpc.4)" = mfpca_tru$values[4],
                                       "s(id,fpc.5)" = mfpca_tru$values[5],
                                       "s(id,fpc.6)" = mfpca_tru$values[6],
                                       "s(id,fpc.7)" = mfpca_tru$values[7],
                                       "s(id,fpc.8)" = mfpca_tru$values[8],
                                       "s(id,fpc.9)" = mfpca_tru$values[9],
                                       "s(id,fpc.10)" = mfpca_tru$values[10],
                                       "s(id,fpc.11)" = mfpca_tru$values[11],
                                       "s(id,fpc.12)" = mfpca_tru$values[12])))
saveRDS(b_multau, file = paste0(server_wd, "scen_mauff/mul/bamlss_multau_",
                                "nonuup.Rds"))

# Traceplots
lambda <- sapply(b_multau$model.stats$optimizer$par_trace, 
                 function(x) x$lambda$s[[1]])
gamma <- sapply(b_multau$model.stats$optimizer$par_trace, 
                function(x) x$gamma$p)
alpha <- sapply(b_multau$model.stats$optimizer$par_trace, 
                function(x) x$alpha$p)
mus1 <- sapply(b_multau$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[1]])
mus2 <- sapply(b_multau$model.stats$optimizer$par_trace, 
               function(x) x$mu$s[[2]])
mup <- sapply(b_multau$model.stats$optimizer$par_trace, 
              function(x) x$mu$p)
sigma <- sapply(b_multau$model.stats$optimizer$par_trace, 
                function(x) x$sigma$p)

ggplot(data = data.frame(t(lambda)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Lambda Parameters")

ggplot(data = data.frame(t(gamma)) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = -5.8, linetype = "dotted", color = "#00BFC4") +
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Gamma Parameters")

ggplot(data = data.frame(t(alpha)) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = 0.75, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Alpha Parameters") +
  ylim(c(-10,10))

ggplot(data = data.frame(t(mus1)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Random Intercept Parameters")

ggplot(data = data.frame(t(mus2)) %>%
         select(-tau21, -edf) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Random Slope Parameters")

ggplot(data = data.frame(t(mup)) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = -0.3, linetype = "dotted", color = "#00BFC4") +
  geom_hline(yintercept = 1.71, linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Mu Parameters")

ggplot(data = data.frame(t(sigma)) %>%
         mutate(it = 1:448) %>%
         pivot_longer(cols = -it),
       aes(x = it, y = value, col = name)) +
  geom_line() +
  geom_hline(yintercept = log(0.57), linetype = "dotted", color = "#F8766D") +
  theme_bw() +
  theme(legend.position = "None") +
  labs(x = "Iteration", y =  "Estimate") +
  ggtitle("Sigma Parameters")

