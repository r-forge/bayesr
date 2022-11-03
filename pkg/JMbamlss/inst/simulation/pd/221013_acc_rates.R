location <- "workstation"


if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}


if (location == "server_linux") {
  results_wd <- paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
                       "simulation/")
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
}



library(bamlss)
library(tidyverse)
library(manipulate)


# Extract Acceptance Probabilities ----------------------------------------

extract_acc_prob_i <- function (i, setting, folder, name) {
  
  load(paste0(results_wd, setting, "/", folder, "/", i))
  
  acc <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)), 
                               fixed = TRUE)]
  su <- summary(acc)
  data.frame("q_2.5" = su$quantiles[, 1],"q_25" = su$quantiles[, 2],
             "q_50" = su$quantiles[, 3], "Mean" = su$statistics[, 1], 
             "q_75" = su$quantiles[, 4], "q_97.5" = su$quantiles[, 5],
             "model" = name)
  
}
extract_acc_prob <- Vectorize(extract_acc_prob_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
acc_info <- function (setting, folder, name = "Est_95") {
  
  m <- list.files(paste0(results_wd, setting, "/", folder,"/"))
  extract_acc_prob(m, setting = setting, folder = folder, name = name)
  
}

b_est_1 <- acc_info(setting = "scen_I_221019", folder = "bamlss_est_1",
                    "Est_1")
b_est_975 <- acc_info(setting = "scen_I_221019", folder = "bamlss_est",
                      "Est_975")
b_est_95 <- acc_info(setting = "scen_I_130922", folder = "bamlss_est",
                     "Est_95")
b_tru_95 <- acc_info(setting = "scen_I_130922", folder = "bamlss_tru",
                     "Tru_95")
b_tru <- acc_info(setting = "scen_I_051022", folder = "bamlss_tru",
                  "Tru")
save(b_est_1, b_est_975, b_est_95, b_tru_95, b_tru,
     file = paste0(results_wd, "acc_props.Rdata"))


# Acceptance Rates for JMbayes2 -------------------------------------------


extract_acc_rate_i <- function (i, setting, folder, name) {
  
  load(paste0(results_wd, setting, "/", folder, "/", i))
  
  acc <- jmb$acc_rates
  data.frame("rate" = c(mean(acc$bs_gammas), drop(acc$gammas), 
                        mean(acc$alphas), mean(acc$b), 
                        mean(c(acc$L, acc$sds)), mean(acc$sigmas)),
             "term" = c("lambda", "gamma", "alpha", "re", "D", "sigma"),
             "model" = name)
  
}
extract_acc_rate <- Vectorize(extract_acc_rate_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
acc_info_jmb <- function (setting, folder, name = "JMB") {
  
  m <- list.files(paste0(results_wd, setting, "/", folder,"/"))
  extract_acc_rate(m, setting = setting, folder = folder, name = name)
  
}

jmb_s <- acc_info_jmb(setting = "scen_I_130922", folder = "jmb")
save(jmb_s, file = paste0(results_wd, "acc_props_jmb.Rdata"))


# Combine Alpha Probabilities ---------------------------------------------


b_tru_its <- as.numeric(substr(names(b_tru), 2, 4))
n_tru1 <- b_tru_its[which(b_tru_its < 150)]
acc_res <- c(b_est_1, b_est_975, b_est_95, b_tru_95, b_tru)
acc_res <- do.call(rbind, Map(cbind, it = names(acc_res), acc_res)) %>%
  rownames_to_column(var = "Term") %>%
  mutate(it = as.numeric(substr(it, 2, 4)),
         model = ifelse(model == "Tru" & it %in% n_tru1, "Tru_1",
                        model),
         Model = factor(model, levels = c("Tru_1", "Tru",  "Tru_95", 
                                          "Est_1", "Est_975", "Est_95"),
                        labels = c("Tru_1", "Tru_975",  "Tru_95", 
                                   "Est_1", "Est_975", "Est_95")),
         term = gsub("[0-9]$", "", Term),
         term = substr(term, 12, nchar(term)-6),
         Predictor = factor(str_extract(term, "[^\\.]+"),
                            levels = c("lambda", "gamma", "alpha", "mu",
                                       "sigma")),
         term = gsub("model\\.matrix", "modelmatrix", term),
         term = gsub("s\\(id,fpc\\.([0-9]+)\\)", "\\1", term),
         Term = factor(sub(".*(\\.[sp]\\.)(.)", "\\2", term)),
         term = NULL) %>%
  pivot_longer(q_2.5:q_97.5, names_to = "metric", values_to = "Alpha") %>%
  mutate(metric = factor(metric, levels = c("q_2.5", "q_25", "q_50", "Mean",
                                            "q_75", "q_97.5"),
                         labels = c("2.5% Quantile", "25% Quantile",
                                    "50% Quantile", "Mean", "75% Quantile", 
                                    "97.5% Quantile")))


acc_res_jmb <- do.call(rbind, Map(cbind, it = names(jmb_s), jmb_s)) %>%
  mutate(Term = factor(term, levels = c("lambda", "gamma", "alpha", "re", 
                                        "D", "sigma")),
         Rate = rate)

# Plot Alpha Probabilities ------------------------------------------------

ggplot(acc_res %>% filter(Term %in% c("modelmatrix", "s(survtime)")), 
       aes(x = Predictor, y = Alpha, colour = Model)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(metric~.) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("Acceptance Probabilities for Non-PCRE Terms")

ggplot(acc_res %>% filter(!Term %in% c("modelmatrix", "s(survtime)")) %>%
         mutate(Term = factor(Term, levels = paste(1:12))), 
       aes(x = Term, y = Alpha, colour = Model)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(metric~.) +
  theme_bw() +
  ggtitle("Acceptance Probabilities for PCRE Terms")
  
ggplot(acc_res %>% filter(!Term %in% c("modelmatrix", "s(survtime)")) %>%
         mutate(Term = factor(Term, levels = paste(1:12))) %>%
         filter(Term %in% paste(1:4)), 
       aes(x = Term, y = Alpha, colour = Model)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(metric~.) +
  theme_bw() +
  ggtitle("Acceptance Probabilities for PCRE Terms")

ggplot(acc_res_jmb, aes(x = Term, y = Rate)) +
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  ggtitle("Acceptance Rate for JMB")


# Some Acceptance Probs are Close to 0 ------------------------------------

acc_res %>% 
  filter(metric == "97.5% Quantile", Alpha < 0.25)
# Iterations 183 and 191 seem problematic

acc_res %>%
  filter(model == "Est_975", metric == "97.5% Quantile", it == 183)
  



# Acceptance rates --------------------------------------------------------


load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/bamlss_tru",
            "/b100.Rdata"))
acc_ful <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                 fixed = TRUE)]
summary(acc_ful)
ful <- cbind(summary(acc_ful)$quantiles[, 2, drop = FALSE],
             summary(acc_ful)$statistics[, 1, drop = FALSE],
             summary(acc_ful)$quantiles[, 4, drop = FALSE])
xtable::xtable(ful)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/bamlss_tru",
            "/b100.Rdata"))
acc_95 <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                fixed = TRUE)]
summary(acc_95)
t95 <- cbind(summary(acc_95)$quantiles[, 2, drop = FALSE],
             summary(acc_95)$statistics[, 1, drop = FALSE],
             summary(acc_95)$quantiles[, 4, drop = FALSE])
print(xtable::xtable(t95), include.rownames = FALSE)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/bamlss_est",
            "/b100.Rdata"))
acc_est <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                 fixed = TRUE)]
summary(acc_est)
e95 <- cbind(summary(acc_est)$quantiles[, 2, drop = FALSE],
             summary(acc_est)$statistics[, 1, drop = FALSE],
             summary(acc_est)$quantiles[, 4, drop = FALSE])
print(xtable::xtable(e95), include.rownames = FALSE)



# Acceptance Rates for JMbayes2 -------------------------------------------


load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_221019/jmb",
            "/jmb_100.Rdata"))
sum(diff(jmb$mcmc$bs_gammas[[1]][, 1]) == 0)
jmb$acc_rates$bs_gammas

