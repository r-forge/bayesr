
# Prediction Function Check UP --------------------------------------------

# JMbamlss predict function does not care what FPC information you supply with 
# the newdata, it will always use the available FPC information from the formula
# to create the original FPC basis






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


load(paste0(server_wd, "scen_II_221209/A/b100.Rdata"))
load(paste0(server_wd, "scen_II_221209/data/d100.Rdata"))

newdata_full <- d_rirs$data_full
newdata_fpc0 <- d_rirs$data_full %>%
  mutate(fpc.1 = 0, fpc.2 = 0, fpc.3 = 0, fpc.4 = 0, fpc.5 = 0, fpc.6 = 0)
newdata_NULL <- d_rirs$data_full %>%
  mutate(fpc.1 = NULL, fpc.2 = NULL, fpc.3 = NULL, fpc.4 = NULL, fpc.5 = NULL,
         fpc.6 = NULL)

p_full <- predict(b_est, newdata = newdata_full, model = "mu")
p_fpc0 <- predict(b_est, newdata = newdata_fpc0, model = "mu")
p_NULL <- predict(b_est, newdata = newdata_NULL, model = "mu")
all.equal(p_full, p_fpc0)
all.equal(p_full, p_NULL)

p_fri <- predict(b_est, newdata = newdata_full, model = "mu", term = "id")

newdata_full <- newdata_full %>%
  mutate(p_fri = p_fri, p_mu = p_full, p_wofri = p_full - p_fri)


ggplot(newdata_full, aes(x = obstime, y = p_fri, colour = id)) +
  facet_wrap(~marker) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

ggplot(newdata_full, aes(x = obstime, y = p_mu, colour = id)) +
  facet_wrap(~marker) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")

ggplot(newdata_full, aes(x = obstime, y = p_wofri, colour = id)) +
  facet_wrap(~marker) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none")
