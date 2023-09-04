
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
                                           "simulation"),
                    "server_linux" = "~/H:/volkmana.hub/JMbamlss/simulation",
                    "server_windows" = "H:/JMbamlss/simulation")


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



# Extract FPCs from models ------------------------------------------------


extract_comptime_info_i <- function (i, results_wd, setting, folder) {
  
  b_est <- readRDS(file.path(results_wd, setting, folder, i))
  
  c(attr(b_est, "comp_time"))
}
extract_comptime_info <- Vectorize(extract_comptime_info_i,
                                   vectorize.args = "i")
comptime_info <- function (results_wd, setting, folder) {
  
  m <- list.files(file.path(results_wd, setting, folder))
  extract_comptime_info(m, results_wd = results_wd, setting = setting,
                        folder = folder)
  
}



# Extract Computation Times -----------------------------------------------

# Simulation Scenario I
simI_tru <- comptime_info(results_wd = server_wd,
                          setting = "scen_I_230719",
                          folder = "bamlss_tru")
simI_est1 <- comptime_info(results_wd = server_wd,
                           setting = "scen_I_230719",
                           folder = "bamlss_est1")
simI_est95 <- comptime_info(results_wd = server_wd,
                            setting = "scen_I_230719",
                            folder = "bamlss_est95")
simI_jmb <- comptime_info(results_wd = server_wd,
                          setting = "scen_I_230719",
                          folder = "jmb")

load(file.path(server_wd, "scen_I_230719", "comp_times.Rdata"))
# bamlss_true: mavis - 100:128, springsteen - 129:199, afflux - 200:249,
# pandia 250:299
summary(simI_tru["elapsed", 1:29])
summary(simI_tru["elapsed", 30:100])
# bamlss_est1: springsteen - 100:299
summary(simI_est1["elapsed", ])
# bamlss_est95: mavis - 100:299
summary(simI_est95["elapsed", ])
# jmb: springsteen - 100:299
summary(simI_jmb["elapsed", ])


# Simulation Scenario II
simII_tru <- comptime_info(results_wd = server_wd,
                          setting = "scen_II_230719",
                          folder = "bamlss_tru")
simII_est1 <- comptime_info(results_wd = server_wd,
                           setting = "scen_II_230719",
                           folder = "bamlss_est1")
simII_est95 <- comptime_info(results_wd = server_wd,
                            setting = "scen_II_230719",
                            folder = "bamlss_est95")
simII_jmb <- comptime_info(results_wd = server_wd,
                          setting = "scen_II_230719",
                          folder = "jmb")