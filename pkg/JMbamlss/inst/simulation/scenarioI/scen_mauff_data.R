# To Create The Data for the MAUFF Simulation
# See File "E1 Data creation for simulation scenario E"


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
  "workstation" = paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.",
                         "hu-berlin.de,share=volkmana.hub/JMbamlss/simulation/"),
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

# Setting for the simulation
start <- 1
stop <- 200
setting <- "/scen_mauff"
Sys.time()
sessionInfo()

for (i in start:stop) {
  
  # Load the data
  load(paste0(server_wd, setting, "/data/data", i, ".Rdata"))
  
  max_longi_time <- max(dat$year)
  max_surv_time <- max(dat$Time)
  
  d_sim <- pivot_longer(dat, y1:y6, names_to = "marker", values_to = "y") %>%
    mutate(marker = factor(marker, labels = paste0("m", 1:6)),
           id = factor(id),
           Time1 = Time/max_surv_time,
           year1 = year/max_surv_time,
           Time1cens = ifelse(Time/max_longi_time > 1, 1, Time/max_longi_time),
           year1cens = year/max_longi_time,
           event1cens = ifelse(Time/max_longi_time > 1, 0, event)) %>%
    arrange(marker, id, year) %>%
    as.data.frame()
  
  saveRDS(d_sim, file = paste0(server_wd, setting, "/data/data", i, ".rds"))
}
