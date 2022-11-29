
# Convergence of Sampler --------------------------------------------------



location <- "workstation"


if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
}

library(tidyverse)


if (location == "server_linux") {
  results_wd <- paste0("/home/RDC/volkmana.hub/H:/volkmana.hub/JMbamlss/", 
                       "simulation/")
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
}

extract_conv <- Vectorize(function (i, setting, folder, name) {
  
  load(paste0(results_wd, setting, "/", folder, "/", i))
  
  b_est$model.stats$optimizer$converged
  
}, vectorize.args = "i", SIMPLIFY = FALSE)

conv_info <- function (setting, folder, name = "Est_95") {
  
  m <- list.files(paste0(results_wd, setting, "/", folder,"/"))
  extract_conv(m, setting = setting, folder = folder, name = name)
  
}

b_est_1 <- conv_info(setting = "scen_I_221019", folder = "bamlss_est_1", 
                     name = "Est_1")
b_est_975 <- conv_info(setting = "scen_I_221019", folder = "bamlss_est", 
                       name = "Est_975")
b_est_95 <- conv_info(setting = "scen_I_130922", folder = "bamlss_est", 
                       name = "Est_95")
