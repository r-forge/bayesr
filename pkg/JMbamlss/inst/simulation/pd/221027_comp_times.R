
# Compare computing times -------------------------------------------------

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


extract_comp_i <- function (i, setting, folder, name) {
  
  load(paste0(results_wd, setting, "/", folder, "/", i))
  
  attr(b_est, "comp_time")
  
}
extract_comp <- Vectorize(extract_comp_i, vectorize.args = "i",
                              SIMPLIFY = FALSE)
comp_info <- function (setting, folder, name = "Est_95") {
  
  m <- list.files(paste0(results_wd, setting, "/", folder,"/"))
  extract_comp(m, setting = setting, folder = folder, name = name)
  
}

b_est_1 <- comp_info(setting = "scen_I_221019", folder = "bamlss_est_1", 
                     name = "Est_1")
b_est_975 <- comp_info(setting = "scen_I_221019", folder = "bamlss_est", 
                       name = "Est_975")
b_est_975 <- comp_info(setting = "scen_I_221019", folder = "bamlss_est", 
                       name = "Est_975")

comp_dat <- data.frame(t = sapply(b_est_975, "[[", 3)) %>%
  rownames_to_column() %>%
  mutate(it = as.numeric(substr(rowname, 2, 4)),
         Code = factor(ifelse(it < 150, "old", "new"), 
                       levels = c("old", "new")),
         Min = t /60)

ggplot(comp_dat, aes(y = Min, x = Code)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Computation Time for Est_975 Before and After SpeedUp")


comp_dat <- data.frame(t = sapply(b_est_1, "[[", 3)) %>%
  rownames_to_column() %>%
  mutate(it = as.numeric(substr(rowname, 2, 4)),
         Code = factor(ifelse(it < 107, "old", "new"), 
                       levels = c("old", "new")),
         Min = t /60)
ggplot(comp_dat, aes(y = Min, x = Code)) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Computation Time for Est_1 Before and After SpeedUp")
