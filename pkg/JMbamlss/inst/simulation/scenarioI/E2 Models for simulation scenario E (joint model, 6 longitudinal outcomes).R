
# Code from MAUFF et al. (2021)

# Specify location
location <- "server_linux"
if(location %in% c("server_linux", "server_windows")){
  .libPaths(if (location == "server_linux") {
    c("~/H:/volkmana.hub/R4_linux_b", "~/H:/volkmana.hub/R4_linux")
  } else "H:/R4_windows")
  setwd(if (location == "server_linux") "~/H:/volkmana.hub/JMbamlss"
        else "H:/JMbamlss")
  results_wd <- if(location == "server_linux") "simulation/"
} else {
  results_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-",
                       "berlin.de,share=volkmana.hub/JMbamlss/simulation/")
}
setting <- "scen_mauff/"

# ----------------------------------------------------------------------------

# Running models on data for simulation scenario E, (joint longitudinal and survival model, 6 longitudinal outcomes) ###################
# Note: in order to run this file, file E1 (data creation) needs to have been run in full. 
# Please note file locations.

# loading libraries

library("JMbayes")
library("MASS")
library("jagsUI")
library("plyr")

# running models and collecting results

for (datanum in c(1:200)) { 
  print(paste("Starting dataset ", datanum))
  
  startTime <- Sys.time()
  seed <- sample(c(7001:8000), 1, replace = FALSE)
  inputDataPath <- paste0(results_wd, setting, "/data/data", datanum,'.RData')
  
  load(inputDataPath)
  
  
  ####################################
  # Run mv model with JMbayes #
  ####################################
  
  # mv model for longitudinal outcomes  
  multMixed <- mvglmer(list(y1 ~  year + (year | id),
                            y2 ~  year + (year | id),
                            y3 ~  year + (year | id),
                            y4 ~  year + (year | id),
                            y5 ~  year + (year | id),
                            y6 ~  year + (year | id)),
                       data = dat, #n.processors = 2,
                       families = list(gaussian, gaussian, gaussian, gaussian, gaussian, gaussian))
  
  print("q1")
  
  cph <- coxph(Surv(Time, event) ~ group, dat.id, model = TRUE)
  
  # joint model based on mv longitudinal model and cox ph
  multJM_wRE <- mvJointModelBayes(multMixed, cph, timeVar = "year",# n_cores = 6, 
                                  update_RE = TRUE)
  print("q2")
  
  
  ####################
  # Collect results #
  ###################
  
  ind_wRE <- grep('Outcome', names(summary(multJM_wRE)), fixed = TRUE)
  
  res <- list(
    'multJM_wRE' = list(
      'survival' = list(
        'unweighted' = summary(multJM_wRE)$Survival,
        'weighted' = summary(multJM_wRE, TRUE)$Survival
      ),
      'longitudinal' = list(
        'unweighted' = do.call('rbind', summary(multJM_wRE)[ind_wRE]),
        'weighted' = do.call('rbind', summary(multJM_wRE, TRUE)[ind_wRE])
      ),
      'Dmat' = list(
        'unweighted' = summary(multJM_wRE)$D,
        'weighted' = summary(multJM_wRE, TRUE)$D
      )
    ),
    'seed' = seed,
    'times' = list(multMixed = multMixed$mcmc.info$elapsed.mins,
                   multJM_wRE = multJM_wRE$mcmc_info$elapsed_mins
    )
  )
  
  print(Sys.time() - startTime)
  
  resultPath <- paste0(results_wd, setting, "jmb/jmb_", datanum,'.rds')
  saveRDS(res,resultPath)
  # clear environment just to be safe
  
  # delete all unused objects
  which_to_delete <- function () {
    objs <- ls(pos = 1L)
    objs[!objs %in% c("results_wd", "setting")]
  }
  
  rm(list = which_to_delete())
}

