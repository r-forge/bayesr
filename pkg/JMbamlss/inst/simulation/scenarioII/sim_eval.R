

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
local_wd

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


# Evaluate bamlss TRUE FPCs -----------------------------------------------

mnames_btru <- list.files(path = file.path(server_wd, "scen_II_230719", 
                                           "bamlss_tru"))
preds_btru <- JMbamlss:::sim_bamlss_predict(mnames_btru, server_wd, 
                                            "/scen_II_230719/bamlss_tru/",
                                            "/scen_II_230719/data/", rds = TRUE)
# saveRDS(preds_btru, 
#         file = file.path(server_wd, "scen_II_230719", "preds_btru.rds"))
saveRDS(preds_btru, file = "preds_btru.rds")

it_list <- JMbamlss:::sim_results(lapply(preds_btru, "[[", "predictions"),
                                  lapply(preds_btru, "[[", "simulations"),
                                  name = "TRU")
eval_btru <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                               it_list))
# saveRDS(eval_btru,
#         file = file.path(server_wd, "scen_II_230719", "eval_btru.rds"))
saveRDS(eval_btru,
        file = "eval_btru.rds")
rm(preds_btru, eval_btru, it_list)



# Evaluate bamlss EST 1 FPCs ----------------------------------------------

mnames_best1 <- list.files(path = file.path(server_wd, "scen_II_230719", 
                                            "bamlss_est1"))
preds_best1 <- JMbamlss:::sim_bamlss_predict(mnames_best1, server_wd, 
                                             "/scen_II_230719/bamlss_est1/",
                                             "/scen_II_230719/data/", 
                                             rds = TRUE)
# saveRDS(preds_best1, 
#         file = file.path(server_wd, "scen_II_230719", "preds_best1.rds"))
saveRDS(preds_best1, 
        file = "preds_best1.rds")

it_list <- JMbamlss:::sim_results(lapply(preds_best1, "[[", "predictions"),
                                  lapply(preds_best1, "[[", "simulations"),
                                  name = "EST1")
eval_best1 <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                                 it_list))
# saveRDS(eval_best1,
#         file = file.path(server_wd, "scen_II_230719", "eval_best1.rds"))
saveRDS(eval_best1,
        file = "eval_best1.rds")
rm(preds_best1, eval_best1, it_list)




# Evaluate bamlss EST 95 FPCs ---------------------------------------------

mnames_best95 <- list.files(path = file.path(server_wd, "scen_II_230719", 
                                            "bamlss_est95"))
preds_best95 <- JMbamlss:::sim_bamlss_predict(mnames_best95, server_wd, 
                                             "/scen_II_230719/bamlss_est95/",
                                             "/scen_II_230719/data/",
                                             rds = TRUE)
# saveRDS(preds_best95, 
#         file = file.path(server_wd, "scen_II_230719", "preds_best95.rds"))
saveRDS(preds_best95, 
        file = "preds_best95.rds")

it_list <- JMbamlss:::sim_results(lapply(preds_best95, "[[", "predictions"),
                                  lapply(preds_best95, "[[", "simulations"),
                                  name = "EST95")
eval_best95 <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                                 it_list))
# saveRDS(eval_best95,
#         file = file.path(server_wd, "scen_II_230719", "eval_best95.rds"))
saveRDS(eval_best95,
        file = "eval_best95.rds")
rm(preds_best95, eval_best95, it_list)



# Evaluate jmbayes  -------------------------------------------------------

mnames_jmb <- list.files(path = file.path(server_wd, "scen_II_230719", "jmb"))
preds_jmb <- JMbamlss:::sim_jmb_predict(mnames_jmb, server_wd, 
                                        "/scen_II_230719/jmb/",
                                        "/scen_II_230719/data/", rds = TRUE)
# saveRDS(preds_jmb, 
#         file = file.path(server_wd, "scen_II_230719", "preds_jmb.rds"))
saveRDS(preds_jmb, 
        file = "preds_jmb.rds")

it_list <- JMbamlss:::sim_results(lapply(preds_jmb, "[[", "predictions"),
                                  lapply(preds_jmb, "[[", "simulations"),
                                  name = "JMB")
eval_jmb <- do.call(rbind, Map(cbind, it = sub("\\.rds", "", names(it_list)),
                               it_list))
# saveRDS(eval_jmb,
#         file = file.path(server_wd, "scen_II_230719", "eval_jmb.rds"))
saveRDS(eval_jmb,
        file = "eval_jmb.rds")
rm(preds_jmb, eval_jmb, it_list)
