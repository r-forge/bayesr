library(bamlss)
library(tidyverse)
library(manipulate)


# Acceptance rates --------------------------------------------------------


load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/bamlss_tru",
            "/b100.Rdata"))
acc_ful <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                 fixed = TRUE)]
summary(acc_ful)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/bamlss_tru",
            "/b100.Rdata"))
acc_95 <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                fixed = TRUE)]
summary(acc_95)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/bamlss_est",
            "/b100.Rdata"))
acc_est <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                 fixed = TRUE)]
summary(acc_est)


load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/jmb",
            "/jmb_100.Rdata"))
