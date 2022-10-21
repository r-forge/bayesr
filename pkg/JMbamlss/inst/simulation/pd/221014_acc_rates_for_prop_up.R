library(bamlss)
library(tidyverse)

# Check update every 10th step --------------------------------------------
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/scen_I_221013/b100.Rdata"))
b_ten <- b_est

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/scen_I_051022/bamlss_tru/",
            "b100.Rdata"))

acc_full <- samples(b_est)[, grep(".alpha", colnames(samples(b_est)),
                                  fixed = TRUE)]
summary(acc_full)

acc_ten <- samples(b_ten)[, grep(".alpha", colnames(samples(b_ten)), 
                                 fixed = TRUE)]
summary(acc_ten)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/scen_I_221013/", 
            "b100_2.Rdata"))
b_two <- b_est
acc_two <- samples(b_two)[, grep(".alpha", colnames(samples(b_two)), 
                                 fixed = TRUE)]
summary(acc_two)
