library(bamlss)

load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_130922/bamlss_tru",
            "/b110.Rdata"))


plot(b_est, which = "samples", model = "mu", ask = FALSE)
sam <- samples(b_est)

# Low autocorrelation
summary(sam[, which(colnames(sam) == "mu.s.s(id,fpc.1).b131")])
summary(sam[, which(colnames(sam) == "mu.s.s(id,fpc.2).b19")])
summary(sam[, which(colnames(sam) == "mu.s.s(id,fpc.2).b21")])

# High autocorrelation
summary(sam[, which(colnames(sam) == "mu.s.s(id,fpc.1).b132")])
summary(sam[, which(colnames(sam) == "mu.s.s(id,fpc.2).b16")])
summary(sam[, which(colnames(sam) == "mu.s.s(id,fpc.2).b20")])


# Comparison to JMbayes2
library(JMbayes2)
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,", 
            "share=volkmana.hub/JMbamlss/simulation/scen_I_221019/jmb/",
            "/jmb_110.Rdata"))
JMbayes2:::traceplot(jmb, parm = "betas", smooth = TRUE)
str(jmb$mcmc)
coda:::traceplot(jmb$mcmc$betas1, smooth = TRUE)
b_mcmc <- as.mcmc(matrix(jmb$mcmc$b[[1]], nrow = 1000, ncol = 150*12,
                         byrow = TRUE))
coda:::traceplot(b_mcmc, smooth = TRUE, ask = TRUE)
coda:::autocorr.plot(b_mcmc, lag.max = 20, ask = TRUE)
coda:::autocorr.plot(jmb$mcmc$betas1, lag.max = 20, ask = TRUE)

coda:::traceplot(sam[, which(colnames(sam) == 
                               "mu.p.model.matrix.markerm6:obstime:x3")],
                 smooth = TRUE)

# Longitudinal bamlss-model
load(paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin.de,",
            "share=volkmana.hub/JMbamlss/simulation/", 
            "scen_I_110_bamlss_long.Rdata"))
plot(b_est, which = "samples", model = "mu", ask = TRUE)
