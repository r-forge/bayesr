

# Check whether the computational speedup changed the results -------------


sim_wd <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin",
                 ".de,share=volkmana.hub/JMbamlss/simulation/scen_I_130922/")
sim_wd1 <- paste0("/run/user/1000/gvfs/smb-share:server=clapton.wiwi.hu-berlin",
                  ".de,share=volkmana.hub/JMbamlss/simulation/scen_I_051022/")


# Tru_95 Model ------------------------------------------------------------

load(paste0(sim_wd, "data/d100.Rdata"))
load(paste0(sim_wd, "bamlss_tru/b100.Rdata"))
mfpca_tru <- attr(b_est, "FPCs")
nfpc <- 7 ## attr(b_est, "nfpc")

mfpca_tru_list <- lapply(seq_len(nfpc), function (i, mfpca = mfpca_tru) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
f_tru <- b_est$formula
d_rirs_tru <- attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)

b_new <- bamlss(f_tru, family = mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, sampler = FALSE)

all.equal(b_new$parameters, b_est$parameters)



# Est_95 Model ------------------------------------------------------------

load(paste0(sim_wd, "bamlss_est/b100.Rdata"))
mfpca_est <- attr(b_est, "FPCs")
nfpc <- 6 ## attr(b_est, "nfpc")
f_est <- b_est$formula

mfpca_est_list <- lapply(vals[seq_len(nfpc)], function (i, mfpca = mfpca_est) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
d_rirs_est <- attach_wfpc(mfpca_est, d_rirs$data, n = nfpc)

b_new <- bamlss(f_est, family = mjm_bamlss, data = d_rirs_est, 
                timevar = "obstime", maxit = 1500, sampler = FALSE)

all.equal(b_new$parameters, b_est$parameters)



# Tru_1 -------------------------------------------------------------------

load(paste0(sim_wd, "bamlss_tru/b100.Rdata"))
mfpca_tru <- attr(b_est, "FPCs")
load(paste0(sim_wd1, "bamlss_tru/b100.Rdata"))
nfpc <- 12 ## attr(b_est, "nfpc")

mfpca_tru_list <- lapply(seq_len(nfpc), function (i, mfpca = mfpca_tru) {
  list(functions = extractObs(mfpca$functions, i),
       values = mfpca$values[i])
})
f_tru <- b_est$formula
d_rirs_tru <- attach_wfpc(mfpca_tru, d_rirs$data, n = nfpc)

b_new <- bamlss(f_tru, family = mjm_bamlss, data = d_rirs_tru, 
                timevar = "obstime", maxit = 1500, sampler = FALSE)

all.equal(b_new$parameters, b_est$parameters)



# New Code for Est_95 Iteration 100 ---------------------------------------

load(paste0(sim_wd, "bamlss_est/b100.Rdata"))
b_old <- b_est

load(paste0(sim_wd, "b100.Rdata"))
b_new <- b_est
identical(b_old, b_new)

attr(b_old, "comp_time")
attr(b_new, "comp_time")
all.equal(b_new$parameters, b_old$parameters)
all.equal(samples(b_new), samples(b_old))
summary(samples(b_new, model = "alpha"))
summary(samples(b_old, model = "alpha"))
