library(JMbayes2)
# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum cholesterol with random ordering
set.seed(123)
pbc2_ro <- pbc2[sample(seq_len(nrow(pbc2)), nrow(pbc2)), ]
fm1 <- lme(log(serChol) ~ year * sex, data = pbc2_ro, random = ~ year | id,
           na.action = na.omit)
fm2 <- lme(log(alkaline) ~ year * sex, data = pbc2_ro, random = ~ year | id,
           na.action = na.omit)

# the joint model throws an error
jointFit1 <- jm(CoxFit, list(fm1), time_var = "year",
               n_iter = 1200L, n_burnin = 200L, n_thin = 2L, 
               control = list(cores = 1))
jointFit3 <- jm(CoxFit, list(fm1, fm2), time_var = "year",
                n_iter = 1200L, n_burnin = 200L, n_thin = 2L, 
                control = list(cores = 3))
jointFit1_1 <- jm(CoxFit, list(fm1), time_var = "year",
                n_iter = 1200L, n_burnin = 200L, n_thin = 2L, 
                control = list(cores = 1, n_chains = 1))
