library(bamlss)
data("pbc2", package = "JMbayes")

## Set up the model formula including
## functional random intercepts using ti().
f_bamlss <- list(
  Surv2(years, status2, obs = log(serBilir)) ~ s(years,k=20),
  gamma ~ s(age,k=20) + drug + sex,
  mu ~ ti(id,bs="re") + 
    ti(year,k=20) + 
    ti(id,year,bs=c("re","cr"),k=c(nlevels(pbc2$id),8)),
  sigma ~ 1,
  alpha ~ s(years,k=20),
  dalpha ~ -1
)
## Set the seed for reproducibility.
set.seed(123)

debug(bamlss:::sam_JM)
## Estimate model.
pbc2 <- pbc2[1:300, ]
b <- bamlss(f_bamlss, data = pbc2, family = "jm",
            timevar = "year", idvar = "id", optimizer = FALSE,
            n.iter = 4, burnin = 1, step = 1)


# --- Code from JM.R
# 312 Individuen und 25 Integrationspunkte = 7800 Zeilen
# für s(years) term, k = 20 -> also 19 Parameter = 19 Spalten
X <- x$fit.fun_timegrid(NULL)

# eta_timegrid: für jeden Integrationspunkt sind die zeitvariierenden etas 
#   ausgewertet und aufaddiert
eeta <- exp(eta_timegrid)

# Integral wird berechnet? Funktion in C
# 
int <- survint(X, eeta, width, exp(eta$gamma))
# returns negative Hessian

xgrad <- drop(t(status) %*% x$XT - int$grad)