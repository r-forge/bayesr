library("bamlss")

data("rent", package = "gamlss.data")

rent$R2 <- rent$R / 1000

f <- R2 ~ Sp + Sm + B + H + L + loc + s(Fl) + s(A) 

b0 <- bamlss(f, family = gamma, data = rent, engine = "JAGS")

f <- R2 ~ sx(Fl) + sx(A)

b1 <- bamlss(f, family = gamma, data = rent, engine = "BayesX")
