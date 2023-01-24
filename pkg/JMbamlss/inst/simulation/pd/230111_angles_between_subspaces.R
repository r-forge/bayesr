# Objects for all clusters ------------------------------------------------

# Number of individuals and other quantities
n <- 300
argvals <- seq(0, 1, by = 0.01)
x <- seq(0, 1, by = 0.1)

# Random covariance matrix
# Set the eigenvalues but the eigenvectors are random
set.seed(1105)
p <- 6
P <- qr.Q(qr(matrix(rnorm(p^2), ncol = p)))
evals <- c(4, 3, 2, 1, 0.5, 0.2)
cov <- crossprod(P, P*(evals))


# Find spline functions
# Marker1
m1sp1 <- splinefun(x, c(0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))
m1sp2 <- splinefun(x, c(0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4, 0.5, 0.5, 0.5))
m1sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.6, 0.5, 0.4))
# Marker2
m2sp1 <- splinefun(x, c(0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05, 0, 0, 0))
m2sp2 <- splinefun(x, c(0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3, 0.2, 0.1, 0.05))
m2sp3 <- splinefun(x, c(0, 0, 0, 0, 0, 0, 0.5, 1, 0.7, 0.5, 0.3))

m1 <- funData(argvals = argvals,
              X = matrix(c(m1sp1(argvals), m1sp2(argvals), m1sp3(argvals)), 
                         nrow = 3, byrow = TRUE))
m2 <- funData(argvals = argvals,
              X = matrix(c(m2sp1(argvals), m2sp2(argvals), m2sp3(argvals)), 
                         nrow = 3, byrow = TRUE))

# True multivariate covariance structure
m <- JMbamlss:::MFPCA_cov(cov = cov, basis_funs = list(m1, m2))



# Angles ------------------------------------------------------------------


prangles <- function(A, B) {
  if (is.vector(A)) A = matrix(A, ncol=1)
  if (is.vector(B)) B = matrix(B, ncol=1)
  
  Qa = qr.Q(qr(A))
  Qb = qr.Q(qr(B))
  
  C = svd(crossprod(Qa, Qb))$d
  C = vapply(C, function(x) min(x, 1), numeric(1))
  angles = sort(acos(C), decreasing=TRUE)
  return(angles)
}


A <- t(mfpca_est$functions@.Data[[1]]@X)
B <- t(m$functions@.Data[[1]]@X)
Qa <- qr.Q(qr(A))
Qb <- qr.Q(qr(B))
C <- svd(crossprod(Qa, Qb))$d
C <- vapply(C, function(x) min(x, 1), numeric(1))
angles <- sort(acos(C), decreasing=TRUE)
acos(Reduce("*", C))
acos(Reduce("*", cos(angles)))

P <- qr.Q(qr(matrix(rnorm(6^2), ncol = p)))
B1 <- t(P %*% t(B))
prangles(B, B1)
acos(Reduce("*", cos(prangles(B, B1))))

# PRACMA!
# pracma::subspace
pracma::subspace(B, B1)



# Rotation ----------------------------------------------------------------

X <- t(cbind(m$functions@.Data[[1]]@X, m$functions@.Data[[2]]@X))
R <- fda::varmx(X)
X_r <- X %*% R

aha <- multiFunData(funData(argvals(m$functions@.Data[[1]]),
                            t(X_r[1:101, ])),
                    funData(argvals(m$functions@.Data[[2]]),
                            t(X_r[102:202, ])))
