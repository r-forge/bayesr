################################################################################
################################################################################
################################################################################
################################################################################

##### Problem of Alex: Sample nice functions that are correlated
##### Idea Paul: Use approach of Monod et al. (Gaussian process on integer lattice)

library(mgcv)
library(mvtnorm)

## auxiliary function ##########################################################

plot.design.matrix <- function(x,X,main="Design matrix",ylim=c(min(X),max(X)),col=FALSE)
{
  if(col==TRUE){
    col=1:dim(X)[2]
  } else {
    col=rep("black",dim(X)[2])
  }
  
  o=order(x)
  x=x[o]
  X=X[o,]
  
  plot(x,X[,1],ylim=ylim,type="l",main=main,col=col[1])
  for(j in 2:dim(X)[2]){
    lines(x,X[,j],col=col[j])
  }
}

## set up B-spline design for x ################################################

n=1000

x=1:n/n

a=0; b=1

numIntKnots=3 ## number of interior knots

fit=smoothCon(s(x,bs="bs",m=c(3,2),k=numIntKnots+4),data=as.data.frame(x))

B=fit[[1]]$X 
D=dim(B)[2] ## number of B-spline basis functions

plot.design.matrix(x,B)

## Gaussian process on integer lattice #########################################

nfunc=6 ## number of functions to be drawn

index=matrix(0,D*nfunc,2)
c=0
for(i in 1:D){
  for(j in 1:nfunc){
    c=c+1
    index[c,]=c(i,j)
  }
}
rm(c)
K=matrix(0,D*nfunc,D*nfunc)

for(i in 1:(D*nfunc)){
  for(j in 1:(D*nfunc)){
    K[i,j]=exp(-0.5*(index[i,1]-index[j,1])^2)*exp(-0.5*1/10*(index[i,2]-index[j,2])^2)
  }
}

## Note: The first length scale controls the wiggliness of the functions (the smaller, the more regular)
## The second length scale controls the similarity of the functions (the smaller, the more similar).

Beta=matrix(rmvnorm(n=1,sigma=K),D,nfunc,byrow=TRUE)

plot(x,B%*%Beta[,1],ylim=c(-2,2),type="l")
for(j in 2:nfunc){
  lines(x,B%*%Beta[,j],col=j)
}

cov((Beta))
k_ord <- c(seq(1, 36, by = 7), seq(2, 37, by = 7), seq(3, 38, by = 7),
           seq(4, 39, by = 7), seq(5, 40, by = 7), seq(6, 41, by = 7),
           seq(7, 42, by = 7))
k_ord
Sigma <- round(K[k_ord, k_ord], 3)
str(B)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

Sigma <- rbind(cbind(ar1_cor(n = 5, rho = 0.7), 0.3*ar1_cor(n = 5, rho = 0.3)),
               0.3*cbind(ar1_cor(n = 5, rho = 0.3), ar1_cor(n = 5, rho = 0.7)))
BetaAR1 <- matrix(rmvnorm(n=1,sigma=Sigma),D,nfunc,byrow=TRUE)
plot(x,B%*%Beta[,1],ylim=c(-2,2),type="l")
for(j in 2:nfunc){
  lines(x,B%*%Beta[,j],col=j)
}
