library(RandomFields)
library(TMB)

#choose number of points + prediction points
n=1000
pred.n = 100

#generate coordinates for points
coord <- array(0,c(n,2))
for(i in 1:n){
  coord[i,] <- c(runif(1,0,10),runif(1,0,10))
}

a = RMmatern(0.5)
b = RFsimulate(a, x=coord[,1], y=coord[,2])

c = as.matrix(b)

plot(coord, cex = c, main = "Matern RF (Random Fields)")

library(geoR)

d = grf(n, grid = coord, xlims = c(0, 10), ylims = c(0, 10), nsim = 1, cov.model = "matern",
        cov.pars = c(1,1), 
        kappa = 0.5, nugget = 0, lambda = 1,
        mean = 0, RF=TRUE)

plot(coord, cex=d$data, main = "Matern RF (geoR)")

cov <- rep(0, n)
response <- rep(0, n)
response2 <- rep(0, n)

for(i in 1:n){
  cov[i] <- runif(1,0,3)
  response[i] <- rnorm(1, 4 + 2*cov[i] + d$data[i], 1)
  #response[i] <- rnorm(1, 4 + 2*cov[i], 1)
  response2[i] <- rnorm(1, 41 + 21*cov[i] + c[i], 1)
}

#TMB stuff
compile("TMBExample.cpp")
dyn.load(dynlib("TMBExample"))

f <- MakeADFun(
      data = list(x=response, cov=cov),
      parameters = list(beta0=0, beta1=0, sigma=1),
      DLL = "TMBExample"
)

fit = nlminb(f$par,f$fn,f$gr,lower=c(-10,-10,0),upper=c(10.0,10.0,10))
print(fit)


