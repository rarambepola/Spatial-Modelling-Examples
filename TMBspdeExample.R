setwd("J:/INLAExamples")

beta0hist = rep(0, (ntotal <- 50))
beta1hist = rep(0, ntotal)

library(RandomFields)
library(TMB)
library(INLA)
  

#TMB stuff

compile("TMBspdeExample.cpp")
dyn.load(dynlib("TMBspdeExample"))


#choose number of points + prediction points
n=200
pred.n = 100


#generate coordinates for points
coord <- array(0,c(n,2))
for(i in 1:n){
  coord[i,] <- c(runif(1,0,10),runif(1,0,10))
}

a <- RMmatern(0.5)
b <- RFsimulate(a, x=coord[,1], y=coord[,2])

c <- as.matrix(b)

plot(coord, cex = c, main = "Matern RF (Random Fields)")

library(geoR)

d <- grf(n, grid = coord, xlims = c(0, 10), ylims = c(0, 10), nsim = 1, cov.model = "matern",
        cov.pars = c(1,1), 
        kappa = 0.5, nugget = 0, lambda = 1.00,
        mean = 0, RF=TRUE)

plot(coord, cex=d$data, main = "Matern RF (geoR)")

cov <- rep(0, n)
response <- rep(0, n)
response2 <- rep(0, n)

for(i in 1:n){
  cov[i] <- runif(1,0,3)
  response[i] <- rnorm(1, 0.5 + 0.5*cov[i] + d$data[i], 1)
  #response[i] <- rnorm(1, 4 + 2*cov[i], 1)
  #response[i] <- rnorm(1, 1 + 2*cov[i] + c[i], 1)
}



mesh <- inla.mesh.2d(loc = coord, max.edge=c(1,2), cutoff=0.25) 
plot(mesh)

spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]


n_s <- nrow(spde$M0)
#n_s <- n


f <- MakeADFun(
      data = list(X=response, cov=cov, spde=spde),
      parameters = list(beta0=0, beta1=0, sigma=1, log_kappa=0.5, x=rep(10.0, n_s)),
      random="x",
      DLL = "TMBspdeExample"
)
fit <- nlminb(f$par,f$fn,f$gr,lower=c(-10,-10,0,0))

compile("TMBExample.cpp")
dyn.load(dynlib("TMBExample"))

g <- MakeADFun(
  data = list(x=response, cov=cov),
  parameters = list(beta0=0, beta1=0, sigma=1),
  DLL = "TMBExample"
)

fit2 <- nlminb(g$par,g$fn,g$gr,lower=c(-10,-10,0))

print(fit$par)
print(fit2$par)
