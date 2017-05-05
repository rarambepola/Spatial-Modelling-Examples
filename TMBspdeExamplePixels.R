setwd("J:/TMB")

ntotal = 10

error =0
errorspde = 0

library(RandomFields)
library(TMB)
library(INLA)

#TMB stuff

compile("TMBspdePixels.cpp")
dyn.load(dynlib("TMBspdePixels"))

for(k in 1:ntotal){

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

rb0 <- 0.25
rb1 <- 0.99
for(i in 1:n){
  cov[i] <- runif(1,0,3)
  #response[i] <- rnorm(1, rb0 + rb1*cov[i] + d$data[i], 1)
  #response[i] <- rnorm(1, 4 + 2*cov[i], 1)
  response[i] <- rnorm(1, rb0 + rb1*cov[i] + c[i], 1)
}



mesh <- inla.mesh.2d(loc = coord, max.edge=c(1,2), cutoff=0.25) 
plot(mesh)

spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]

idx <- mesh$idx$loc
n_s <- nrow(spde$M0)


f <- MakeADFun(
      data = list(X=response, cov=cov, spde=spde, idx=idx),
      parameters = list(beta5=0, beta1=0, sigma=1, log_kappa=2.5, x=runif(n_s,0,10)),
      random="x",
      DLL = "TMBspdePixels"
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

error = error + abs(as.numeric(fit2$par[1]) - rb0) + abs(as.numeric(fit2$par[2]) - rb1)
errorspde = errorspde + abs(as.numeric(fit$par[1]) - rb0) + abs(as.numeric(fit$par[2]) - rb1)
}

print(fit)
print(fit2)

print(errorspde)
print(error)

