library(RandomFields)

#choose number of points + prediction points
n=1000
pred.n = 100

#generate coordinates for points
coord <- array(0,c(n,2))
for(i in 1:n){
  coord[i,] <- c(runif(1,0,1),runif(1,0,1))
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
  response[i] <- rnorm(1, 4 + 2*cov[i] + d$data[i])
  response2[i] <- rnorm(1, 41 + 21*cov[i] + c[i])
}


library(INLA)


mesh <- inla.mesh.2d(loc = coord, max.edge=c(35,100.0)) 

spde <- inla.spde2.matern(mesh=mesh, alpha=2) 


A <- inla.spde.make.A(mesh=mesh, loc=coord) 

s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)

#add data to stack
stack.norm <- inla.stack(data=list(Response=response), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                                 list(Covariate=cov)), tag="all")


stack.norm2 <- inla.stack(data=list(Response=response2), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                                 list(Covariate=cov)), tag="all")
#define formula for model
formula <- Response ~ -1 + Intercept + Covariate + f(spatial.field, model=spde)

output.pred <- inla(formula, data=inla.stack.data(stack.norm, spde=spde), family="gaussian", 
                    control.predictor=list(A=inla.stack.A(stack.norm), compute=TRUE),
                    control.compute=list(cpo=TRUE, dic=TRUE))

output.pred$summary.fixed

output.pred2 <- inla(formula, data=inla.stack.data(stack.norm2, spde=spde), family="gaussian", 
                    control.predictor=list(A=inla.stack.A(stack.norm2), compute=TRUE),
                    control.compute=list(cpo=TRUE, dic=TRUE))

output.pred2$summary.fixed