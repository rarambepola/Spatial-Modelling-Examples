setwd("J:/INLAExamples")

library(RandomFields)
library(TMB)
library(INLA)

###coefficients
b0 <- 1.0
b1 <- 0.7
b2 <- 1.0

##Functions
#function to find which box number a point is in
#(boxes are numbered left to right, bottom to top, e.g. lowest row
#is 1,2,3,4 next is 5,6,7,8 etc. 
#but how they're numbered doesn't really matter)
boxfind <- function(coord, box.x.length, box.y.length, box.y.n){
  
  x.num <- max(1, ceiling(coord[1]/box.x.length))
  y.num <- max(1, ceiling(coord[2]/box.y.length))
  
  return(y.num + (x.num-1)*box.y.n)
}

#function to plot and colour pixels based on a vector of values
pixelplot <- function(values, total.x.size, total.y.size, x.n, y.n, main=NULL){

  x.length <- total.x.size/x.n
  y.length <- total.y.size/y.n
  
  rbPal <- colorRampPalette(c('red','blue'))
  col <- rbPal(10)[as.numeric(cut(values,breaks = 10))]
  
  plot(x=NULL, y=NULL, xlim=range(0:total.x.size), ylim=range(0:total.y.size), main = main) 
  
  m <- 1
  for(i in 1:x.n){
    for(j in 1:y.n){
      rect(xleft = (i-1)*x.length, ybottom = (j-1)*y.length, xright = i*x.length, ytop = j*y.length, density = NULL, angle = 45,
           col = col[m], border = col[m], lty = par("lty"), lwd = par("lwd"))
      m <- m+1
    }
  }
  
  return()
}
  

###Set up
#TMB stuff
compile("TMBPoisson.cpp")
dyn.load(dynlib("TMBPoisson"))

compile("TMBPoissonNoPixels.cpp")
dyn.load(dynlib("TMBPoissonNoPixels"))

#choose size
total.x.size <- 70
total.y.size <- 70

#choose number of pixels
x.n <- 70
y.n <- 70

x.length <- total.x.size/x.n
y.length <- total.y.size/y.n

n.total <- x.n*y.n


#generate coordinates
x <- seq(0, total.x.size, length.out=x.n)
y <- seq(0, total.y.size, length.out=y.n)

m <- 1
coord <- array(0,c(n.total,2))
for(i in 1:x.n){
  for(j in 1:y.n){
    coord[m,] = c(x[i], y[j])
    m <- m+1
  }
}

#create coordinates for border 
#(actually need just the corners, so border.n=1 and this is unecessary, it turns out...)
border.n <- 1
border <- array(0, c(4*(border.n+1), 2))
for(i in 1:(border.n + 1)){
  border[i,] <- c(0, (i-1)*total.y.size/border.n)
  border[i + (border.n+1), ] <- c(total.x.size, (i-1)*total.y.size/border.n)
  border[i + 2*(border.n+1), ] <- c((i-1)*total.x.size/border.n, 0)
  border[i + 3*(border.n+1), ] <- c((i-1)*total.x.size/border.n, total.y.size)
}

plot(border)
plot(coord)

###Generate fake data
#generate random field
a <- RMmatern(3)
as <- RFsimulate(a, x=coord[,1], y=coord[,2])
a.field <- as.matrix(as)

a2 <- RMmatern(4, scale=3, var=1)
as2 <- RFsimulate(a2, x=coord[,1], y=coord[,2])
cov <- as.matrix(as2)

a3 <- RMmatern(0.2, scale=5, var=6)
as3 <- RFsimulate(a3, x=coord[,1], y=coord[,2])
cov2 <- as.matrix(as3)


#visualise random field
pixelplot(c, total.x.size, total.y.size, x.n, y.n, main="Matern RF (Random Fields)")

plot(coord, cex = c, main = "Matern RF (Random Fields)")

#another way to generate random field
if(FALSE){
library(geoR)

d <- grf(n, grid = coord, xlims = c(0, 10), ylims = c(0, 10), nsim = 1, cov.model = "matern",
        cov.pars = c(1,1), 
        kappa = 0.5, nugget = 0, lambda = 1.00,
        mean = 0, RF=TRUE)

plot(coord, cex=d$data, main = "Matern RF (geoR)")
}

#generate covariate and responses
#cov <- rep(0, n.total)
response <- rep(0, n.total)


for(i in 1:n.total){
  #cov[i] <- rnorm(1,coord[i,1]+coord[i,2], 20)/20
  #cov[i] <- runif(1,0,3)
  response[i] <- rpois(1, exp(b0 + b1*cov[i] + b2*cov2[i] + a.field[i]))
}

#visualise response
pixelplot(response, total.x.size, total.y.size, x.n, y.n, main="Response")


###Create spde object
mesh <- inla.mesh.2d(loc.domain = border, max.edge=c(3,2), offset=c(0.03, 0.5), cutoff=1, max.n = 500) 
plot(mesh)

spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]

idx <- mesh$idx$loc
n_s <- nrow(spde$M0)

x <- rep(0, n_s)

###Create "box" (i.e. aggregate) data
box.x.n <- 10
box.y.n <- 10
box.n.total <- box.x.n*box.y.n

box.x.length <- total.x.size/box.x.n
box.y.length <- total.y.size/box.y.n

box.cov <- rep(0, box.n.total)
box.response <- rep(0, box.n.total)
box.total <- rep(0, box.n.total)
box.index <- rep(c(), box.n.total)

coord.box.number <- rep(0, n.total)

#for each pixel, find which box it's in, then sum the pixel covariates and responses for each box

for(i in 1:n.total){
  box.number <- boxfind(coord[i,], box.x.length, box.y.length, box.x.n)
  coord.box.number[i] <- box.number
  box.cov[box.number] <- box.cov[box.number] + cov[i]
  box.response[box.number] <- box.response[box.number] + response[i]
  box.total[box.number] <- box.total[box.number] + 1
}

ordered.coord <- array(0, c(n.total, 2))
ordered.index <- rep(0, n.total)
ordered.cov <- rep(0, n.total)
ordered.cov2 <- rep(0, n.total)

m <- 1
for(i in 1:box.n.total){
  for(j in 1:n.total){
    if(coord.box.number[j] == i){
      ordered.coord[m, ] <- coord[j, ]
      ordered.cov[m] <- cov[j]
      ordered.cov2[m] <- cov2[j]
      ordered.index[m] <- j
      m <- m+1
    }
  }
}

ordered.index.inverse <- rep(0, n.total)

for(i in 1:n.total){
  ordered.index.inverse[ordered.index[i]] = i
}



#get an average value for covariates and responses - not sure we actually want this
for(i in 1:box.n.total){
  box.index[[i]] <- rep(0, box.total[i])
  #box.cov[i] = box.cov[i]/box.total[i]
  #box.response[i] = box.response[i]/box.total[i]
}

box.count <- rep(1, box.n.total)
for(i in 1:n.total){
  box.number <- boxfind(coord[i,], box.x.length, box.y.length, box.y.n)
  box.index[[box.number]][box.count[box.number]] = i
  box.count[box.number] = box.count[box.number] + 1
}

#visualise aggregated data
pixelplot(box.cov, total.x.size, total.y.size, box.x.n, box.y.n, main="Aggregated (averaged) covariate")
pixelplot(box.response, total.x.size, total.y.size, box.x.n, box.y.n, main="Aggregated (averaged) response")

A <- inla.spde.make.A(mesh=mesh, loc=ordered.coord) 

if(TRUE){

f <- MakeADFun(
      data = list(X=box.response, cov=ordered.cov, cov2 = ordered.cov2, spde=spde, Apixel = A, box_total = box.total),
      parameters = list(beta0=0, beta1=0, beta2=0, log_kappa=2.5, log_tau=0.0, x=runif(n_s,0,10)),
      random="x",
      DLL = "TMBPoisson"
)

#fit <- nlminb(f$par,f$fn,f$gr,lower=c(-10,-10,0,0))
fit.box <- nlminb(f$par,f$fn,f$gr)


if(FALSE){
f3 <- MakeADFun(
  data = list(X=box.response, cov=cov, cov2 = cov2, spde=spde, Apixel = A, box_total = box.total),
  parameters = list(beta0=0, beta1=0, beta2=0, log_kappa=2.5, log_tau=0.0, x=runif(n_s,0,10)),
  random="x",
  DLL = "TMBPoisson"
)



#fit <- nlminb(f$par,f$fn,f$gr,lower=c(-10,-10,0,0))
fit.box.wrong <- nlminb(f3$par,f3$fn,f3$gr)
}

A2 <- inla.spde.make.A(mesh=mesh, loc=coord)  

f2 <- MakeADFun(
  data = list(X=response, cov=cov, cov2 = cov2, spde=spde, Apixel = A2),
  parameters = list(beta0=0, beta1=0, beta2=0, log_kappa=2.5, log_tau=0.0, x=runif(n_s,0,10)),
  random="x",
  DLL = "TMBPoissonNoPixels"
)


#fit <- nlminb(f$par,f$fn,f$gr,lower=c(-10,-10,0,0))
fit.points <- nlminb(f2$par,f2$fn,f2$gr)

#print(fit.box.wrong$par)
print(fit.box$par)
print(fit.points$par)
cat("beta0", b0, "beta1", b1, "beta2", b2,"\n")


}

if(FALSE){

if(FALSE){
compile("TMBExample.cpp")
dyn.load(dynlib("TMBExample"))

g <- MakeADFun(
  data = list(x=response, cov=cov),
  parameters = list(beta0=0, beta1=0, sigma=1),
  DLL = "TMBExample"
)

fit2 <- nlminb(g$par,g$fn,g$gr,lower=c(-10,-10,0))


print(fit2)
}

pred <- rep(0, n.total)
beta0 <- fit$par[1]
beta1 <- fit$par[2]

error1 = 0
for(i in 1:n.total){
  pred[i] <- rpois(1, exp(beta0 + beta1*cov[i]))
  error1 = error1 + abs(pred[i] - response[i])
}

pixelplot(response, total.x.size, total.y.size, x.n, y.n, main="Response")
pixelplot(pred, total.x.size, total.y.size, x.n, y.n, main="Predicted Response")

####

field <- rep(0, n.total)
out <- sdreport(f, getJointPrecision = 1)
a = array(0, c(1035, 1))
for(i in 1:1035){
  a[i] <- out$par.random[i]
}

error = 0
for(i in 1:n.total){
  pred[i] <- exp(beta0 + beta1*cov[i] + (A%*%a)[ordered.index.inverse[i]])
  field[i] <- (A%*%a)[ordered.index.inverse[i]]
  error = error + abs(pred[i] - response[i])
}


pixelplot(response, total.x.size, total.y.size, x.n, y.n, main="Response")
pixelplot(pred, total.x.size, total.y.size, x.n, y.n, main="Predicted Response")
pixelplot(field, total.x.size, total.y.size, x.n, y.n, main="Field")

pixelplot(c, total.x.size, total.y.size, x.n, y.n, main="Field")


##confidence intervals...
b <- summary(out, "fixed")
beta0max <- b[1,1] + 2*b[1,2]
beta0min <- b[1,1] - 2*b[1,2]
beta1max <- b[2,1] + 2*b[2,2]
beta1min <- b[2,1] - 2*b[2,2]

count <- 0
predmax <- rep(0, n.total)
predmin <- rep(0, n.total)
for(i in 1:n.total){
  predmax[i] <- exp(beta0max + beta1max*cov[i] + (A%*%a)[ordered.index.inverse[i]])
  predmin[i] <- exp(beta0min + beta1min*cov[i] + (A%*%a)[ordered.index.inverse[i]])
  
  if(response[i]<predmax[i] && response[i]>predmin[i]){
    count <- count+1
  }
  
}

}

