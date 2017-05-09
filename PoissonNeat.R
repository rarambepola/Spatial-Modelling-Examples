setwd("J:/INLAExamples")

library(RandomFields)
library(TMB)
library(INLA)

##coefficients and pixel/box numbers
b0 <- 1.0
b1 <- 0.7
b2 <- 1.0

#choose number of pixels
x.n <- 70
y.n <- 70

#choose number of boxes
box.x.n <- 10
box.y.n <- 10

#choose size of grid
total.x.size <- 70
total.y.size <- 70

###Set up
#TMB stuff
compile("TMBPoisson.cpp")
dyn.load(dynlib("TMBPoisson"))

compile("TMBPoissonNoPixels.cpp")
dyn.load(dynlib("TMBPoissonNoPixels"))


##Functions
#function to find which box number a point is in
#(boxes are numbered left to right, bottom to top, e.g. lowest row
#is 1,2,3,4 next is 5,6,7,8 etc. )
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


#generate coordinates
x.length <- total.x.size/x.n
y.length <- total.y.size/y.n

n.total <- x.n*y.n

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


###Generate fake data
#generate random field
a <- RMmatern(3)
as <- RFsimulate(a, x=coord[,1], y=coord[,2])
a.field <- as.matrix(as)


#generate covariate and responses
a2 <- RMmatern(4, scale=3, var=1)
as2 <- RFsimulate(a2, x=coord[,1], y=coord[,2])
cov <- as.matrix(as2)

a3 <- RMmatern(0.2, scale=5, var=6)
as3 <- RFsimulate(a3, x=coord[,1], y=coord[,2])
cov2 <- as.matrix(as3)

response <- rep(0, n.total)
for(i in 1:n.total){
  response[i] <- rpois(1, exp(b0 + b1*cov[i] + b2*cov2[i] + a.field[i]))
}


###Create spde object
mesh <- inla.mesh.2d(loc.domain = border, max.edge=c(3,2), offset=c(0.03, 0.5), cutoff=1, max.n = 500) 


spde <- (inla.spde2.matern(mesh=mesh, alpha=2)$param.inla)[c("M0","M1","M2")]

idx <- mesh$idx$loc
n_s <- nrow(spde$M0)

x <- rep(0, n_s)

###Create "box" (i.e. aggregate) data
box.n.total <- box.x.n*box.y.n

box.x.length <- total.x.size/box.x.n
box.y.length <- total.y.size/box.y.n

box.response <- rep(0, box.n.total)
box.total <- rep(0, box.n.total)
box.index <- rep(c(), box.n.total)

coord.box.number <- rep(0, n.total)

#for each pixel, find which box it's in, then sum the pixel covariates and responses for each box

for(i in 1:n.total){
  box.number <- boxfind(coord[i,], box.x.length, box.y.length, box.x.n)
  coord.box.number[i] <- box.number
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

box.count <- rep(1, box.n.total)
for(i in 1:n.total){
  box.number <- boxfind(coord[i,], box.x.length, box.y.length, box.y.n)
  box.count[box.number] = box.count[box.number] + 1
}

A <- inla.spde.make.A(mesh=mesh, loc=ordered.coord) 



f <- MakeADFun(
      data = list(X=box.response, cov=ordered.cov, cov2 = ordered.cov2, spde=spde, Apixel = A, box_total = box.total),
      parameters = list(beta0=0, beta1=0, beta2=0, log_kappa=2.5, log_tau=0.0, x=runif(n_s,0,10)),
      random="x",
      DLL = "TMBPoisson"
)

#fit <- nlminb(f$par,f$fn,f$gr,lower=c(-10,-10,0,0))
fit.box <- nlminb(f$par,f$fn,f$gr)


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


