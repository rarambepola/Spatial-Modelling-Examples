library(RandomFields)

#choose number of points
n=100


#choose pixel size or pixel number
size = 0.1
#pixel.n = 10/size
pixel.n = 100

#generate pixe; coordinates and covariates
x <-seq(0,10,length.out=pixel.n)
y <-seq(0,10,length.out=pixel.n)

pixel.coord <- array(0, c(pixel.n**2, 2))

m = 1
for(i in 1:pixel.n){
  for(j in 1:pixel.n){
    pixel.coord[m,] <- c(x[i], y[j])
    m = m+1
  }
}

pixel.cov <- runif(pixel.n**2,0,3)



#generate coordinates for other points
coord <- array(0,c(n,2))
for(i in 1:n){
  coord[i,] <- c(runif(1,0,10), runif(1,0,10))
}

#generate matern random field
a = RMmatern(0.5)
b = RFsimulate(a, x=coord[,1], y=coord[,2])

c = as.matrix(b)

#visualise data so far
plot(coord, cex = c, main = "Matern RF (Random Fields)")
plot(coord, main="Location of data points")
points(pixel.coord, col="red")

#generate covariates for data points and responses
cov <- rep(0, n)
response <- rep(0, n)


for(i in 1:n){
  cov[i] <- runif(1,0,3)
  response[i] <- rnorm(1, 4 + 2*cov[i] + d$data[i])
}

#set up INLA spde object
library(INLA)


mesh <- inla.mesh.2d(loc = coord, max.edge=c(35,100)) 

spde <- inla.spde2.matern(mesh=mesh, alpha=2) 

A <- inla.spde.make.A(mesh=mesh, loc=coord) 

A.pred <- inla.spde.make.A(mesh=mesh, loc=pixel.coord)

s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)

#add data to stack
stack.norm <- inla.stack(data=list(Response=response), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                                 list(Covariate=cov)), tag="all")

stack.pred <- inla.stack(data=list(Response=NA), A=list(A.pred,1), effects=list(c(s.index, list(Intercept=1)),
                                                                                   list(Covariate=pixel.cov)), tag="pred")
#define formula for model
formula <- Response ~ -1 + Intercept + Covariate + f(spatial.field, model=spde)
#formula <- Response ~ -1 + Intercept + f(spatial.field, model=spde)

stack.join <- inla.stack(stack.norm, stack.pred)

output.pred <- inla(formula, data=inla.stack.data(stack.join, spde=spde), family="gaussian", 
                    control.predictor=list(A=inla.stack.A(stack.join), compute=TRUE),
                    control.compute=list(cpo=TRUE, dic=TRUE))

output.pred$summary.fixed

#find the pixel point indices
index <- inla.stack.index(stack.join, tag='pred')$data

#extract the predictions for nu at the pixels and generate predictions for the response
nu.pred <- output.pred$summary.linear.pred[index,]$mean

response.pred = array(0,c(length(nu.pred),1))
for(i in 1:length(nu.pred)){
  response.pred[i] = rnorm(1, nu.pred[i], 0.1)
}

#visualise predicted responses against real responses
plot(coord, cex=response/(r<-7), main = "Responses")
#points(pixel.coord, cex=response.pred/r, col="red")

#visualise the predicted responses as a heatmap
testm <- array(0, c(pixel.n, pixel.n))
for(i in 1:length(response.pred)){
  testm[i] <- response.pred[i]
}

library(fields)
image.plot(testm, col=two.colors(n=256, start="white", end="black", middle="grey" ))
