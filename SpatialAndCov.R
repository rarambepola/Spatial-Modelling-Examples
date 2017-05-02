library(INLA)

#choose number of points + prediction points
n=250
pred.n = 100

#generate coordinates for points
coord <- array(0,c(n,2))
for(i in 1:n){
  coord[i,] <- c(runif(1,0,10),runif(1,0,10))
}


pred.coord <- array(0,c(pred.n,2))
for(i in 1:pred.n){
  pred.coord[i,] <- c(runif(1,0,10), runif(1,0,10))
}

#initialize arrays
w <- rep(0, n)
cov <- rep(0, n)
response <- rep(0, n)

#generate 'spatial covariate' - should be replaced by actual random field later
for(i in 1:n){
  w[i] <- rnorm(1,(5*coord[i,1]+5*coord[i,2])**2,10)
}

pred.w <- rep(0, pred.n)

for(i in 1:pred.n){
  pred.w[i] <- rnorm(1,(5*pred.coord[i,1]+5*pred.coord[i,2])**4,10)
}

#generate some random values for a covariate
cov <- rep(0, n)
for(i in 1:(n)){
  cov[i] <- runif(1,0,100)
}

pred.cov <- rep(0, pred.n)

for(i in 1:pred.n){
  pred.cov[i] <- runif(1,0,100)
}

#plot the points and visualise the values of the 'random field' w (colour) and the covariate (size)
rbPal <- colorRampPalette(c('red','blue'))
w.col <- rbPal(10)[as.numeric(cut(w,breaks = 10))]
pred.w.col <- rbPal(10)[as.numeric(cut(pred.w,breaks = 10))]
plot(coord, col=w.col, main="W (colour), cov(size)", pch=19, cex=cov/50)
points(pred.coord, col=pred.w.col, pch=17, cex=pred.cov/50)


#generate the response based on the field and the covariate and visualise
for(i in 1:n){
  response[i] <- rnorm(1, 10*cov[i] + 0.07*w[i],1)
}

plot(coord, main="Response", cex=(resize <- 3)*response/max(response))


#create mesh, spde, projector matrix etc.
mesh <- inla.mesh.2d(loc = coord, max.edge=c(35,100)) 

spde <- inla.spde2.matern(mesh=mesh, alpha=2) 


A <- inla.spde.make.A(mesh=mesh, loc=coord) 

s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)

#add data to stack
stack.norm <- inla.stack(data=list(Response=response), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                                 list(Covariate=cov)), tag="all")
#define formula for model
formula <- Response ~ -1 + Intercept + Covariate + f(spatial.field, model=spde)

#repeat with the prediction points
Apred <- inla.spde.make.A(mesh, loc=pred.coord)


stack.pred <- inla.stack(data=list(Response=NA), A=list(Apred,1), effects=list(c(s.index, list(Intercept=1)),
                                                                             list(Covariate=pred.cov)), tag="pred")


stack.join <- inla.stack(stack.norm, stack.pred)


#fit the model using INLA
output.pred <- inla(formula, data=inla.stack.data(stack.join, spde=spde), family="gaussian", 
                    control.predictor=list(A=inla.stack.A(stack.join), compute=TRUE),
                    control.compute=list(cpo=TRUE, dic=TRUE))

#print the output - check if the coefficients are correct
output.pred$summary.fixed

#find out index for predction points and generate responses based on the predictions of nu
index <- inla.stack.index(stack.join, tag='pred')$data

nu.pred = output.pred$summary.linear.pred[index,]$mean
response.pred = array(0,c(length(nu.pred),1))
for(i in 1:length(nu.pred)){
  response.pred[i] = rnorm(1, nu.pred[i], 1)
}

plot(pred.coord, col=pred.w.col, pch=17, cex=pred.cov/50)
plot(pred.coord, pch=2, cex = resize*response.pred/max(response.pred))

