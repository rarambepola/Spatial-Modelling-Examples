boxfinder <- function(coords, edgelength){
  n <- 1
  for(i in 1:sqrt(nboxes)){
    if(coords[1]>i*edgelength){
      n <- n+1
    }
  }
  
  m <- 1 
  for(i in 1:sqrt(nboxes)){
    if(coords[2]>i*edgelength){
      m <- m+1
    }
  }
  return(n + (m-1)*4)
}

n=15

x <- runif(n*n, 0, 10)
y <- runif(n*n, 0, 10)

m=1
coord <- array(0,c(n*n,2))
for(i in 1:n){
  for(j in 1:n){
    coord[m,] <- c(runif(1,0,10),runif(1,0,10))
    m=m+1
  }
}



w <- rep(0, n*n)


for(i in 1:(n*n)){
  w[i] <- rnorm(1,5*coord[i,1]+5*coord[i,2],10)
}

plot(coord)
rbPal <- colorRampPalette(c('red','blue'))
w.col <- rbPal(10)[as.numeric(cut(w,breaks = 10))]
plot(coord, col=w.col, main="W", pch=19)

#plot(w,coord[,1])
#plot(w,coord[,2])

b.index <- c()
wh.index <- c()
b.count <- 1 
wh.count <- 1

b.yesno <- rep(0,n*n)

for(i in 1:(n*n)){
  if(runif(1,min(w)-1,w[i])<15){
    b.index[[b.count]] <- i
    b.count <- b.count + 1
    b.yesno[i] <- 1
  }
  else{
    wh.index[[wh.count]] <- i
    wh.count <- wh.count + 1
  }
}

b <- coord[b.index,]
wh <- coord[wh.index,]

plot(b, col="black", main="", pch=19)
points(wh, col="darkorange", main="wh", pch=19)


rlevel <- rep(0, n*n)

for(i in 1:(n*n)){
  rlevel[i] <- rnorm(1, 10+ 5*w[i], 1)
  
}

rbPal <- colorRampPalette(c('yellow','black'))
rlevel.col <- rbPal(10)[as.numeric(cut(rlevel,breaks = 10))]
plot(coord, col=rlevel.col, main="rlevel", pch=19)

sqrtnumberofboxes = 7
nboxes = sqrtnumberofboxes**2

edgelength = 10/sqrtnumberofboxes

boxcoords <- array(0,c(nboxes,2))
box.count <- 1
box.rlevel <- rep(0, nboxes)
box.number <- rep(0, nboxes)
box.wealth <- rep(0, nboxes)
box.totalnumber <- rep(0, nboxes)

for(i in 1:(sqrtnumberofboxes)){
  for(j in 1:sqrtnumberofboxes){
  boxcoords[box.count,] = c(i*edgelength, j*edgelength)
  box.count <- box.count + 1}
}


#points(boxcoords, col = "red", pch=19)

for(i in 1:(n*n)){
  boxnum <- boxfinder(coord[i,], edgelength)
  box.rlevel[boxnum] <- box.rlevel[boxnum] + rlevel[i]
  box.number[boxnum] <- box.number[boxnum] + b.yesno[i]
  box.totalnumber[boxnum] <- box.totalnumber[boxnum] + 1
  box.wealth[boxnum] <- box.wealth[boxnum] + w[i]
}

for(i in 1:nboxes){
  box.rlevel[i] = box.rlevel[i]/box.totalnumber[i]
  box.number[i] = box.number[i]/box.totalnumber[i]
  box.wealth[i] = box.wealth[i]/box.totalnumber[i]
}

library(INLA)
mesh <- inla.mesh.2d(loc = coord, max.edge=c(35,100)) 

spde <- inla.spde2.matern(mesh=mesh, alpha=2) 


A <- inla.spde.make.A(mesh=mesh, loc=boxcoords) 

s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde)


stack.norm <- inla.stack(data=list(rlevel=box.rlevel), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                         list(Number=box.number)), tag="all")

formula <- rlevel ~ -1 + Intercept + Number + f(spatial.field, model=spde)

#formula <- rlevel ~ -1 + Intercept + Number 


Apred <- inla.spde.make.A(mesh, loc=coord)
stack.pred <- inla.stack(data=list(rlevel=NA), A=list(Apred,1), effects=list(c(s.index, list(Intercept=1)),
                                                                                 list(Number=b.yesno)), tag="pred")

stack.join <- inla.stack(stack.norm, stack.pred)

output.pred <- inla(formula, data=inla.stack.data(stack.join, spde=spde), family="gaussian", 
                          control.predictor=list(A=inla.stack.A(stack.join), compute=TRUE),
                          control.compute=list(cpo=TRUE, dic=TRUE))

output.pred$summary.fixed


index <- inla.stack.index(stack.join, tag='pred')$data


nu.pred = output.pred$summary.linear.pred[index,]$mean
rlevel.pred = array(0,c(length(nu.pred),1))
for(i in 1:length(nu.pred)){
  rlevel.pred[i] = rnorm(1, nu.pred[i], 1)
}

#plot(b, col="black", main="", pch=19)
#points(wh, col="darkorange", main="wh", pch=19)


rbPal <- colorRampPalette(c('yellow','black'))
rlevel.pred.col <- rbPal(10)[as.numeric(cut(rlevel.pred,breaks = 10))]
plot(coord, col=rlevel.pred.col, main="Predicted rlevel", pch=19)

