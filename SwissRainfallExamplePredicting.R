

library(geoR) #package containing the rainfall data
library(INLA)

#rescale coordinates
sic.all$coords[,1] <- (sic.all$coords[,1] - apply(sic.borders,2,mean)[1])
sic.all$coords[,2] <- (sic.all$coords[,2] - apply(sic.borders,2,mean)[2])
sic.borders <- apply(sic.borders,2,scale,scale=F)

#extract coordinates, rainfall data and elevation data
coord <- sic.all$coords
datn <- sqrt(sic.all$data)
datp <- sqrt(sic.all$data)
elevation <- sic.all$altitude/1000

#create fake data based on elevation
for(i in 1:length(elevation)){
    datn[i] <- rnorm(1, 1+ 12*elevation[i],1)
    datp[i] <- rpois(1, exp(0.5 + 2*elevation[i]))
}

#create mesh, spde object, projector matrix etc.
Swiss.mesh <- inla.mesh.2d(loc.domain=sic.borders, max.edge=c(35,100)) 

Swiss.spde <- inla.spde2.matern(mesh=Swiss.mesh, alpha=2) 

A <- inla.spde.make.A(mesh=Swiss.mesh, loc=coord) 

s.index <- inla.spde.make.index(name="spatial.field", n.spde=Swiss.spde$n.spde)

#create data stack for normal and poisson examples
stack.norm <- inla.stack(data=list(rain=datn), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                    list(Elevation=elevation)), tag="all")


stack.pois <- inla.stack(data=list(rain=datp), A=list(A,1), effects=list(c(s.index, list(Intercept=1)),
                                                                       list(Elevation=elevation)), tag="all")

#specify linear model
formula <- rain ~ -1 + Intercept + Elevation + f(spatial.field, model=spde)

#run estimation
Swiss.output.pois <- inla(formula, data=inla.stack.data(stack.pois, spde=Swiss.spde), family="poisson", 
                     control.predictor=list(A=inla.stack.A(stack.pois), compute=TRUE),
                     control.compute=list(cpo=TRUE, dic=TRUE))

Swiss.output.norm <- inla(formula, data=inla.stack.data(stack.norm, spde=Swiss.spde), family="gaussian", 
                         control.predictor=list(A=inla.stack.A(stack.norm), compute=TRUE),
                         control.compute=list(cpo=TRUE, dic=TRUE))

Swiss.output.pois$summary.fixed
Swiss.output.norm$summary.fixed

plot(elevation, datn)



if(FALSE){
plot(coord, cex=elevation, main="Elevation")
#plot(coord, cex=datp/40, main = "Poisson response")
plot(coord, cex=datn/20, main = "Normal response")
  
rbPal <- colorRampPalette(c('red','blue'))
elevation.col <- rbPal(10)[as.numeric(cut(elevation,breaks = 10))]
plot(coord, col=elevation.col, main="Elevation", pch=19)
#points(coord, pch = 2, cex = datn/20)


datn.col <- rbPal(10)[as.numeric(cut(max(elevation)*datn/max(datn),breaks = 10))]
plot(coord, col=datn.col, main="Normal response", pch=19)


datp.col <- rbPal(10)[as.numeric(cut(datp,breaks = 10))]
plot(coord, col=datp.col, main="Poisson response", pch=19)
}


