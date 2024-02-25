#' ---
#' title: Diffusion equation
#' author: Leonor
#' date: 30 January 2024
#' ---

##### Packages ####
library(deSolve)
library(MetBrewer)
library(ggplot2)
library(fields)

n <- 100
d <- 100
param.deltaZ <- d/n
t <- 100
param.u <- 0.04* 24 #[/day]
param.D <- 43.2 #[m^2/day]
z <- c()
param.kw <- 0.2 #[/m]
param.kp <- 15*10^(-12) #[m^2 / cell]
param.Io <- 350 * 60 * 60 * 24 #[µmol photons /m^2 /day]
param.HI <- 30 * 60 * 60 *24 # [µmol photons / m^2 / day]
param.gmax <-0.04* 24
param.m <- 0.01* 24 #[/day]
param.HN <- 0.0425 # [mmol nutrient / m^3]
param.Nb <- 5 # [mmol nutrient / m^3]
param.Y <- 10^(-9) # [mmol nutrient / cell]

#Defining the grid

z = seq(param.deltaZ/2, by=param.deltaZ, to=d)

#Defining phi

#phi <- seq(1,n,1)
phi <- rep(1000000,n)
#phi[1] <- 10

#P <- phi[1]

#Defining N

N <- rep(1,n)
#N[1] <-100

Y <- c(phi, N)

#Integral of light intensity


func.I = function (param, P){
  integral = cumsum(param.kw + param.kp*P)*param.deltaZ - 0.5*param.deltaZ*(param.kw + param.kp*P)
  I <- param.Io*exp(- integral)
  
  return(I)
}

#I.out <- func.I(param,phi)

#plot(I.out,z)

#Growth function

func.growth = function(param,I, N){
  g.I <- I/(I+param.HI)
  g.N <- N/(N+param.HN) 
  return(c(g.I,g.N))
}


#Defining the flux

derivative = function(t, Y, params) {
  
  phi <- Y[1:n]
  N <- Y[(n+1):(2*n)]
  
  Ja <- rep(0, n+1)
  Jd <- rep(0, n+1)
  JN <- rep(0,n+1)
  
  for(i in 2:n){
    Ja[i] <- param.u*phi[i-1]
    Jd[i] <- -param.D*(phi[i]-phi[i-1])/param.deltaZ
    JN[i] <- -param.D*(N[i]-N[i-1])/param.deltaZ 
  }

#Boundary conditions  
  Ja[1] <- 0
  Ja[n+1] <- 0
  Jd[1] <- 0
  Jd[n+1] <- 0
  JN[n+1] <- -param.D*(param.Nb-N[n])/param.deltaZ
  JN[1] <- 0
  
  J <- Ja + Jd 
  
  dphidt <- seq (1,n,1)
  dNdt <- seq (1,n,1)
  
  # Calculate I
  I = func.I(param,phi)
  
  
# calculate r
  #g = func.growth(param, I, N)
  g = param.gmax*pmin(I/(I+param.HI), N/(N+param.HN))
  r = g - param.m
  
  
  for (i in 1:n){
    dphidt[i] <- r[i]*phi[i]- (J[i+1]- J[i])/param.deltaZ
    dNdt[i] <- -param.Y*r[i]*phi[i] - (JN[i+1]- JN[i])/param.deltaZ
  }
  
  return( list(c(dphidt, dNdt)) )
}


times <- seq(0,1000,1)  

derivative.out <- ode(Y, times, derivative)
derivative.out

phi.out <- derivative.out[,2:(n+1)]
N.out <- derivative.out[,(n+2):(2*n+1)]

I.out <- func.I(param, phi.out[length(times),]) 
g.out <- func.growth(param, I.out, N.out[length(times),])

#plot(I.out[100,],z, type ="l", lty = 1)

par(mar=c(5, 4, 4, 8), xpd=TRUE) # adds space next to the plot so we can put a legend there
matplot(cbind(g.out[n:1],g.out[(2*n):(n+1)]),z, type = "l", lty = 1, 
        col = c("red", "blue"), xlab="Growth rate", ylab= "Distance to seabed", main="Limiting factor")
legend("topright", inset=c(-0.5,0), legend = c("Light", "Nutrient"), 
       col = c("red", "blue"),
       lty = 1)


time = derivative.out[,1]
#ph = derivative.out[,ncol(derivative.out):2]

image.plot(time,z,phi.out[,ncol(phi.out):1],col = hcl.colors(50, "viridis"), ylab="Distance from seabed", main="Phytoplankton")
image.plot(time,z,N.out[,ncol(phi.out):1],col = hcl.colors(50, "viridis"), ylab="Distance from seabed", main= "Nutrients")
