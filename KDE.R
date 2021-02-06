#KDE in the multivariate case
library(ks)
library(viridis)
library(MASS)
library(plotly)
library(bivariate) # produces the density functions
library(tidyverse)
library(ggExtra) #compare marginal densities

#Study parameters
set.seed(123)

#Preparation of study

#Study parameters
#Number of draws per sample
n = 100

#Number of repetitions in simulation
rep = 1000

#Number of dimensions of data
dim = 2

#Create an empty data array with the correct dimensions
#1. dimension: draws, 2. dimension: variables, 3. dimension: simulations
data = array(data = rep(NA, n * dim * rep), c(n, dim, rep))

### Create bivariate density ###
#bimodal bivariate normal mixture
bi_mix <- bmbvpdf(mean.X1 = 0, mean.Y1 = -3, sd.X1 = 1, sd.Y1 = 1.5, mean.X2 = -3, mean.Y2 = 3, sd.X2 = 1.2, sd.Y2 = 2.1)
plot(bi_mix)

### Sample from the function ###
#Acceptance-Rejection sampling
#Bivariate uniform with maximum of target
x1_low <- -10
x1_up <- 10
x2_low <- -10
x2_up <- 10
upper <- cubvpdf(a.X = x1_low, b.X = x1_up, a.Y = x2_low, b.Y = x2_up) # interval

#Add constant
sampled <- data.frame(x = runif(100000, x1_low, x1_up), y = runif(100000,x2_low,x2_up))
sampled$density <- bi_mix(x = sampled$x, y = sampled$y)

max_value <- max(sampled$density) #Take this or a slightly higher value as constant to add to the uniform
max_value <- max_value * 1.1 #Increase by 10% to make sure, we actually are higher then the peak

rm(sampled)

#Sampling function
ARS <- function(target, approx, x_low, x_up, y_low, y_up, n, c){
  output <- data.frame(x=rep(NA,n),y=rep(NA, n))
  counter <- 1
  while(counter < (n+1)) {
    x_test <- runif(1, x_low, x_up)
    y_test <- runif(1, y_low, y_up)
    z_test <- target(x_test, y_test)
    g_z <- c * approx(x_test, y_test)
    if(g_z == 0) {
      next
    }
    ratio <- z_test/g_z
    u <- runif(1)
    if(u <= ratio){
      output$x[counter] <- x_test
      output$y[counter] <- y_test
      counter <- counter + 1
    }
  }
  return(output)
}

for (i in 1:rep){
  sample <- ARS(target = bi_mix, approx = upper, x_low = x1_low,
                x_up = x1_up, y_low = x2_low, y_up = x2_up, n = n, c = max_value)
  data[,1,i] = sample$x
  data[,2,i] = sample$y
}

#Prepare matrix where estimates will be saved
#First dimension: x-coordinates
#Second dimension: y-coordinates
#Kernel density estimates will be values
#Third dimension: simulations

#Grid resolution
nx = 501
ny = nx
xmin = -10
ymin = -10
xmax = 10
ymax = 10
xcoordinates = seq(xmin, xmax, length.out = nx)
ycoordinates = seq(ymin, ymax, length.out = ny)

#Save the Kernel density estimates
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))
#Kernel density estimation
for (i in 1:rep){
  Hpi <- ks::Hpi(x = data[,,i])
  kdeHpi <- ks::kde(x = data[,,i], H = Hpi,
                    xmin=c(xmin,ymin), xmax=c(xmax,ymax),
                    gridsize = c(nx,ny), bgridsize=c(nx,ny))
  estimates[,,i] = kdeHpi$estimate
}
#Average estimate
aver_est = array(data = rep(NA, nx * ny), c(nx, ny),
                 dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx) {
  for (j in 1:ny){
    aver_est[i,j] = mean(estimates[i,j,])
  }
}

#Theoretical density values on grid points
true_values = array(data = rep(NA, nx * ny), c(nx, ny),
                    dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

#Calculate true density values
for (i in 1:nx) {
  xvalue = xcoordinates[i]
  for (j in 1:ny){
    yvalue = ycoordinates[j]
    true_values[i,j] = bi_mix(xvalue, yvalue)
  }
}

#Difference of true values and estimates
diff_est_true = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                      dimnames = list(as.character(xcoordinates),
                                      as.character(ycoordinates), as.character(c(1:rep))))
for (i in 1:rep){
  diff_est_true[,,i] = estimates[,,i] - true_values[,]
}

#Equalize average difference and bias
bias_estimates = aver_est - true_values

#Calculate relative bias
rel_bias_estimates = bias_estimates / true_values

#Calculate variances of estimates
var_estimates = array(data = rep(NA, nx * ny), c(nx, ny),
                      dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx){
  for (j in 1:ny){
    var_estimates[i,j] = var(estimates[i,j,])
  }
}

#RMSE of estimates
rmse_estimates = array(data = rep(NA, nx * ny), c(nx, ny),
                       dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))
for (i in 1:nx){
  for (j in 1:ny){
    rmse_estimates[i,j] = sqrt(mean(diff_est_true[i,j,]^2))
  }
}
