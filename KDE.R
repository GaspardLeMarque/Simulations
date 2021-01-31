#KDE in the multivariate case
library(ks)
library(viridis)
library(rgl)
library(misc3d)
library(MASS)

#Study parameters
#Number of draws per sample
n = 100

#Number of repetitions in simulation
rep = 100

#Number of dimensions of data
dim = 2

#Create an empty data array with the given parameters
#1. dimension: draws, 2. dimension: variables, 3. dimension: simulations
data = array(data = rep(NA, n * dim * rep), c(n, dim, rep))

#Bivariate normal distribution
#Data generation process
mu = c(0,0)
sigma = matrix(data = c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
for (i in 1:rep){
  sample = MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  data[,,i] = sample
}

#Prepare matrix where estimates will be saved
#First dimension: x-coordinates
#Second dimension: y-coordinates
#Values will be kernel density estimates
#Third dimension: simulations

#Grid resolution
nx = 151
ny = 151
xmin = -10
ymin = -10
xmax = 10
ymax = 10
xcoordinates = seq(xmin, xmax, length.out = 151)
ycoordinates = seq(ymin, ymax, length.out = 151)

#Save the Kernel density estimates
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))
#Kernel density estimation
for (i in 1:rep){
  Hpi <- ks::Hpi(x = data[,,i])
  kdeHpi <- ks::kde(x = data[,,i], H = Hpi, xmin=c(xmin,ymin), xmax=c(xmax,ymax), bgridsize=c(nx,ny))
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

#Function to calculate bivariate normal density
bivariate_normal_density = function(x, y, mu = c(0,0), Sigma = matrix(data = c(1, 0, 0, 1),
                                                                      nrow = 2, ncol = 2)){
  sigma1 = Sigma[1,1]^0.5
  sigma2 = Sigma[2,2]^0.5
  cov = Sigma[1,2]
  corr = cov/(sigma1 * sigma2)
  mu1 = mu[1]
  mu2 = mu[2]
  z = (x - mu1)^2/sigma1^2 - 2 * corr * (x - mu1) * (y - mu2) / (sigma1 * sigma2) +
    (y - mu2)^2 / sigma2^2
  f = 1/(2 * pi * sigma1 * sigma2 * (1 - corr^2)^0.5) * exp(-z/(2 * (1 - corr^2)))
  return(f)
}

#True density values
for (i in 1:nx) {
  xvalue = xcoordinates[i]
  for (j in 1:ny){
    yvalue = ycoordinates[j]
    true_values[i,j] = bivariate_normal_density(xvalue, yvalue, mu, Sigma = sigma)
  }
}

#Difference of true values and estimates
diff_est_true = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
        dimnames = list(as.character(xcoordinates), as.character(ycoordinates), as.character(c(1:rep))))
for (i in 1:rep){
  diff_est_true[,,i] = estimates[,,i] - true_values[,]
}

#Equalize average difference and bias
bias_estimates = aver_est - true_values

#Relative bias
rel_bias_estimates = bias_estimates / true_values

#Variances of estimates
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
