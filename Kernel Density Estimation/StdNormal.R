#KDE in the multivariate case
library(ks) #KDE
library(MASS)
library(plotly)
library(bivariate) #produces the density functions
library(tidyverse)
library(tictoc) #measures performance
library(gsl)
library(copula)
library(reshape2)
library(bivariate)
library(sparr) #Adaptive Kernel Density Estimate
library(TDA) #kNN Density



### Functions ###
#Silverman's rule of thumb formula
calc_sil <- function(data){
  len <- length(data)/2
  sd_1 <- sd(data[,1])
  sd_2 <- sd(data[,2])
  H_11 <- len^(-1/(6)) * sd_1
  H_22 <- len^(-1/(6)) * sd_2
  H <- diag(x=c(H_11, H_22))
  return(H)
}

#Sampling function (Acception-Rejection method)
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

#Function for Riemann integral approximation
apxrint <- function(x, y = NULL, f) {
  stopifnot(is.matrix(f)|is.vector(f))
  stopifnot(is.vector(x))
  stopifnot(is.vector(y)| is.null(y))
  
  if (is.vector(f) & is.vector(x)){
    stopifnot(length(f) == length(x))
    nonan <- !is.na(f)
    nn <- sum(nonan)
    if(nn < 2L) return(0)
    Y <- f[nonan]
    X <- x[nonan]
    return(0.5 * sum( (Y[-1] + Y[-nn]) * diff(X)))
  }
  else if (is.matrix(f) & is.vector(x) & is.vector(y)){
    Lf <- dim(f); Lx <- length(x); Ly <- length(y)
    stopifnot((Lf == c(Lx, Ly)) | (t(Lf) == c(Lx, Ly)))
    nan <- is.na(f)
    f[nan] <- 0
    f[2:(Lx-1), ] = f[2:(Lx-1), ] * 2;
    f[, 2:(Ly-1)] = f[ , 2:(Ly-1)] * 2;
    return(sum(f) * diff(x)[1] * diff(y)[1] / 4)
  }
}

#### Standard Normal #####

#Study preparation
set.seed(1)

#Number of draws per sample
n = 100

#Number of repetitions in simulation
rep = 1000

#Number of dimensions of data
dim = 2

#Grid resolution
nx = 128
ny = nx
xmin = -5
ymin = -5
xmax = 5
ymax = 5
xcoordinates = seq(xmin, xmax, length.out = nx)
ycoordinates = seq(ymin, ymax, length.out = ny)
Grid <- expand.grid(xcoordinates, ycoordinates)

#Create an empty data array with the correct dimensions
#1. dimension: draws, 2. dimension: variables, 3. dimension: simulations
data = array(data = rep(NA, n * dim * rep), c(n, dim, rep))

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

#Data generation
mu = c(0,0)
sigma = matrix(data = c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
for (i in 1:rep){
  sample = MASS::mvrnorm(n = n, mu = mu, Sigma = sigma)
  data[,,i] = sample
}

#Theoretical density values on grid points
true_values = array(data = rep(NA, nx * ny), c(nx, ny),
                    dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

#True density values
for (i in 1:nx) {
  xvalue = xcoordinates[i]
  for (j in 1:ny){
    yvalue = ycoordinates[j]
    true_values[i,j] = bivariate_normal_density(xvalue, yvalue, mu, Sigma = sigma)
  }
}


#Plug-in KDE and Silverman's Rule

#Choose selector
selector = "sil" #OR selector = "default"

#Save the Kernel density estimates
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))

#Kernel density estimation
tic.clearlog() # to benchmark
if (selector == "sil"){#Silverman's rule of thumb case
  for (i in 1:rep){
    tic()
    H <- calc_sil(data[,,i])
    kdeHpi <- ks::kde(x = data[,,i], H = H, xmin=c(xmin,ymin), xmax=c(xmax,ymax), gridsize=c(nx,ny))
    estimates[,,i] = kdeHpi$estimate
    toc(log = TRUE, quiet = TRUE)
  }
} else { #Standard case
  for (i in 1:rep){
    tic()
    Hpi <- ks::Hpi(x = data[,,i])
    kdeHpi <- ks::kde(x = data[,,i], H = Hpi, xmin=c(xmin,ymin), xmax=c(xmax,ymax), gridsize=c(nx,ny))
    estimates[,,i] = kdeHpi$estimate
    toc(log = TRUE, quiet = TRUE)
  }
}
log.txt_normal_silver <- tic.log(format = TRUE)
# log.txt_normal <- tic.log(format = TRUE) when using default selector

#Average estimate
aver_est = array(data = rep(NA, nx * ny), c(nx, ny),
                 dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx) {
  for (j in 1:ny){
    aver_est[i,j] = mean(estimates[i,j,])
  }
}

#Difference of true values and estimates
diff_est_true = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                      dimnames = list(as.character(xcoordinates), as.character(ycoordinates), as.character(c(1:rep))))
for (i in 1:rep){
  diff_est_true[,,i] = estimates[,,i] - true_values[,]
}

#Average difference = bias
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

#Squared differences of estimate and true value
sq_diff_est_true = diff_est_true^2

if(selector == 'sil'){
  IMSE_silver_standard_normal = vector(mode = 'numeric', length = rep)
  for (i in 1:rep){
    IMSE_silver_standard_normal[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
  }
} else{
  IMSE_standard_normal = vector(mode = 'numeric', length = rep)
  for (i in 1:rep){
    IMSE_standard_normal[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
  }
}

#Save the data (Change for different selectors)
aver_est_normal_silver <- aver_est
#OR aver_est_normal <- aver_est
bias_estimates_normal_silver <- bias_estimates 
#OR bias_estimates_normal <- bias_estimates 
rmse_estimates_normal_silver <- rmse_estimates
#OR rmse_estimates_normal <- rmse_estimates
true_est_normal_silver <- true_values
#OR true_est_normal <- true_values
IMSE_normal_silver <- IMSE_silver_standard_normal
#OR IMSE_normal <- IMSE_standard_normal
save(list=c("IMSE_normal_silver", 
            "aver_est_normal_silver", 
            "bias_estimates_normal_silver", 
            "rmse_estimates_normal_silver", 
            "true_est_normal_silver", 
            "log.txt_normal_silver"), file = "normal_silver.RDA")
#OR save(list=c("IMSE_normal", "aver_est_normal", "bias_estimates_normal", "rmse_estimates_normal", "true_est_normal", "log.txt_normal"), file = "normal.RDA")

### Adaptive KDE ### 
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))
#AKDE
h0 <- 0.20 
tic.clearlog()
for (i in 1:rep) {
  tic()
  ppp = as.ppp(X = data[,,i], c(xmin, xmax, ymin, ymax))
  pilot = bivariate.density(pp = ppp, h0=h0, parallelise = TRUE, resolution = nx)
  AKDE = bivariate.density(pp = ppp, h0=h0, pilot.density = pilot$z, adapt = TRUE, parallelise = TRUE, resolution = nx, verbose = FALSE)
  estimates[,,i] = t(AKDE$z$v)
  toc(log = TRUE, quiet = TRUE)
}

log.txt_normal_AKDE <- tic.log(format = TRUE)

#Average estimate
aver_est = array(data = rep(NA, nx * ny), c(nx, ny),
                 dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx) {
  for (j in 1:ny){
    aver_est[i,j] = mean(estimates[i,j,])
  }
}

#Difference of true values and estimates
diff_est_true = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                      dimnames = list(as.character(xcoordinates), as.character(ycoordinates), as.character(c(1:rep))))
for (i in 1:rep){
  diff_est_true[,,i] = estimates[,,i] - true_values[,]
}

#Average difference = bias
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

#Squared differences of estimate and true value
sq_diff_est_true = diff_est_true^2
IMSE = vector(mode = 'numeric', length = rep)
for (i in 1:rep){
  IMSE[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
}

#Save the data
aver_est_normal_AKDE <- aver_est
IMSE_normal_AKDE <- IMSE
rmse_estimates_normal_AKDE <- rmse_estimates
bias_estimates_normal_AKDE <- bias_estimates
save(list=c("aver_est_normal_AKDE", 
            "IMSE_normal_AKDE", 
            "rmse_estimates_normal_AKDE", 
            "bias_estimates_normal_AKDE", 
            "log.txt_normal_AKDE"), file = "standard_normal_AKDE.RDA")

### kNN Density ###

estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))

tic.clearlog()
for(i in 1:rep){
  tic()
  KNN <- knnDE(X=data[,,i], Grid = Grid, k = 10) #k selected with k=sqrt(n)
  estimates[,,i] <- matrix(data = KNN, nrow = length(xcoordinates))
  toc(log = TRUE, quiet = TRUE)
}
log.txt_standard_normal_KNN <- tic.log(format = TRUE)

#Average estimate
aver_est = array(data = rep(NA, nx * ny), c(nx, ny),
                 dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx) {
  for (j in 1:ny){
    aver_est[i,j] = mean(estimates[i,j,])
  }
}

#Difference of true values and estimates
diff_est_true = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                      dimnames = list(as.character(xcoordinates), as.character(ycoordinates), as.character(c(1:rep))))
for (i in 1:rep){
  diff_est_true[,,i] = estimates[,,i] - true_values[,]
}

#Average difference = bias
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
sq_diff_est_true = diff_est_true^2
IMSE = vector(mode = 'numeric', length = rep)
for (i in 1:rep){
  IMSE[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
}

#Save the data
IMSE_standard_normal_KNN <- IMSE
aver_est_standard_normal_KNN <- aver_est
rmse_estimates_standard_normal_KNN <- rmse_estimates
bias_estimates_standard_normal_KNN <- bias_estimates

save(list=c("aver_est_standard_normal_KNN", 
            "IMSE_standard_normal_KNN", 
            "rmse_estimates_standard_normal_KNN", 
            "bias_estimates_standard_normal_KNN", 
            "log.txt_standard_normal_KNN"), file = "standard_normal_KNN.RDA")
