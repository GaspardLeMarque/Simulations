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

##### Gamma Beta Distribution #####
#Study preparation
set.seed(1)

#Number of draws per sample
n = 100

#Number of repetitions in simulation
rep = 1000

#Number of dimensions of data
dim = 2

#Grid resolution for estimates
nx = 128
ny = nx
xmin = 0
ymin = 0
xmax = 10
ymax = 1
xcoordinates = seq(xmin, xmax, length.out = nx)
ycoordinates = seq(ymin, ymax, length.out = ny)
Grid <- expand.grid(xcoordinates, ycoordinates)

#Create an empty data array with the correct dimensions
#1. dimension: draws, 2. dimension: variables, 3. dimension: simulations
data = array(data = rep(NA, n * dim * rep), c(n, dim, rep))

#Combining gamma and beta distributions
gamma_beta <- normalCopula(param = c(0.4), dim = 2)
gamma_beta <- mvdc(copula = gamma_beta, 
                   margins = c("gamma", "beta"), 
                   paramMargins = list(list(shape=2, scale=1), list(shape1=2, shape2=2)))

#Data generation
for (i in 1:rep){
  sample <- rMvdc(n, gamma_beta)
  data[,1,i] = sample[,1]
  data[,2,i] = sample[,2]
}

#True density values on specified grid points
true_values = array(data = rep(NA, nx * ny), c(nx, ny),
                    dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx) {
  xvalue = xcoordinates[i]
  points = cbind(c(rep(xvalue, nx)), ycoordinates)
  true_values[i,] = dMvdc(x = points, gamma_beta)
}

#Plug-in KDE and Silverman's Rule

#Choose selector
selector = "sil" #OR selector = "default"

#Save the Kernel density estimates
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))

#Kernel density estimation
tic.clearlog() #to benchmark
if (selector == "sil"){ #Silverman's rule of thumb case
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

log.txt_gamma_beta_sil <- tic.log(format = TRUE)
#log.txt_gamma_beta <- tic.log(format = TRUE)

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
  IMSE_silver_beta_gamma = vector(mode = 'numeric', length = rep)
  for (i in 1:rep){
    IMSE_silver_beta_gamma[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
  }
} else{
  IMSE_beta_gamma = vector(mode = 'numeric', length = rep)
  for (i in 1:rep){
    IMSE_beta_gamma[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
  }
}

#Save the data (Change for different selectors)
aver_est_gamma_beta_silver <- aver_est
#aver_est_gamma_beta <- aver_est
bias_estimates_gamma_beta_silver <- bias_estimates 
#bias_estimates_gamma_beta <- bias_estimates 
rmse_estimates_gamma_beta_silver <- rmse_estimates
#rmse_estimates_gamma_beta <- rmse_estimates
true_est_gamma_beta_silver <- true_values
#true_est_gamma_beta <- true_values
IMSE_gamma_beta_silver <- IMSE_silver_beta_gamma
#IMSE_gamma_beta <- IMSE_beta_gamma
save(list=c("IMSE_gamma_beta_silver", 
			"aver_est_gamma_beta_silver", 
			"bias_estimates_gamma_beta_silver", 
			"rmse_estimates_gamma_beta_silver", 
			"true_est_gamma_beta_silver", 
			"log.txt_gamma_beta_sil"), file = "gamma_beta_silver.RDA")
#save(list=c("IMSE_gamma_beta", "aver_est_gamma_beta", "bias_estimates_gamma_beta", "rmse_estimates_gamma_beta", "true_est_gamma_beta", "log.txt_gamma_beta_sil"), file = "gamma_beta.RDA")

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

log.txt_gamma_beta_AKDE <- tic.log(format = TRUE)

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
aver_est_gamma_beta <- aver_est
IMSE_gamma_beta <- IMSE
rmse_estimates_gamma_beta <- rmse_estimates
bias_estimates_gamma_beta <- bias_estimates
log.txt_gamma_beta <- log.txt

save(list=c("aver_est_gamma_beta", 
			"IMSE_gamma_beta", 
			"rmse_estimates_gamma_beta", 
			"bias_estimates_gamma_beta", 
			"log.txt_gamma_beta"), file = "beta_gamma_akde.RDA")

### kNN Density ###

estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))

#Kernel density estimation
tic.clearlog()
for(i in 1:rep){
  tic()
  KNN <- knnDE(X=data[,,i], Grid = Grid, k = 10)
  estimates[,,i] <- matrix(data = KNN, nrow = length(xcoordinates))
  toc(log = TRUE, quiet = TRUE)
}
log.txt_gamma_beta_KNN <- tic.log(format = TRUE)

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

IMSE_gamma_beta_KNN <- IMSE
aver_est_gamma_beta_KNN <- aver_est
rmse_estimates_gamma_beta_KNN <- rmse_estimates
bias_estimates_gamma_beta_KNN <- bias_estimates

save(list=c("aver_est_gamma_beta_KNN", 
			"IMSE_gamma_beta_KNN", 
			"rmse_estimates_gamma_beta_KNN", 
			"bias_estimates_gamma_beta_KNN", 
			"log.txt_gamma_beta_KNN"), file = "gamma_beta_KNN.RDA")
