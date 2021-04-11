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

##### Bimodal Normal #####
# Study preparation
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
xmin = -10
ymin = -10
xmax = 10
ymax = 10
xcoordinates = seq(xmin, xmax, length.out = nx)
ycoordinates = seq(ymin, ymax, length.out = ny)
Grid <- expand.grid(xcoordinates, ycoordinates)

#Creature mixture distribution (package "bivariate")
bi_mix <- bmbvpdf(mean.X1 = 0, mean.Y1 = -3, sd.X1 = 1, sd.Y1 = 1.5, mean.X2 = -3, mean.Y2 = 3, sd.X2 = 1.2, sd.Y2 = 2.1)

#Upper function. Bivariate uniform with maximum of target. We find the constant to add via simulation
x1_low <- -10
x1_up <- 10
x2_low <- -10
x2_up <- 10
upper <- cubvpdf(a.X = x1_low, b.X = x1_up, a.Y = x2_low, b.Y = x2_up) # interval

#Add a constant 
sampled <- data.frame(x = runif(100000, x1_low, x1_up), y = runif(100000,x2_low,x2_up))
sampled$density <- bi_mix(x = sampled$x, y = sampled$y)

max_value <- max(sampled$density) # we take this or a slightly higher value as constant to add to the uniform
max_value <- max_value * 1.1 # increase by 10% to make sure, we actually are higher then the peak

rm(sampled)

data = array(data = rep(NA, n * dim * rep), c(n, dim, rep))

#Data generation
for (i in 1:rep){
  sample <- ARS(target = bi_mix, approx = upper, x_low = x1_low,
                x_up = x1_up, y_low = x2_low, y_up = x2_up, n = n, c = max_value)
  data[,1,i] = sample$x
  data[,2,i] = sample$y
}

#Calculate theoretical density values on grid points
true_values = array(data = rep(NA, nx * ny), c(nx, ny),
                    dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

#True density values
for (i in 1:nx) {
  xvalue = xcoordinates[i]
  for (j in 1:ny){
    yvalue = ycoordinates[j]
    true_values[i,j] = bi_mix(xvalue, yvalue)
  }
}

## Plug-in KDE and Silvermans Rule

#Choose selector
selector <- "sil"
#selector <- "default"

#Save the Kernel density estimates
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))
#Kernel density estimation
tic.clearlog()
if (selector == "sil"){
  for (i in 1:rep){
    tic()
    H <- calc_sil(data[,,i])
    kdeHpi <- ks::kde(x = data[,,i], H = H, xmin=c(xmin,ymin), xmax=c(xmax,ymax), gridsize=c(nx,ny))
    estimates[,,i] = kdeHpi$estimate
    toc(log = TRUE, quiet = TRUE)
  }
} else { 
  for (i in 1:rep){
    tic()
    Hpi <- ks::Hpi(x = data[,,i])
    kdeHpi <- ks::kde(x = data[,,i], H = Hpi, xmin=c(xmin,ymin), xmax=c(xmax,ymax), gridsize=c(nx,ny))
    estimates[,,i] = kdeHpi$estimate
    toc(log = TRUE, quiet = TRUE)
  }
}
log.txt_bimodal_silver <- tic.log(format = TRUE)
#log.txt_bimodal <- tic.log(format = TRUE)

#Calculate average estimate
aver_est = array(data = rep(NA, nx * ny), c(nx, ny),
                 dimnames = list(as.character(xcoordinates), as.character(ycoordinates)))

for (i in 1:nx) {
  for (j in 1:ny){
    aver_est[i,j] = mean(estimates[i,j,])
  }
}

#Difference of true values and estimates
diff_est_true = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                      dimnames = list(as.character(xcoordinates),
                                      as.character(ycoordinates), as.character(c(1:rep))))
for (i in 1:rep){
  diff_est_true[,,i] = estimates[,,i] - true_values[,]
}

#Average difference = bias
bias_estimates = aver_est - true_values

#Relative bias
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

#Squared differences of estimate and true value
sq_diff_est_true = diff_est_true^2

if(selector == 'sil'){
  IMSE_silver_bimodal_normal = vector(mode = 'numeric', length = rep)
  for (i in 1:rep){
    IMSE_silver_bimodal_normal[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
  }
} else{
  IMSE_bimodal_normal = vector(mode = 'numeric', length = rep)
  for (i in 1:rep){
    IMSE_bimodal_normal[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
  }
}

#Save the data (Change for different selectors)
aver_est_bimodal_normal_silver <- aver_est
bias_estimates_bimodal_normal_silver <- bias_estimates 
rmse_estimates_bimodal_normal_silver <- rmse_estimates
true_est_bimodal_normal_silver <- true_values
IMSE_bimodal_normal_silver <- IMSE_silver_bimodal_normal
save(list=c("IMSE_bimodal_normal_silver", 
            "aver_est_bimodal_normal_silver", 
            "bias_estimates_bimodal_normal_silver", 
            "rmse_estimates_bimodal_normal_silver", 
            "true_est_bimodal_normal_silver", 
            "log.txt_bimodal_silver"), file = "bimodal_silver.RDA")

#aver_est_bimodal_normal <- aver_est
#bias_estimates_bimodal_normal <- bias_estimates 
#rmse_estimates_bimodal_normal <- rmse_estimates
#true_est_bimodal_normal <- true_values
#IMSE_bimodal_normal <- IMSE_bimodal_normal

#OR save(list=c("IMSE_bimodal_normal", "aver_est_bimodal_normal", "bias_estimates_bimodal_normal", "rmse_estimates_bimodal_normal", "true_est_bimodal_normal", "log.txt_bimodal"), file = "bimodal.RDA")

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

log.txt_bimodal_AKDE <- tic.log(format = TRUE)

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

#Squared differences of estimate and true
sq_diff_est_true = diff_est_true^2
IMSE = vector(mode = 'numeric', length = rep)
for (i in 1:rep){
  IMSE[i] <- apxrint(xcoordinates, ycoordinates, sq_diff_est_true[,,i])
}

aver_est_bimodal_AKDE <- aver_est
IMSE_bimodal_AKDE <- IMSE
rmse_estimates_bimodal_AKDE <- rmse_estimates
bias_estimates_bimodal_AKDE <- bias_estimates

save(list=c("aver_est_bimodal_AKDE", 
			"IMSE_bimodal_AKDE", 
			"rmse_estimates_bimodal_AKDE", 
			"bias_estimates_bimodal_AKDE", 
			"log.txt_bimodal_AKDE"), file = "bimodal_AKDE.RDA")

### kNN Density ###

data = array(data = rep(NA, n * dim * rep), c(n, dim, rep))
estimates = array(data = rep(NA, nx * ny * rep), c(nx, ny, rep),
                  dimnames = list(as.character(xcoordinates), as.character(ycoordinates),
                                  as.character(c(1:rep))))

tic.clearlog()
for(i in 1:rep){
  tic()
  KNN <- knnDE(X=data[,,i], Grid = Grid, k = 10)
  estimates[,,i] <- matrix(data = KNN, nrow = length(xcoordinates))
  toc(log = TRUE, quiet = TRUE)
}

log.txt_bimodal_KNN <- tic.log(format = TRUE)

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

#True density values
for (i in 1:nx) {
  xvalue = xcoordinates[i]
  for (j in 1:ny){
    yvalue = ycoordinates[j]
    true_values[i,j] = bi_mix(xvalue, yvalue)
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
IMSE_bimodal_KNN <- IMSE
aver_est_bimodal_KNN <- aver_est
rmse_estimates_bimodal_KNN <- rmse_estimates
bias_estimates_bimodal_KNN <- bias_estimates

save(list=c("aver_est_bimodal_KNN", 
			"IMSE_bimodal_KNN", 
			"rmse_estimates_bimodal_KNN", 
			"bias_estimates_bimodal_KNN", 
			"log.txt_bimodal_KNN"), file = "bimodal_KNN.RDA")
