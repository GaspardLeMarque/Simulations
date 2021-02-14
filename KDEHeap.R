#KDE with Kernelheaping library

library(ggplot2) 
library(dplyr)
library(Kernelheaping)
library(bivariate)
library(sparr)
library(tidyverse)
library(ggExtra)

set.seed(123) 

#Simulate 5000 observations from a log-normal distribution
dat <- data.frame(x = rlnorm(5000, meanlog = 7.5, sdlog = 0.5)) 
ggplot(dat, aes(x, y = ..density..)) +  
  geom_histogram(position = "identity", color = "black", fill = "white") +
  geom_density(fill = "blue", alpha = 0.7, linetype = 0) +
  ylim(c(0, 0.00055)) +
  xlim(c(0, 10000))

#Split the data into 14 classes, imitating a typical income survey
classes <- c(0, 500, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000,
             10000, 15000, Inf) 
dat$xclass <- cut(dat$x, breaks = classes)

#Use class centers as input density
classCenters <- (classes[-1] - classes[-length(classes)]) /
  2 + classes[-length(classes)] 
classCenters[length(classCenters)] <- 2 * classCenters[length(classCenters) - 1] 
levels(dat$xclass) <- classCenters 
dat$xclassCenter <- as.numeric(as.character(dat$xclass)) 

#This plot is spiky and unrealistic
ggplot(dat, aes(xclassCenter, y = ..density..)) +
  geom_histogram(breaks = classes, color = "black", fill = "white") +
  geom_density(fill = "blue", alpha = 0.7, linetype = 0)

#Adjusting the density by changing the parameter (try to avoid oversmoothing)
ggplot(dat, aes(xclassCenter, y = ..density..)) +
  geom_histogram(breaks = classes, color = "black", fill = "white") +
  geom_density(adjust = 2.5, fill = "blue", alpha = 0.7, linetype = 0) +
  ylim(c(0, 0.00055)) +
  xlim(c(0, 10000))

###Implement an iterative Stochastic Expectation Maximization###

#Any bandwidth selector is possible to use here
#Use KDE for classified data
densityEst <- dclass(xclass = dat$xclass, classes = classes, burnin = 50,
                     samples = 100, evalpoints = 1000)
classes <- data.frame(xclass = densityEst$xclass) 
estimates <- data.frame(Mestimates = densityEst$Mestimates, 
                        gridx = densityEst$gridx) 

#Plot is smooth and close to the estimate on the original unclassified data 
ggplot() +
  geom_histogram(data = classes, 
                 aes(xclass, y = ..density..),
                 breaks = densityEst$classes, 
                 color = "black", fill = "white") +
  geom_area(data = estimates, 
            aes(gridx, Mestimates),
            fill = "blue", 
            alpha = 0.7, 
            size = 1) +
  ylim(c(0, 0.00055)) +
  xlim(c(0, 10000))

#Do Adaptive Kernel Estimate for bimodal bivariate normal mixture
#Bimodal bivariate normal mixture
bi_mix <- bmbvpdf(mean.X1 = 0, mean.Y1 = -3, sd.X1 = 1, sd.Y1 = 1.5, mean.X2 = -3, mean.Y2 = 3, sd.X2 = 1.2, sd.Y2 = 2.1)

#Acceptance-Rejection sampling
#Bivariate uniform with maximum of target as an upper function
x1_low <- -5
x1_up <- 5
x2_low <- -5
x2_up <- 5
upper <- cubvpdf(a.X = x1_low, b.X = x1_up, a.Y = x2_low, b.Y = x2_up) # interval

#Add a constant to the uniform distro
sampled <- data.frame(x = runif(100000, x1_low, x1_up), y = runif(100000,x2_low,x2_up))
sampled$density <- bi_mix(x = sampled$x, y = sampled$y)

#Increase the maximum by 10% to get higher value for sure
max_value <- max(sampled$density) 
max_value <- max_value * 1.1 

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

#Create a sample using bimodal bivariate distro as a target
#and upper function as an approximation
sample <- ARS(target = bi_mix, approx = upper, 
              x_low = x1_low, x_up = x1_up, 
              y_low = x2_low, y_up = x2_up, 
              n = 5000, c = max_value)

#First attempt to use Adaptive Kernel
new_sample <- as.ppp(sample, c(-5, 5, -5, 5))
est <- bivariate.density(new_sample,h0=1.5,hp=1,adapt=TRUE)
plot(est)
