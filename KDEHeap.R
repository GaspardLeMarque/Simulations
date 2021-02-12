#KDE with Kernelheaping library

library(ggplot2) 
library(dplyr)
library(Kernelheaping)
library(bivariate)

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
